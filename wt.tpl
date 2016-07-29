DATA_SECTION
  init_int styr
  init_int endyr
  init_int age_st
  init_int age_end
  int nages;
  !! nages = age_end - age_st +1;
  init_matrix wt_obs(styr,endyr,age_st,age_end)
  init_matrix sd_obs(styr,endyr,age_st,age_end)
	int ret
 LOCAL_CALCS
  ret=0;
  if (ad_comm::argc > 1)
  {
    int on=0;
    if ( (on=option_match(ad_comm::argc,ad_comm::argv,"-ret"))>-1)
    {
      if (on>ad_comm::argc-2 | ad_comm::argv[on+1][0] == '-')
      {
        cerr << "Invalid number of iseed arguements, command line option -ret ignored" << endl;
      }
      else
      {
        ret = atoi(ad_comm::argv[on+1]);
      }
    }
		endyr= endyr -ret;
  }
  cout <<endyr<<endl;//exit(1);
 END_CALCS
INITIALIZATION_SECTION

PARAMETER_SECTION
  // Predicted weight matrix
  matrix wt_pre(styr,endyr+3,age_st,age_end);
  init_bounded_number log_sd_coh(-5,10,3);
  init_bounded_number log_sd_yr(-5,10,4);
  init_vector mnwt(age_st,age_end);
  vector wt_inc(age_st+1,age_end);
  sdreport_number sig_coh;
  sdreport_number sig_yr;
  sdreport_vector wt_last(age_st,age_end);
  sdreport_vector wt_cur(age_st,age_end);
  sdreport_vector wt_next(age_st,age_end);
  sdreport_vector wt_yraf(age_st,age_end);


  // init_bounded_vector coh_eff(styr-nages-age_st+1,endyr-age_st,-5,5,2);
  // init_bounded_vector yr_eff(styr,endyr,-5,5,3);
  random_effects_vector coh_eff(styr-nages-age_st+1,endyr-age_st+3,3);
  random_effects_vector yr_eff(styr,endyr+3,3);

  objective_function_value nll;

PROCEDURE_SECTION
  dvariable sigma_coh = (mfexp(log_sd_coh));
  dvariable sigma_yr  = (mfexp(log_sd_yr ));
  wt_inc = mnwt(age_st+1,age_end) - ++mnwt(age_st,age_end-1);
  // cout <<wt_inc <<endl;exit(1);
  // Initialize first year
  wt_pre(styr) = mnwt;
  wt_pre(styr)(age_st+1,age_end) = ++wt_pre(styr)(age_st,age_end-1) + wt_inc*exp(sigma_yr*yr_eff(styr));
  for (int j=age_st;j<=age_end;j++)
  {
    wt_pre(styr,j) *= exp(sigma_coh*coh_eff(styr-j));
    nll += square(wt_obs(styr,j)-wt_pre(styr,j))/(2.*square(sd_obs(styr,j)));
  }

  // subsequent years
  for (int i=styr+1;i<=endyr+3;i++)
  {
    wt_pre(i)(age_st+1,age_end) = ++wt_pre(i-1)(age_st,age_end-1) + wt_inc*exp(sigma_yr*yr_eff(i));
    // cout <<i<<" "<<wt_pre(i) <<endl;// exit(1);
    for (int j=age_st;j<=age_end;j++)
    {
      wt_pre(i,j) *= exp(sigma_coh*coh_eff(i-j));
      if (i <= endyr)
        nll += square(wt_obs(i,j)-wt_pre(i,j))/(2.*square(sd_obs(i,j)));
    }
  }
  // if (current_phase()>2)
  {
		nll += 0.5*norm2(coh_eff);
		nll += 0.5*norm2( yr_eff);
    // nll +=  yr_eff.size()*log_sd_yr  + norm2(yr_eff )/(2*sigma_yr );
  }
  if (sd_phase())
  {
    sig_coh = exp(log_sd_coh);
    sig_yr  = exp(log_sd_yr );
    wt_last = wt_pre(endyr)  *exp(sigma_coh*sigma_coh/2. );// + sigma_yr*sigma_yr/2.);
    wt_cur  = wt_pre(endyr+1)*exp(sigma_coh*sigma_coh/2. );// + sigma_yr*sigma_yr/2.);;
    wt_next = wt_pre(endyr+2)*exp(sigma_coh*sigma_coh/2. );// + sigma_yr*sigma_yr/2.);;
    wt_yraf = wt_pre(endyr+3)*exp(sigma_coh*sigma_coh/2. );// + sigma_yr*sigma_yr/2.);;
  }

REPORT_SECTION
 report << wt_pre<<endl;
