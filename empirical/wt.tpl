DATA_SECTION
  init_int retyr ;// Retrospective years to test over (0 if current year)
  init_int styr;
  init_int lyr;
  init_int sage;
  init_int lage;
  !! cout<< styr <<" "<< lyr <<" "<< sage <<" "<< lage<<endl;
  int lyr_r; 
  int nages;
  !! nages = lage - sage +1;
  init_matrix wt_obs(styr,lyr,sage,lage)
  init_matrix sd_obs(styr,lyr,sage,lage)
  vector ages(sage,lage)
 LOCAL_CALCS
  for (int i=sage;i<=lage;i++) ages(i) = double(i);
  retyr=0;
  if (ad_comm::argc > 1)
  {
    int on=0;
    if ( (on=option_match(ad_comm::argc,ad_comm::argv,"-ret"))>-1) {
      if (on>ad_comm::argc-2 | ad_comm::argv[on+1][0] == '-') {
        cerr << "Invalid number of iseed arguements, command line option -ret ignored" << endl;
      }
      else {
        retyr = atoi(ad_comm::argv[on+1]);
      }
    }
		lyr_r= lyr - retyr;
  }
 END_CALCS
INITIALIZATION_SECTION
  log_winf 0.75
  log_K -2.
  log_t0 .4
  log_sd_coh -1.
  log_sd_yr -.77

PARAMETER_SECTION
  // Predicted weight matrix
  init_number log_winf(1);
  init_number log_K(1);
  init_number log_t0(1);
  init_bounded_number log_sd_coh(-5,10,4);
  init_bounded_number log_sd_yr(-5,10,4);
  !! cout<< styr <<" "<< lyr <<" "<< sage <<" "<< lage<<endl;
  matrix wt_pre(styr,lyr_r+3,sage,lage);
  !! cout<< wt_pre.indexmin() <<" "<< wt_pre.indexmax() <<" "<< sage <<" "<< lage<<endl;
  vector mnwt(sage,lage);
  vector wt_inc(sage,lage-1);
  sdreport_number sigma_coh;
  sdreport_number sigma_yr;
  sdreport_number winf;
  sdreport_number K;
  sdreport_number t0;
  sdreport_vector wt_last(sage,lage);
  sdreport_vector wt_cur(sage,lage);
  sdreport_vector wt_next(sage,lage);
  sdreport_vector wt_yraf(sage,lage);

  random_effects_vector coh_eff(styr-nages-sage+1,lyr_r-sage+3,2);
  random_effects_vector yr_eff(styr,lyr_r+3,3);
  // init_vector coh_eff(styr-nages-sage,lyr_r-sage+3,2);
  // init_vector yr_eff(styr,lyr_r+3,3);

  objective_function_value nll;

PROCEDURE_SECTION
  sigma_coh    = exp(log_sd_coh);
  sigma_yr     = exp(log_sd_yr );
  winf         = exp(log_winf);
  K            = exp(log_K);
  t0           = (log_t0);
  mnwt         = winf*(pow(1.-exp(-K*(ages-t0)),3));
  wt_inc       = --mnwt(sage+1,lage) - mnwt(sage,lage-1);

  // Initialize first year
  wt_pre(styr) = mnwt;
  for (int j=sage;j<=lage;j++)
  {
    if (j>sage) 
      wt_pre(styr,j) = wt_pre(styr,j-1) + wt_inc(j-1)*exp(sigma_yr*yr_eff(styr));
    int icoh        = styr-j;
    wt_pre(styr,j) *= mfexp(sigma_coh*coh_eff(icoh));
    // Fit global mean to all years...
    nll += square(wt_obs(styr,j)-mnwt(j))/(2.*square(sd_obs(styr,j)));
    // Fit annual predictions to all years...
    nll += square(wt_obs(styr,j)-wt_pre(styr,j))/(2.*square(sd_obs(styr,j)));
  }

  // subsequent years
  for (int i=styr+1;i<=lyr_r+3;i++)
  {
    wt_pre(i,sage)            = mnwt(sage)*exp(sigma_coh*coh_eff(i-sage));
    wt_pre(i)(sage+1,lage) = ++(wt_pre(i-1)(sage,lage-1) + wt_inc*exp(sigma_yr*yr_eff(i)));
    for (int j=sage;j<=lage;j++)
    {
      if (i <= lyr_r)
      {
        nll += square(wt_obs(i,j)-mnwt(j))/(2.*square(sd_obs(i,j)));
        nll += square(wt_obs(i,j) - wt_pre(i,j))/(2.*square(sd_obs(i,j)));
      }
    }
  }
	nll += 0.5*norm2(coh_eff);
	nll += 0.5*norm2( yr_eff);

  if (sd_phase())
  {
    sigma_coh = exp(log_sd_coh);
    sigma_yr  = exp(log_sd_yr );
    wt_last   = wt_pre(lyr_r)  ; // *exp(sigma_coh*sigma_coh/2. );// + sigma_yr*sigma_yr/2.);
    wt_cur    = wt_pre(lyr_r+1); // *exp(sigma_coh*sigma_coh/2. );// + sigma_yr*sigma_yr/2.);;
    wt_next   = wt_pre(lyr_r+2); // *exp(sigma_coh*sigma_coh/2. );// + sigma_yr*sigma_yr/2.);;
    wt_yraf   = wt_pre(lyr_r+3); // *exp(sigma_coh*sigma_coh/2. );// + sigma_yr*sigma_yr/2.);;
  }

REPORT_SECTION
  /* 
  switch ( retyr ) 
  {
    case 0 : 
      report << "lof1"<<endl; 0;
      report << "lof2"<<endl; 0;
      report << "lof3"<<endl; 0;
    break;

    case 1 : 
      report << "lof1"<<endl; pow(norm2(wt_obs(lyr_r+1)-wt_cur),.5);
      report << "lof2"<<endl; 0;
      report << "lof3"<<endl; 0;
    break;

    case 2 : 
      report << "lof1"<<endl; 0;
      report << "lof2"<<endl; 0;
      report << "lof3"<<endl; 0;
    break;

    default : 
    // Process for all other cases.
    break;
   }
  */ 

  report<<"data"<<endl;
  report<<wt_obs<<endl;
  report <<"W"<<endl;
  report << wt_pre<<endl;
  if (last_phase())
  {
    dmatrix resid(styr,lyr_r,sage,lage);
    for (int i=styr;i<=lyr_r;i++)
      resid(i) = value(wt_obs(i)-wt_pre(i));
    report <<"residuals"<<endl;
    report << resid << endl;
  }
  REPORT(sigma_yr);
  REPORT(yr_eff);
  REPORT(sigma_coh);
  report <<"cohort"<<endl;
  for (int j=lage;j>=sage;j--)
    report << styr-j <<" ";
  for (int i=styr+1;i<=lyr_r+3;i++)
    report << i-sage <<" ";
  report<<endl;
  report <<"coh_eff"<<endl;
  for (int j=lage;j>=sage;j--)
    report << sigma_coh*coh_eff(styr-j) <<" ";
  for (int i=styr+1;i<=lyr_r+3;i++)
    report << sigma_coh*coh_eff(i-sage) <<" ";
  report<<endl;
  REPORT(ages);
  REPORT(mnwt);
  REPORT(winf);
  REPORT(K);
  REPORT(t0);

FINAL_SECTION
  /* for (int i=lyr+1;i<=lyr+3;i++) 
  {
    cout << (W_sd(i)) <<endl;
    cout << (W_sd.sd(i))<<endl;
  }
  */ 

GLOBALS_SECTION
  #undef REPORT
  #define REPORT(object) report << #object "\n" << setw(8)  << setprecision(6) << setfixed() << object << endl;

  #undef COUT
  #define COUT(object) cout << #object "\n" << setw(6) << setprecision(5) << setfixed() << object << endl;

