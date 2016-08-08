DATA_SECTION
  init_int endyr_r ;// Retrospective years to test over (0 if current year)
  init_int styr
  init_int endyr
  init_int age_st
  init_int age_end
  int nages;
  init_matrix wt_obs(styr,endyr,age_st,age_end)
  init_matrix sd_obs(styr,endyr,age_st,age_end)
  vector ages(age_st,age_end)
  init age_st_s
  init age_end_s
  !! nages = age_end - age_st +1;
  !! for (int i=age_st;i<=age_end;i++) ages(i) = double(i);
  !! cout<< styr <<" "<< endyr <<" "<< age_st <<" "<< age_end<<endl;
  
INITIALIZATION_SECTION
  L1 30
  L2 55
  log_K -.3
  log_alpha -11.
  // log_t0 -3.4
  log_sd_coh -2.
  log_sd_yr -1.2

PARAMETER_SECTION
  init_bounded_number L1(10,50,1);
  init_bounded_number L2(30,90,2);
  init_number log_alpha(-1);
  init_number log_K(2);
  // init_number log_t0(3);
  init_bounded_number log_sd_coh(-5,10,5);
  init_bounded_number log_sd_yr(-5,10,4);
  // Predicted weight matrix
  matrix wt_pre(styr,endyr_r+3,age_st,age_end);
  number alpha;
  vector mnwt(age_st,age_end);
  vector wt_inc(age_st,age_end-1);
  sdreport_number sig_coh;
  sdreport_number sig_yr;
  sdreport_number sigma_coh;
  sdreport_number sigma_yr;
  // sdreport_number winf;
  sdreport_number K;
  // sdreport_number t0;
  sdreport_vector wt_last(age_st,age_end);
  sdreport_vector wt_cur(age_st,age_end);
  sdreport_vector wt_next(age_st,age_end);
  sdreport_vector wt_yraf(age_st,age_end);


  random_effects_vector coh_eff(styr-nages-age_st+1,endyr_r-age_st+3,3);
  random_effects_vector yr_eff(styr,endyr_r+3,4);

  objective_function_value nll;

PROCEDURE_SECTION
  sigma_coh    = exp(log_sd_coh);
  sigma_yr     = exp(log_sd_yr );
  // winf         = exp(log_winf);
  K            = exp(log_K);
  alpha        = exp(log_alpha);
  // t0           = (log_t0);
  // mnwt         = winf*(pow(1.-exp(-K*(ages-t0)),3));
  for (int j=age_st;j<=age_end;j++)
  {
    mnwt(j)    = alpha * pow(L1 + (L2-L1)*(1.-pow(K,double(j-age_st))) / (1.-pow(K,double(nages-1))) ,3);
    // cout << j << " " << mnwt(j)<<endl;
  }
  wt_inc       = --mnwt(age_st+1,age_end) - mnwt(age_st,age_end-1);

  // Initialize first year
  wt_pre(styr)    = mnwt;
  for (int j=age_st;j<=age_end;j++)
  {
    if (j>age_st) 
      wt_pre(styr,j) = wt_pre(styr,j-1) + wt_inc(j-1)*exp(sigma_yr*yr_eff(styr));
    int icoh        = styr-j;
    wt_pre(styr,j) *= mfexp(sigma_coh*coh_eff(icoh));
    // Fit global mean to all years...
    nll += square(wt_obs(styr,j)-mnwt(j))/(2.*square(sd_obs(styr,j)));
    // Fit annual predictions to all years...
    nll += square(wt_obs(styr,j)-wt_pre(styr,j))/(2.*square(sd_obs(styr,j)));
  }

  // subsequent years
  for (int i=styr+1;i<=endyr_r+3;i++)
  {
    wt_pre(i,age_st) = mnwt(age_st)*exp(sigma_coh*coh_eff(i-age_st));
    wt_pre(i)(age_st+1,age_end) = ++(wt_pre(i-1)(age_st,age_end-1) + wt_inc*exp(sigma_yr*yr_eff(i)));
    for (int j=age_st;j<=age_end;j++)
    {
      if (i <= endyr_r)
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
    sig_coh = exp(log_sd_coh);
    sig_yr  = exp(log_sd_yr );
    wt_last = wt_pre(endyr_r)  ; // *exp(sigma_coh*sigma_coh/2.  + sigma_yr*sigma_yr/2.);
    wt_cur  = wt_pre(endyr_r+1); // *exp(sigma_coh*sigma_coh/2.  + sigma_yr*sigma_yr/2.);
    wt_next = wt_pre(endyr_r+2); // *exp(sigma_coh*sigma_coh/2.  + sigma_yr*sigma_yr/2.);
    wt_yraf = wt_pre(endyr_r+3); // *exp(sigma_coh*sigma_coh/2.  + sigma_yr*sigma_yr/2.);
  }
  for (int i=1;i<=3;i++) wt_pre(endyr_r+i) *=  exp(sigma_coh*sigma_coh/2.  + sigma_yr*sigma_yr/2.);

REPORT_SECTION
  if (last_phase())
  {
  int retyr ;
  retyr = int (endyr - endyr_r);
  switch ( retyr ) 
  {
    case 0 : 
      report << "lof1"<<endl; report << 0 <<endl;
      report << "lof2"<<endl; report << 0 <<endl;
      report << "lof3"<<endl; report << 0 <<endl;
    break;

    case 1 : 
      report << "lof1"<<endl; report << pow(norm2(wt_obs(endyr_r+1)-wt_cur),.5) <<endl;
      report << "lof2"<<endl; report << 0 <<endl;
      report << "lof3"<<endl; report << 0 <<endl;
    break;

    case 2 : 
      report << "lof1"<<endl; report << pow(norm2(wt_obs(endyr_r+1)-wt_cur),.5) <<endl;
      report << "lof2"<<endl; report << pow(norm2(wt_obs(endyr_r+2)-wt_next),.5) <<endl;
      report << "lof3"<<endl; report << 0 <<endl;
    break;

    default : 
      report << "lof1"<<endl; report << pow(norm2(wt_obs(endyr_r+1)-wt_cur),.5) <<endl;
      report << "lof2"<<endl; report << pow(norm2(wt_obs(endyr_r+2)-wt_next),.5) <<endl;
      report << "lof3"<<endl; report << pow(norm2(wt_obs(endyr_r+3)-wt_yraf),.5) <<endl;
    break;
  }
  }
    /* */ 
  report <<"cur_yr" << endl;
  report << endyr_r << endl;
  report <<"data"   << endl;
  report <<wt_obs   << endl;
  report <<"W1"     << endl;
  dvector W1(age_st,age_end) ;
  W1 = wt_obs(endyr_r);
  REPORT(W1);
  report <<"W3"     << endl;
  dvector W3(age_st,age_end) ;
  W3.initialize();
  for (int i=1;i<=3;i++)
    W3 += wt_obs(endyr_r-i+1);
  report << W3/3. <<endl;
  W3.initialize();
  for (int i=1;i<=10;i++)
    W3 += wt_obs(endyr_r-i+1);
  report << "W10"<<endl<<W3/10. <<endl;
  report <<"W"<<endl;
  report << wt_pre<<endl;
  if (last_phase())
  {
    dmatrix resid(styr,endyr_r,age_st,age_end);
    for (int i=styr;i<=endyr_r;i++)
      resid(i) = value(wt_obs(i)-wt_pre(i));
    report <<"residuals"<<endl;
    report << resid << endl;
  }
  REPORT(sigma_yr);
  REPORT(yr_eff);
  REPORT(sigma_coh);
  report <<"cohort"<<endl;
  for (int j=age_st;j>=age_end;j--)
    report << styr-j <<" ";
  for (int i=styr+1;i<=endyr_r+3;i++)
    report << i-age_st <<" ";
  report<<endl;
  report <<"coh_eff"<<endl;
  for (int j=age_end;j>=age_st;j--)
    report << sigma_coh*coh_eff(styr-j) <<" ";
  for (int i=styr+1;i<=endyr_r+3;i++)
    report << sigma_coh*coh_eff(i-age_st) <<" ";
  report<<endl;
  REPORT(ages);
  REPORT(mnwt);
  // REPORT(winf);
  REPORT(K);
  REPORT(L1);
  REPORT(L2);
  // REPORT(t0);

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

