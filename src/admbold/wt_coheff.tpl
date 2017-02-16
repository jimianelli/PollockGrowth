//----------------------------
// Goal to predict fishery weights at age
// Data include survey and fishery data
//----------------------------
//----------------------------

DATA_SECTION
  init_int cur_yr ;// Retrospective years to test over (0 if current year)
  int retro;
  init_int styr;
  init_int endyr;
  init_int ndat;
  int nscale_parm;
  !! nscale_parm = ndat-1;
  init_ivector fyr(1,ndat);
  init_ivector lyr(1,ndat);
  init_int age_st;
  init_int age_end;
  int nages;
  init_3darray wt_obs(1,ndat,fyr,lyr,age_st,age_end);
  init_3darray sd_obs(1,ndat,fyr,lyr,age_st,age_end);
  int phase_d_scale;
  vector ages(age_st,age_end);
 LOCAL_CALCS
  // Need to reset fyr and lyr to be cur_yr covered (for retrospective fitting)
  retro = endyr-cur_yr-2;
  for (int i=1;i<=ndat;i++) lyr(i) -= retro;
  // here
  if (ndat>1) phase_d_scale = 1; else phase_d_scale = -1;
  nages = age_end - age_st + 1;
  for (int i=age_st;i<=age_end;i++) ages(i) = double(i);
  // COUT(sd_obs); 
  cout<< ndat <<" "<<fyr<<" "<<styr <<" "<< endyr <<" "<< age_st <<" "<< age_end<<endl;
  //  exit(1);
 END_CALCS
  
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
  init_number log_K(3);
  init_matrix d_scale(1,nscale_parm,age_st,age_end,phase_d_scale);
  // init_number log_t0(3);
  // Predicted weight matrix
  3darray wt_hat(1,ndat,fyr,lyr,age_st,age_end);
  matrix  wt_pre(styr,endyr,age_st,age_end);
  number alpha;
  vector mnwt(age_st,age_end);
  vector wt_inc(age_st,age_end-1);
  sdreport_number sigma_coh;
  sdreport_number sigma_yr;
  // sdreport_number winf;
  sdreport_number K;
  // sdreport_number t0;
  sdreport_vector wt_cur(age_st,age_end);
  sdreport_vector wt_next(age_st,age_end);
  sdreport_vector wt_yraf(age_st,age_end);

  init_bounded_number log_sd_coh(-5,10,5);
  init_bounded_number log_sd_yr(-5,10,-6);
  random_effects_vector coh_eff(styr,endyr,4);
  random_effects_vector yr_eff(styr,endyr,-4);
  //init_vector coh_eff(styr,endyr,3);
  //init_vector yr_eff(styr,endyr,4);

  objective_function_value nll;

PROCEDURE_SECTION
  sigma_coh    = mfexp(log_sd_coh);
  sigma_yr     = mfexp(log_sd_yr );
  // winf         = exp(log_winf);
  K            = mfexp(log_K);
  alpha        = mfexp(log_alpha);
  // t0           = (log_t0);
  // mnwt         = winf*(pow(1.-mfexp(-K*(ages-t0)),3));
  for (int j=age_st;j<=age_end;j++)
  {
    mnwt(j)    = alpha * pow(L1 + (L2-L1)*(1.-pow(K,double(j-age_st))) / (1.-pow(K,double(nages-1))) ,3);
  }
  wt_inc       = --mnwt(age_st+1,age_end) - mnwt(age_st,age_end-1);

  // Initialize first year
  wt_pre(styr)    = mnwt;
  /* for (int j=age_st;j<=age_end;j++)
  {
    if (j>age_st) 
      wt_pre(styr,j) = wt_pre(styr,j-1) + wt_inc(j-1)*exp(sigma_yr*yr_eff(styr));
    int icoh        = styr-j;
    wt_pre(styr,j) *= mfexp(sigma_coh*coh_eff(icoh));
  }
  */

  // subsequent years
  for (int i=styr+1;i<=endyr;i++)
  {
    wt_pre(i,age_st) = mnwt(age_st)*mfexp(square(sigma_coh)/2.+sigma_coh*coh_eff(i));
    if (last_phase())
      wt_pre(i)(age_st+1,age_end) = ++(wt_pre(i-1)(age_st,age_end-1) + wt_inc*mfexp(square(sigma_yr)/2. + sigma_yr*yr_eff(i)));
    else
      wt_pre(i)(age_st+1,age_end) = ++(wt_pre(i-1)(age_st,age_end-1) + wt_inc*mfexp(sigma_yr*yr_eff(i)));
  }
    // Fit global mean to all years...
  for (int h = 1;h<=ndat;h++)
  {
  	if (h>1) 
  	  for (int i=fyr(h);i<=lyr(h);i++)
  		  wt_hat(h,i) = elem_prod(d_scale(h-1) , wt_pre(i));
  	else
  	  for (int i=fyr(h);i<=lyr(h);i++)
  		  wt_hat(h,i) = wt_pre(i);

    /*
    for (int j=age_st;j<=age_end;j++)
    {
      nll += square(wt_obs(h,fyr(h),j)-mnwt(j))/(2.*square(sd_obs(h,fyr(h),j)));
    // Fit annual predictions to all years...
      nll += square(wt_obs(h,fyr(h),j)-wt_hat(h,fyr(h),j))/(2.*square(sd_obs(h,fyr(h),j)));
    }
    */
    for (int i=fyr(h);i<=lyr(h);i++)
    {
      // cout <<i<<" "<<h<<" "<<nll<<" "<<endl;
      for (int j=age_st;j<=age_end;j++)
      {
        nll += square(wt_obs(h,i,j) - mnwt(j))      /(2.*square(sd_obs(h,i,j)));
        nll += square(wt_obs(h,i,j) - wt_hat(h,i,j))/(2.*square(sd_obs(h,i,j)));
      }
    }
  }
  nll += 0.5*norm2(coh_eff);
  nll += 0.5*norm2( yr_eff);

  if (sd_phase())
  {
    wt_cur  = wt_pre(cur_yr); // *mfexp(sigma_coh*sigma_coh/2.  + sigma_yr*sigma_yr/2.);
    wt_next = wt_pre(endyr-1); // *exp(sigma_coh*sigma_coh/2.  + sigma_yr*sigma_yr/2.);
    wt_yraf = wt_pre(endyr); // *exp(sigma_coh*sigma_coh/2.  + sigma_yr*sigma_yr/2.);
  }
  // for (int i=3;i>=0;i--) wt_pre(endyr-i) *=  exp(sigma_coh*sigma_coh/2.  + sigma_yr*sigma_yr/2.);

REPORT_SECTION
  if (last_phase())
  {
  int retyr ;
  retyr = int (endyr - endyr);
  switch ( retyr ) 
  {
    case 0 : 
      report << "lof1"<<endl; report << 0 <<endl;
      report << "lof2"<<endl; report << 0 <<endl;
      report << "lof3"<<endl; report << 0 <<endl;
    break;

    case 1 : 
      report << "lof1"<<endl; report << pow(norm2(wt_obs(1,endyr-3)-wt_cur),.5) <<endl;
      report << "lof2"<<endl; report << 0 <<endl;
      report << "lof3"<<endl; report << 0 <<endl;
    break;

    case 2 : 
      report << "lof1"<<endl; report << pow(norm2(wt_obs(1,endyr-3)-wt_cur),.5) <<endl;
      report << "lof2"<<endl; report << pow(norm2(wt_obs(1,endyr-2)-wt_next),.5) <<endl;
      report << "lof3"<<endl; report << 0 <<endl;
    break;

    default : 
      report << "lof1"<<endl; report << pow(norm2(wt_obs(1,endyr-3)-wt_cur),.5) <<endl;
      report << "lof2"<<endl; report << pow(norm2(wt_obs(1,endyr-2)-wt_next),.5) <<endl;
      report << "lof3"<<endl; report << pow(norm2(wt_obs(1,endyr-1)-wt_yraf),.5) <<endl;
    break;
  }
  }
    /* */ 
  report <<"endyr" << endl;
  report << endyr << endl;
  REPORT(retro);
  REPORT(lyr);
  REPORT(cur_yr);
  report <<"data"   << endl;
  report <<wt_obs   << endl;
  report << "yr"<<endl; 
  for (int iyr=styr;iyr<=endyr;iyr++) report <<iyr<<" "<<endl;
  report <<"wt_pre"<<endl;
  report << wt_pre<<endl;
  if (last_phase())
  {
    for (int h = 1;h<=ndat;h++)
    {
      dmatrix resid(fyr(h),lyr(h),age_st,age_end);
      for (int i=fyr(h);i<=lyr(h);i++)
        resid(i) = value(wt_obs(h,i)-wt_hat(h,i));
      report <<"residuals_"<<h<<endl;
      report << resid << endl;
    }
  }
  REPORT(sigma_yr);
  REPORT(yr_eff);
  REPORT(sigma_coh);
  report <<"cohort"<<endl;
  // for (int j=age_st;j>=age_end;j--)
    // report << styr-j <<" ";
  for (int i=styr+1;i<=endyr;i++)
    report << i-age_st <<" ";
  report<<endl;
  report <<"coh_eff"<<endl;
  // for (int j=age_end;j>=age_st;j--)
    // report << sigma_coh*coh_eff(styr-j) <<" ";
  for (int i=styr;i<=endyr;i++)
    report << sigma_coh*coh_eff(i) <<" ";
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


TOP_OF_MAIN_SECTION
  arrmblsize = 1000000;
  gradient_structure::set_MAX_NVAR_OFFSET(10000);
