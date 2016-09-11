#include <TMB.hpp>
template <class Type> 
Type square(Type x){return x*x;}

template<class Type>
Type objective_function<Type>::operator()()
{
 DATA_INTEGER(cur_yr);
 DATA_INTEGER(styr);
 DATA_INTEGER(endyr);
 DATA_INTEGER(fyr_fsh);
 DATA_INTEGER(lyr_fsh);
 DATA_INTEGER(f_age);
 DATA_INTEGER(l_age);
 DATA_VECTOR(ages);
 DATA_MATRIX(wt_obs);
 DATA_MATRIX(wt_sd);
 int nyrs;
 int nyrs_fsh;
 int retr;
 int nages ;
 retr     = (endyr-cur_yr)-2;
 nyrs     = endyr - styr + 1 - retr;
 lyr_fsh  = lyr_fsh - retr;
 nyrs_fsh = lyr_fsh - fyr_fsh + 1;
 nages    = l_age - f_age + 1;

 // Rcout << lyr_fsh<<" "<<nyrs_fsh<<" "<<retr;
 // These are the parameters (three are vectors; one is a scalar)
 PARAMETER(L1);
 PARAMETER(L2);
 PARAMETER(log_K);
 PARAMETER(log_alpha);
 PARAMETER(log_sigma_coh);
 PARAMETER(log_sigma_yr);
 PARAMETER_VECTOR(coh_eff);
 PARAMETER_VECTOR(yr_eff);
 Type sigma_coh = exp(log_sigma_coh);
 Type sigma_yr  = exp(log_sigma_yr);
 Type K         = exp(log_K);
 Type alpha     = exp(log_alpha);
 matrix<Type> wt_prj(3,nages);
 matrix<Type> wt_pre(nyrs,nages);
 matrix<Type> wt_hat(nyrs_fsh,nages);
 vector<Type> mnwt(nages);
 vector<Type> mnlen(nages);
 vector<Type> likecomp(4);
 vector<Type> wt_inc(nages-1);
 // vector<Type> yrs(nyrs);

 // Predictions and likelihoods
  for (int j=0;j<nages;j++)
  {
    Type age_tmp = Type (j);
    // Rcout <<mnwt.size()<<"\n";
    mnlen[j]   =  L1 + (L2-L1)*(1.-pow(K,age_tmp)) / (1.-pow(K,Type(nages-1))) ;
    mnwt[j]    = alpha * pow(mnlen(j),3);

    if(j > 0) wt_inc(j-1)  = mnwt(j) - mnwt(j-1);
    wt_pre(0,j)    = mnwt(j);
  }
  for (int i=0;i<4;i++) likecomp(i) = 0.; 
  for (int i=0;i<nyrs;i++)
  {
    wt_pre(i,0) = mnwt(0)*exp(square(sigma_coh)/2. + sigma_coh*coh_eff(i));
    if (i>0)
    {
      for (int j=1;j<nages;j++)
      {
        wt_pre(i,j) = wt_pre(i-1,j-1) + wt_inc(j-1)*exp(square(sigma_yr)/2. + sigma_yr*yr_eff(i));
      }
    }
    likecomp(0) += 0.5*square(coh_eff(i));
    likecomp(1) += 0.5*square( yr_eff(i));
  }
  // Now assign wt_pre to years of observations...
  for (int i=fyr_fsh;i<=lyr_fsh;i++)
  {
    int itmp   = i-styr;
    int i2     = i-fyr_fsh;
    for (int j=0;j<nages;j++)
    {
      wt_hat(i2,j) = wt_pre(itmp,j);
      // Rcout<<wt_hat(i2,j)<<" "<<itmp<<" "<<i2<<" "<<j<<"\n";
    }
  }
  /* */
 // -sum(dnorm(wt_obs,wt_pre,wt_sd,true));
  for (int i=0;i<nyrs_fsh;i++)
  {
    for (int j=0;j<nages;j++)
    {
      likecomp(2) += square(wt_obs(i,j) - mnwt(j))    /(2.*square(wt_sd(i,j)));
      likecomp(3) += square(wt_obs(i,j) - wt_hat(i,j))/(2.*square(wt_sd(i,j)));
    }
  }
  for (int i=0;i<3;i++)
  {
    for (int j=0;j<nages;j++)
    {
      int itmp   = nyrs - (3 - i);
      wt_prj(i,j) = wt_pre(itmp,j);
    }
  }
  /*
  */
  REPORT(wt_pre);
  REPORT(wt_obs);
  REPORT(wt_hat);
  REPORT(L1);
  REPORT(L2);
  REPORT(mnlen);
  REPORT(mnwt);
  REPORT(likecomp);
  REPORT(wt_prj  );
  REPORT(retr    );   
  REPORT(nyrs    );   
  REPORT(fyr_fsh );  
  REPORT(lyr_fsh );  
  REPORT(nyrs_fsh); 
  REPORT(nages   );


 // Provide the standard error of Sigma
 ADREPORT(sigma_coh);
 ADREPORT(sigma_yr);
 ADREPORT(wt_prj);


  Type nll = sum(likecomp);
  return nll; 
 }


