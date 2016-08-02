DATA_SECTION
	init_number binit;
	init_int syr;
	init_int nyr;
	init_int sage;
	init_int nage;
	init_matrix data(syr,nyr,sage,nage);
	init_matrix obs_std(syr,nyr,sage,nage);	

	int byr;
	!! byr = syr - (nage-sage+1);
	!! cout<<byr<<endl;
	!! cout<<syr<<endl;
	!! cout<<nyr<<endl;
	

	vector age(sage,nage);
	!!age.fill_seqadd(sage,1);

INITIALIZATION_SECTION
	log_winf  2.302585093;
	log_vonb -1.049822124;
	log_to   -0.693147181;
	b binit;

PARAMETER_SECTION
	init_number log_winf(1);
	init_number log_vonb(1);
	init_number log_to(2);
	init_bounded_number b(2,4,-4);


	init_bounded_number log_sd_year(-30,10,3);
	init_bounded_number log_sd_cohort(-30,10,3);

	sdreport_matrix W_sd(nyr+1,nyr+3,sage,nage);


	// init_bounded_vector year_effect(syr,nyr,-5.0,5.0,2);
	// init_bounded_vector cohort_effect(byr,nyr,-5.0,5.0,2);
	random_effects_vector year_effect(syr,nyr+3,2);
	random_effects_vector cohort_effect(byr,nyr+3,2);

	objective_function_value f;

	number winf;
	number vonb;
	number to;
	number sd_year;
	number sd_cohort;

	number nll_year_effect;
	number nll_cohort_effect;

	vector wa(sage,nage);

	matrix W(syr,nyr+3,sage,nage);	// Predicted Weight
	matrix E(syr,nyr+3,sage,nage);	// Random effects
	matrix R(syr,nyr,sage,nage); 	// Residuals


PROCEDURE_SECTION

	winf = exp(log_winf);
	vonb = exp(log_vonb);
	to   = - exp(log_to);
	sd_year = mfexp(log_sd_year);
	sd_cohort = mfexp(log_sd_cohort);

	wa   = winf * pow(1.0 - exp( -vonb*(age-to)),b );
	//COUT(wa);

	// Initial Cohort effects 
	int iyr = byr;
	for( int j = sage; j <= nage; j++)
	{
		E(syr,j) = exp(sd_cohort * cohort_effect(iyr++) + sd_year * year_effect(syr));
	}
	

	for( int i = syr+1; i <= nyr+3; i++)
	{
		for( int j = sage; j <= nage; j++)
		{
			if( j == sage )
			{
				E(i,j) = exp(sd_cohort * cohort_effect(i));
			} 
			else 
			{
				E(i,j) = E(i-1,j-1) * exp(sd_year * year_effect(i));
			}
		}
	}

	// predicted data where the weight increment for the vb model is given by:
	// dw=%e^(-b*delta*k-a1*b*k) * ((%e^(delta*k+a1*k)-%e^(k*to))^b
	//    -%e^(b*delta*k)*(%e^(a1*k)-%e^(k*to))^b)*winf
	double delta = 1;  // annual time step.
	dvariable t1,t2,t3,t4,t5,t6;
	dvariable k,dw;
	k = vonb;
	W(syr) = elem_prod(wa,E(syr));
	R(syr) = data(syr) - W(syr);
	//cout<<W<<endl;
	for( int i = syr+1; i <= nyr+3; i++)
	{
		for( int j = sage; j <= nage; j++)
		{
			if( j == sage )
			{
				W(i,j) = wa(j) * E(i,j);
			} 
			else if(j > sage)
			{
				double a1 = double(j)-1.0;
				t1 = exp(-b*delta*k-a1*b*k);
				t2 = pow(exp(delta*k+a1*k)-exp(k*to),b);
				t3 = exp(b*delta*k) * pow(exp(a1*k)-exp(k*to),b);
				t5 = t2 - t3;
				dw = winf * t1 * t5;
				
				//cout << j << "\t" << dw<< endl;
				W(i,j) = W(i-1,j-1) + dw * E(i,j);
			}
			// Residual
			if (i<= nyr) 
		  	R(i,j) = data(i,j) - W(i,j);
		}
	}

	// Compute the objective function
	double log2pi = log(2.0*M_PI);
	dvariable constant;
	dvariable nloglike;
	dvariable variance;
	nloglike=0;
	int n = size_count(age);
	for( int i = syr+1; i <= nyr; i++){
		for( int j = sage; j <= nage; j++){
			variance  = square(obs_std(i,j));
			constant  = 0.5 * log2pi + 0.5 * log(variance);
			nloglike += constant + 1.0/(2*variance)*square(R(i,j));
		}
	}

	// pdf for year effect
	dvariable nll_year_effect = 0.0;
	n = 0.5*size_count(year_effect);
	nll_year_effect = n * log2pi + n * log(square(sd_year)) 
					+ 0.5 * norm2(year_effect) / square(sd_year);
	// nll_year_effect = 0.5 * norm2(year_effect);

	// pdf for cohort effect
	dvariable nll_cohort_effect = 0.0;
	n = 0.5*size_count(cohort_effect);
	nll_cohort_effect = n * log2pi + n * log(square(sd_cohort))
					  + 0.5 * norm2(cohort_effect) / square(sd_cohort);
	// nll_cohort_effect = 0.5 * norm2(cohort_effect);

	f = nloglike + nll_year_effect + nll_cohort_effect;

	if (sd_phase()) for (int i=nyr+1;i<=nyr+3;i++) W_sd(i) = W(i);

	//cout<<E<<endl;
	//cout<<W<<endl;
	//exit(1);

REPORT_SECTION

	REPORT(winf);
	REPORT(vonb);
	REPORT(to);
	REPORT(sd_year);
	REPORT(sd_cohort);
	REPORT(data);
	REPORT(W);
	REPORT(R);
	REPORT(E);


FINAL_SECTION
	for (int i=nyr+1;i<=nyr+3;i++) 
	{
		cout << (W_sd(i)) <<endl;
		cout << (W_sd.sd(i))<<endl;
	}


GLOBALS_SECTION
	#undef REPORT
	#define REPORT(object) report << #object "\n" << setw(8) \
	<< setprecision(6) << setfixed() << object << endl;

	
	 #undef COUT
	 #define COUT(object) cout << #object "\n" << setw(6) \
	 << setprecision(5) << setfixed() << object << endl;


