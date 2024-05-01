
 data {
  int<lower=1> N_obs; //number of observations
  int<lower=1> N_site; //number of sites
  int<lower=1> N_trt_groups; //number of treatment groups

  array[N_obs] int site; //site ids
  array[N_obs] int trt; //trt group
  array[N_obs] int resp; //outcome 

}
transformed data { 
    matrix[N_obs, N_trt_groups-1] X_matrix = rep_matrix(0, N_obs, N_trt_groups-1);
    matrix[N_obs, N_site] Z_matrix = rep_matrix(0, N_obs, N_site);
    for (i in 1:N_obs){
      for (j in 2:N_trt_groups){
        if (j==trt[i]) X_matrix[i,j-1] = 1;
      }  
      for(k in 1:N_site)    
      if (k==site[i]) Z_matrix[i,k] = 1;
    }


}
parameters {
  real b0;
  vector[N_trt_groups-1] beta_trt; //beta coefficients for each treatment group
  vector[N_site] alpha_site_raw; //random effect coefficients for site
  cholesky_factor_corr[N_site] lkj_corr; //random effect correlation 
}

transformed parameters {   
 vector[N_site] b0_site = b0 + 0.2*lkj_corr * alpha_site_raw;//implies: random intercept for site is sampled from multivariate normal with mean 0 and 100 SD
 
 
}
model {
  b0~normal(0,2);
  beta_trt[N_trt_groups-1]~normal(0,2);
  lkj_corr ~ lkj_corr_cholesky(5); 
  alpha_site_raw~std_normal();
  
  resp ~ bernoulli_logit( X_matrix * beta_trt + Z_matrix*b0_site);
}

generated quantities {
    

    vector[N_trt_groups-1] pred_prob_trt;
        pred_prob_trt = inv_logit(mean(b0_site) + beta_trt);

array[N_obs] int ypred;

ypred = bernoulli_logit_rng(mean(b0_site) + X_matrix * beta_trt);

	int diff;
		if(max(beta_trt) == beta_trt[3]) diff=1;
		else diff=0;
	

}
