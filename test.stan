
 data {
  int<lower=1> N_obs; //number of observations
  int<lower=1> N_site; //number of sites
  int<lower=1> N_ttt_groups; //number of treatment groups

  array[N_obs] int site; //site ids
  array[N_obs] int ttt; //ttt group
  array[N_obs] int resp; //outcome 

}
transformed data { 
    matrix[N_obs, N_ttt_groups] X_matrix = rep_matrix(0, N_obs, N_ttt_groups);
    for (i in 1:N_obs){
      for (j in 2:N_ttt_groups){
        if (j==ttt[i]) X_matrix[i,j] = 1;
      
      }      
    }
    

}
parameters {
  real b0;
  vector[N_ttt_groups-1] beta_ttt; //beta coefficients for each treatment group
  vector[N_site] alpha_site_raw; //random effect coefficients for site
  cholesky_factor_corr[N_site] lkj_corr; //random effect correlation 
}

transformed parameters {   
 vector[N_site] b0_site = b0 +  lkj_corr * alpha_site_raw;//implies: random intercept for site is sampled from multivariate normal with mean 0 and 100 SD
 vector[N_ttt_groups] beta_trt;
    beta_trt[1] = 0;
    beta_trt[2:N_ttt_groups] = beta_ttt;

}
model {
  b0~normal(0,2);
  beta_trt[N_ttt_groups-1]~normal(0,2);
  lkj_corr ~ lkj_corr_cholesky(10); 
  alpha_site_raw~std_normal();
  
  
  resp ~ bernoulli_logit( b0_site[site] + X_matrix*beta_trt );
}

generated quantities {
    vector[N_ttt_groups] pred_prob_trt;
        pred_prob_trt = inv_logit(mean(b0_site) + beta_trt);

}
