#AUTHOR: Erin Nolan
#TITLE: PROJECT 2 FIT BAYESIAN MODELS ADAPTIVE
#PURPOSE: FUNCTION TO MAKE THE GENERATED DATASET AND FIT MODEL

#t = number of treatment groups
#expdat = the generated dataset (without response) to read in
#rho = intra-cluster correlation
#mod = model set by cmdstan_model
#... = the proportion for each treatment group i.e t1=0.2, t2=0.5 etc
testFull <- function(t,expdat,rho,mod,outdir,int_dat,...){
  tic()
  prop <- list(...)
  comp <- (1:t)
  o <- list()
  beta <- vector()
  for(i in 1:length(comp)){
    o[[i]] <- prop[[comp[i]]]/(1-prop[[comp[i]]])
  }
  beta <- append(beta,log(o[[1]]))
  for(i in 2:length(comp)){
    beta <- append(beta,log(o[[i]]/o[[1]]))
  }
  #intercept 01, beta = log(01), log(02/01), log(03/01),log(04/01)
  #when modelling, factor variable of trt
  sigma2 <- (pi ^ 2) / 3
  theta <- sqrt((rho*sigma2)/(1-rho))

  names(theta)<-c("ascendsite.(Intercept)")
  #fitting model
  results <- vector()
  #make a unique site number
  expdat$siteunique <- paste0(expdat$trt,expdat$site)
  expdat$ascendsite <- as.integer(factor(expdat$siteunique,levels=unique(expdat$siteunique)))
  
  resp <- suppressMessages(simulate.formula( ~ factor(trt) + (1|ascendsite), nsim = 1, family = binomial, 
                                             newdata = expdat,newparams = list(beta=beta, theta=theta)))
  
  resp_dat <- cbind(expdat,resp) 
  #need to remove them to merge, then regenerate
  resp_dat <- resp_dat %>% dplyr::select(-siteunique,-ascendsite)
  resp_dat <- resp_dat %>% rename(resp = sim_1)
  resp_dat <- merge(resp_dat,int_dat,by=c("iid","site","trt"),all.x=TRUE) %>%
    within(., resp <- ifelse(!is.na(resp.y), resp.y, resp.x)) %>%
    dplyr::select(-resp.x,-resp.y) %>% 
    arrange(trt,site)
  
  #making unique site again to make the merging above work (we drop them previously)
  resp_dat$siteunique <- paste0(resp_dat$trt,resp_dat$site)
  resp_dat$ascendsite <- as.integer(factor(resp_dat$siteunique,levels=unique(resp_dat$siteunique)))
  #preparing the data to read into stan
  N_obs <- dim(resp_dat)[1]
  N_site <- length(unique(resp_dat$ascendsite))
  N_trt_groups <- length(unique(resp_dat$trt))
  data <- list(N_obs = N_obs, N_site = N_site, N_trt_groups = N_trt_groups, 
               site = resp_dat$ascendsite, trt = as.numeric(resp_dat$trt), resp = resp_dat$resp)
  
    res <- mod$sample(
      data = data, 
      init = 0,
      iter_warmup = 750,
      iter_sampling = 750,
      chains = 4, 
      parallel_chains = 1,
      adapt_delta = 0.99,
      refresh = 0, 
      max_treedepth=12,
      output_dir=outdir
      
    )
    print(j)
    time <- toc()
    time <- time$toc - time$tic
    
    results <- list(data.frame(res$summary(variables=c("pred_prob_trt","beta_trt","pp_trt2","pp_trt3","pp_trt4","probd_trt2","probd_trt3","probd_trt4"))))
  return(results)
  
}
