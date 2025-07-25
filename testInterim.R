
#t = number of treatment groups
#expdat = the generated dataset (without response) to read in
#rho = intra-cluster correlation
#mod = model set by cmdstan_model
#... = the proportion for each treatment group i.e t1=0.2, t2=0.5 etc
testInterim <- function(t,expdat,rho,mod,outdir,int_dat,...){
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
  
  resp <- as.vector(resp[,1])
  N_obs <- dim(expdat)[1]
  N_site <- length(unique(expdat$ascendsite))
  N_trt_groups <- length(unique(expdat$trt))
  data <- list(N_obs = N_obs, N_site = N_site, N_trt_groups = N_trt_groups, 
               site =  expdat$ascendsite, trt = as.numeric(expdat$trt), resp = resp)   
  
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
    #pp_trt1 etc are the predictive prob that the treatment has the largest beta
    results <- list(data.frame(res$summary(variables=c("pred_prob_trt","beta_trt","pp_trt2","pp_trt3","pp_trt4","probd_trt2","probd_trt3","probd_trt4")),time=time),resp=list(resp))
  return(results)
  
}
