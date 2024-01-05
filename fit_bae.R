#AUTHOR: Erin Nolan
#TITLE: PROJECT 2 FIT BAYESIAN MODELS
#PURPOSE: FUNCTION TO MAKE THE GENERATED DATASET

#t = number of treatment groups
#expdat = the generated dataset (without response) to read in
#rho = intra-cluster correlation
#mod = model set by cmdstan_model
#... = the proportion for each treatment group i.e t1=0.2, t2=0.5 etc
testss <- function(t,expdat,rho,mod,outdir,...){
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
  #beta <- c(log(o[[1]]),log(o[[2]]/o[[1]]),log(o[[3]]/o[[1]]),log(o[[4]]/o[[1]]))
  #intercept 01, beta = log(01), log(02/01), log(03/01),log(04/01)
  #when modelling, factor variable of trt
  sigma2 <- (pi ^ 2) / 3
  theta <- sqrt((rho*sigma2)/(1-rho))
  names(theta)<-c("site.(Intercept)")
  #fitting ss
  comp <- expand.grid(1:t,1:t) %>% filter(Var1 != Var2 , Var1 < Var2)
  names(comp) <- c("Grp1","Grp2")
  results <- vector()
  
  
  return <- tryCatch({

    resp <- suppressMessages(simulate.formula( ~ factor(trt) + (1|site), nsim = 1, family = binomial, 
                                              newdata = expdat,newparams = list(beta=beta, theta=theta)))
    resp <- as.vector(resp[,1])
    N_obs <- dim(expdat)[1]
    N_site <- length(unique(expdat$site))
    N_trt_groups <- length(unique(expdat$trt))
    data <- list(N_obs = N_obs, N_site = N_site, N_trt_groups = N_trt_groups, 
                 site = expdat$site, trt = as.numeric(expdat$trt), resp = resp)
    
    res <- mod$sample(
      data = data, 
      init = 0,
      iter_warmup = 250,
      iter_sampling = 250,
      chains = 4, 
      parallel_chains = 1,
      adapt_delta = 0.8,
      refresh = 0, 
      max_treedepth=12,
      output_dir=outdir
      
    )
    print(j)
    time <- toc()
    time <- time$toc - time$tic
    results <- list(data.frame(res$summary(variables=c("pred_prob_trt","pp_trt2","pp_trt3","pp_trt4","beta_trt")),time=time))
    
  },
  
  error=function(e) { message(conditionMessage(e)) 
    res <- data.frame(variable=NA,mean=NA,median=NA,sd=NA,mad=NA,q5=NA,q95=NA,rhat=NA,ess_bulk=NA,ess_tail=NA,time=NA)
    results <- list(res)
  })
  
  return(results)
  
}

