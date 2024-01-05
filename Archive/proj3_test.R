
#t = number of treatment groups
#expdat = the generated dataset (without response) to read in
#rho = intra-cluster correlation
#mod = model set by cmdstan_model
#... = the proportion for each treatment group i.e t1=0.2, t2=0.5 etc

j <- 91
expdat=intclusters
t=4
mod=mod
outdir=outdir

properties <- expand.grid(trt_eff_scen = c(1,2,3), ctrl_prop = c(0.1), icc = c(0.05,0.2), n_per_k = c(25,50,75,100), k = c(5,10),int = c(0,0.5))
#bind to properties
properties <- rbind(properties) %>%
  mutate(f1 = case_when(trt_eff_scen == 1 ~ 0.2,
                        trt_eff_scen == 2 ~ 0.15,
                        trt_eff_scen == 3 ~ 0),
         f2 = case_when(trt_eff_scen == 1 ~ 0.1,
                        trt_eff_scen == 2 ~ 0.05,
                        trt_eff_scen == 3 ~ 0),
         t1 = ctrl_prop+f1+f2+(int*(f1+f2)),
         t2 = ctrl_prop+f1,
         t3 = ctrl_prop+f2,
         t4 = ctrl_prop)

rho=properties$icc[j]
t1=properties$t1[j]
t2=properties$t2[j]
t3=properties$t3[j]
t4=properties$t4[j]


resp_clus <- cbind(intclusters,resp)
resp_clus <- resp_clus %>%
  mutate(f1 = ifelse(trt %in% c(1, 2), 1, 0),
         f2 = ifelse(trt %in% c(1, 3), 1, 0))

testInterim <- function(t,expdat,rho,mod,outdir,int_dat,...){
  tic()
  prop <- list(t1,t2,t3,t4)
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
  names(theta)<-c("site.(Intercept)")
  #fitting model
  results <- vector()
  return <- tryCatch({
    
    resp <- suppressMessages(simulate.formula( ~ factor(trt)  + (1|site), nsim = 1, family = binomial, 
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
      adapt_delta = 0.99,
      refresh = 0, 
      max_treedepth=12,
      output_dir=outdir
      
    )
    print(j)
    time <- toc()
    time <- time$toc - time$tic
    #pp_trt1 etc are the predictive prob that the treatment has the largest beta
    results <- list(data.frame(res$summary(variables=c("pred_prob_trt","pp_trt1","pp_trt2","pp_trt3",
                                                       "fu_trt1","fu_trt2","fu_trt3")),time=time),resp=list(resp))
    
  },
  
  error=function(e) { message(conditionMessage(e)) 
    res <- list(data.frame(variable=NA,mean=NA,median=NA,sd=NA,mad=NA,q5=NA,q95=NA,rhat=NA,ess_bulk=NA,ess_tail=NA,time=NA))
  })
  
  return(results)
  
}

