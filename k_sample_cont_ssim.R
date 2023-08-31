MakeSSNorm <- function(nreps,t,expdat,rho,sigma,...){
  mean <- list(...)
  comp <- (1:t)
  ss <- list()
  means <- list()
  beta <- vector()
  beta <- append(beta,mean[[1]])
  for(i in 2:length(comp)){
    beta <- append(beta,(mean[[i]]-mean[[1]]))
  }
  #when modelling, factor variable of ttt
  sigma2 <- sigma^2
  theta <- sqrt((rho*sigma2)/(1-rho))
  names(theta)<-c("site.(Intercept)")
  ss <- simulate( ~ factor(ttt) + (1|site), nsim = nreps, family = "gaussian", 
                  newdata = outdat,newparams = list(beta=beta, sigma=sigma,theta=theta))
  return(ss)
}

#maybe need to have sigma specified in the function for gaussian models?

testdat <- MakeSSNorm(10,5,outdat,0.02,sigma=2,t1=3.1,t2=5.2,t3=4.3,t4=4,t5=1.3)


fitSS <- function(nreps,t,dat,ss){
  comp <- expand.grid(1:t,1:t) %>% filter(Var1 != Var2 , Var1 < Var2)
  names(comp) <- c("Grp1","Grp2")
  results <- list()
  dat <- outdat
  for(i in 1:nreps){
    dat$resp <- testdat[,i]
    fit1 <- lme4::lmer(resp ~ factor(ttt) + (1 | site),
                        data=dat)
    results[[i]] <- summary(lsmeans(fit1,pairwise~ttt,adjust="none"))$contrasts
    results[[i]]$type3 <- car::Anova(fit1)$"Pr(>Chisq)"
  }
  return(results)
}

fits <- fitSS(10,5,dat,testdat)
result <- selectRes(fits,"1 - 2")
