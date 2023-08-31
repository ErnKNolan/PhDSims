
# Making the outcome data function----------------------
# MakeSS function
# nreps = number of simulations
# t = number of treatment arms
# expdat = the generated dataset to simulate the results for
# rho = the intraclass correlation
# ... = the proportion of participants that have the event of interest 
# in each treatment arm, denoted as t1=, t2= etc
makeSS <- function(nreps,t,expdat,rho,...){
  prop <- list(...)
  comp <- (1:t)
  ss <- list()
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
  #when modelling, factor variable of ttt
  sigma2 <- (pi ^ 2) / 3
  theta <- sqrt((rho*sigma2)/(1-rho))
  names(theta)<-c("site.(Intercept)")
  ss <- simulate( ~ factor(ttt) + (1|site), nsim = nreps, family = binomial, 
                  newdata = expdat,newparams = list(beta=beta, theta=theta))
  return(ss)
}

#ss <- makeSS(100,4,outdat,0.02,t1=0.2,t2=0.35,t3=0.5,t4=0.65)

# fitting func-----------------
# fitSS function
# nreps = number of simulations==
# t = number of treatment arms
# dat = the generated dataset to model
# ss = the matrix of responses
fitSS <- function(nreps,t,dat,ss){
  comp <- expand.grid(1:t,1:t) %>% filter(Var1 != Var2 , Var1 < Var2)
  names(comp) <- c("Grp1","Grp2")
  results <- list()
  for(i in 1:nreps){
    dat$resp <- ss[,i]
    fit1 <- lme4::glmer(resp ~ factor(ttt) + (1 | site),family = binomial, 
                        data=dat,control=glmerControl(optimizer="bobyqa",
                                                      check.conv.singular = .makeCC(action = "ignore",  tol = 1e-4)))
    results[[i]] <- summary(lsmeans(fit1,pairwise~ttt,adjust="none"))$contrasts
    results[[i]]$type3 <- car::Anova(fit1)$"Pr(>Chisq)"
  }
  return(results)
}

#ssfit <- fitSS(100,4,outdat,ss)


#THE NEW FITSS
# fitfunc <- function(i,t,dat,ss,results){    
#   fit1 <- lme4::glmer(ss[,i] ~ factor(ttt) + (1 | site),family = binomial, 
#            data=dat,control=glmerControl(optimizer="bobyqa",
#            check.conv.singular = .makeCC(action = "ignore",  tol = 1e-4)))
#   res <- summary(lsmeans(fit1,pairwise~ttt,adjust="none"))$contrasts
#   res$type3 <- car::Anova(fit1)$"Pr(>Chisq)"
#   return(res)
# }
# fitSS <- function(i,t,dat,ss){
#   comp <- expand.grid(1:t,1:t) %>% filter(Var1 != Var2 , Var1 < Var2)
#   names(comp) <- c("Grp1","Grp2")
#   results <- fitfunc(i,t=4,dat=dat,ss=ss,results=results)
#   return(results)
# }
# 