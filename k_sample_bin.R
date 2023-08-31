#------------------------------------------------------------------------------------------------------------
#File Name: k_sample_bin
#Purpose: Power calculation for a k sample binomial test

#Authors: Erin Nolan
#Organisation: Hunter Medical Research Institute
#------------------------------------------------------------------------------------------------------------

#Dependencies------------------------------------------------------------------------------------------------
  
pacman::p_load(tidyverse,purrr,lme4,MASS)

#Function(s)-------------------------------------------------------------------------------------------------

#k = number of arms
#proportions for each arm must be specified
#nperarm = sample size per arm
#numsims = number of simulations
kArmBin <- function(nperarm,k,numsims,...){
prop <- list(...) #the proportion in each arm must be specified
y <- list()
phat <- list()
power <- vector()
mean_pval <- vector()
#make the distribution for each group
for(i in 1:k){
  y[[i]] <- rbinom(numsims,nperarm,prop[[i]])
  phat[[i]] <- y[[i]]/nperarm
}
#get the power for each group, all comparisons done
test <- expand.grid(1:k,1:k) %>% filter(Var1 != Var2 , Var1 < Var2)
names(test) <- c("Grp1","Grp2")
for(i in 1:nrow(test)){
  name <- paste(test[i,]$Grp1,"vs",test[i,]$Grp2)
  ctrl <- as.numeric(test[i,]$Grp1)
  trt <- as.numeric(test[i,]$Grp2)
  SEdiff <- sqrt(((phat[[ctrl]]*(1-phat[[ctrl]]))/nperarm) + 
                   ((phat[[trt]]*(1-phat[[trt]]))/nperarm))
  Z <- abs(phat[[trt]] - phat[[ctrl]])/SEdiff
  p <-  2*pnorm(q=Z, lower.tail=FALSE)
  x <- mean(p<0.05)
  power[i] <- x
  mean_pval[i] <- mean(p)
}
res <- data.frame(test,power,mean_pval)
return(res)

}
x <- kArmBin(nperarm=50,numsims=1000,k=5,p1=0.2,p2=0.3,p3=0.4,p4=0.5,p5=0.7)

#---------------------------
#t = number of arms
#nid = number of participants per arm - can make this specified
#per arm if wanted
#specify the number of clusters for each t
makeClusters <- function(t,nid,...){
int <- list()
clusters <- list(...) #the number of clusters in each group
nclusters <- cumsum(clusters) #cumulative n for the for loop
for(i in 1:t){
  if (i==1){ #first group has cluster starting at 1
int[[i]] <- expand.grid(iid = c(1:nid), site = seq(1,nclusters[i],1), ttt = i)
  }
  else {
int[[i]] <- expand.grid(iid = c(1:nid), site = seq(nclusters[i-1]+1,nclusters[i],1), ttt = i)
  }
  }
groups <- map_df(int, ~as.data.frame(.)) #unlist them and make them a dataframe
return(groups)
}
#outdat <- makeClusters(t=3,nid=50,t1=3,t2=4,t3=3)

#---------------------------
# specify the parameters the simulate function:   
# "theta" - in the case of single variance terms, 
# that's just the standard deviation (on the logit scale) 
# of each random effect: e.g. theta=c(1,0.5) (among-individual variation is 4x among-observation variance). 
# "beta" is the fixed-effects parameters, in this case (intercept,treat) also all on the logit scale.
makeSS <- function(t,expdat,rho,nreps,...){
prop <- list(...)
comp <- expand.grid(1:t,1:t) %>% filter(Var1 != Var2 , Var1 < Var2)
names(comp) <- c("Grp1","Grp2")
ss <- list()
for(i in 1:nrow(comp)){
  dat <- expdat[which(expdat$ttt %in% c(comp[i,]$Grp1,comp[i,]$Grp2)),]
  o2 <- prop[[comp[i,]$Grp2]]/(1-prop[[comp[i,]$Grp2]])
  o1 <- prop[[comp[i,]$Grp1]]/(1-prop[[comp[i,]$Grp1]])
  beta <- c(log(o2),log(o1/o2))
  #intercept 01, beta = log(01), log(02/01), log(03/01),log(04/01)
  #when modelling, factor variable of ttt
  names(beta)<-c("(Intercept)","ttt")
  sigma2 <-(pi ^ 2) / 3
  theta <-sqrt((rho*sigma2)/(1-rho))
  names(theta)<-c("site.(Intercept)")
  ss[[i]] <- simulate(~ttt + (1|site), nsim = nreps, family = binomial, 
                 newdata = dat, newparams = list(theta = theta, beta = beta))
  names(ss[i]) <- paste(comp[i,]$Grp1,"vs",comp[i,]$Grp2)
}
return(ss)
}
#ss <- makeSS(3,outdat,0.05,100,t1=0.1,t2=0.2,t3=0.3)

#----------------- fitting func
fitSS <- function(nreps,t,dat,ss){
  test <- dat
  comp <- expand.grid(1:t,1:t) %>% filter(Var1 != Var2 , Var1 < Var2)
  names(comp) <- c("Grp1","Grp2")
  results <- list()
  test2 <- list()
  for(j in 1:nrow(comp)){
    test2[[j]] <- test[which(test$ttt %in% c(comp[j,]$Grp1,comp[j,]$Grp2)),]
    results[[j]] <- list()
    nclus <- bind_rows(test2[[j]]) %>% group_by(ttt) %>% summarise(nclus = n_distinct(site))
    for(i in 1:nreps){
      test2[[j]]$resp <- ss[[j]][,i]
      if(nclus[1,2] < 4 | nclus[2,2] < 4){
        res<-coef(summary(MASS::glmmPQL(resp ~ ttt ,
                                        random = ~ 1 | site, family = binomial,
                                        data=test2[[j]],verbose=F)))["ttt", ]
        names(res)<-c("est","se","df","t","p")
        res <- res[c(1,2,4,5)]
        results[[j]][[i]] <- res
      }
      else{
        fit1 <- lme4::glmer(resp ~ ttt + (1 | site), 
                            family = binomial, data=test2[[j]], 
                            control=glmerControl(optimizer="bobyqa", 
                                                 check.conv.singular = .makeCC(action = "ignore",  tol = 1e-4)))
        res<-coef(summary(lme4::refit(fit1, ss[[j]][[i]])))["ttt", ]
        names(res)<-c("est","se","z","p")
        results[[j]][[i]] <- res
      }
    }
    results[[j]] <- bind_rows(results[[j]])
  }
  return(results)
}

#ss2 <- bind_rows(ss)
##for(i in 1:nreps){
#outdat$resp <- ss2[,i]
#MASS::glmmPQL(resp ~ ttt ,
#              random = ~ 1 | site, family = binomial,
#              data=outdat,verbose=F)
#}

outdat <- makeClusters(t=4,nid=30,t1=13,t2=14,t3=13,t4=15)
ss <- makeSS(4,outdat,0.02,100,t1=0.2,t2=0.35,t3=0.5,t4=0.65)
results <- fitSS(100,4,outdat,ss)
mean(results[[1]]$p < 0.05)

#--------------------------------
#Old fitSS
fitSS <- function(nreps,t,dat,ss){
  test <- dat
  comp <- expand.grid(1:t,1:t) %>% filter(Var1 != Var2 , Var1 < Var2)
  names(comp) <- c("Grp1","Grp2")
  results <- list()
  test2 <- list()
  for(j in 1:nrow(comp)){
    test2[[j]] <- test[which(test$ttt %in% c(comp[j,]$Grp1,comp[j,]$Grp2)),]
    results[[j]] <- list()
    for(i in 1:nreps){
      test2[[j]]$resp <- ss[[j]][,i]
      res<-coef(summary(MASS::glmmPQL(resp ~ ttt ,
                                      random = ~ 1 | site, family = binomial,
                                      data=test2[[j]],verbose=F)))["ttt", ]
      names(res)<-c("est","se","df","t","p")
      results[[j]][[i]] <- res
    }
    results[[j]] <- bind_rows(results[[j]])
  }
  return(results)
}
#results <- fitSS(100,3,outdat,ss)

