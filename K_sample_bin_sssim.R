pacman::p_load(tidyverse,purrr,lme4,MASS,lsmeans)

# Make the generated dataset---------------------------
# makeClusters function
# t = number of arms
# nid = number of participants per arm
# specify the number of clusters for each t
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

#outdat <- makeClusters(t=4,nid=30,t1=10,t2=10,t3=10,t4=10)

# Making the outcome data function----------------------
# MakeSS function
# nreps = number of simulations
# t = number of treatment arms
# expdat = the generated dataset to simulate the results for
# rho = the intraclass correlation
# ... = the proportion of participants that have the event of interest 
# in each treatment arm, denoted as t1=, t2= etc
testss <- function(t,expdat,rho,...){
  
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
  #when modelling, factor variable of ttt
  sigma2 <- (pi ^ 2) / 3
  theta <- sqrt((rho*sigma2)/(1-rho))
  names(theta)<-c("site.(Intercept)")
  #fitting ss
  comp <- expand.grid(1:t,1:t) %>% filter(Var1 != Var2 , Var1 < Var2)
  names(comp) <- c("Grp1","Grp2")
  results <- vector()

    
    return <- tryCatch({

    expdat$resp <- suppressMessages(simulate( ~ factor(ttt) + (1|site), nsim = 1, family = binomial, 
                             newdata = expdat,newparams = list(beta=beta, theta=theta))) %>% as.matrix()
    
    fit1 <- lme4::glmer(resp ~ factor(ttt) + (1 | site),family = binomial, 
                        data=expdat,control=glmerControl(optimizer="bobyqa",
                                                         check.conv.singular = .makeCC(action = "ignore",  tol = 1e-4)))
    res <- data.frame(summary(lsmeans(fit1,pairwise~ttt,adjust="none"))$contrasts)
    res$type3 <- car::Anova(fit1)$"Pr(>Chisq)"
    res$warn <- paste0("0",summary(fit1)$optinfo$conv$lme4$messages)
    results <- list(res)
    },
    
    error=function(e) { message(conditionMessage(e)) 
      res <- data.frame(contrast=NA,estimate=NA,SE=NA,df=NA,z.ratio=NA,p.value=NA)
      res$type3 <- NA
      res$warn <- paste0("0")
      results <- list(res)
    })
    
    return(return)

}


#ssout <- testss(nreps=1000,t=4,expdat=outdat,rho=0.05,t1=0.35,t2=0.45,t3=0.55,t4=0.75)


# Select the result wanted, type 3 or the specific contrast-------------
# selectRes function
# ssfit = results from the models
# result = whether you want "type3" for the overall difference between treatment
# groups or "1 - 2", "1 - 3" etc for a comparison between two groups
selectRes <- function(ssfit,result){
ssout <- bind_rows(ssfit,.id="sim")
if(result=="type3"){
type3 <- ssout %>% 
  group_by(sim) %>% 
  filter(row_number()==1) %>% 
  ungroup() %>%
  dplyr::select(type3)
return(type3)
}
else{
onecont <- ssout %>%
  filter(contrast==result) %>%
  dplyr::select(contrast,estimate,p.value)
return(onecont)
}
}

#result <- selectRes(ssfit,"type3")

#AN example how to run the functions
#outdat <- makeClusters(t=5,nid=20,t1=15,t2=10,t3=8,t4=10,t5=14)
#ss <- makeSS(100,5,outdat,0.04,t1=0.28,t2=0.35,t3=0.5,t4=0.65,t5=0.65)
#ssfit <- fitSS(100,5,outdat,ss)
#result <- selectRes(ssfit,"type3")
