
lme4::glmer(resp ~ ttt + (1 | site),family = binomial, 
            data=expdat) %>% summary()


expdat=outdat[[j]]
t=4
rho=properties$icc[j]
t1=properties$t1[j]
t2=properties$t2[j]
t3=properties$t3[j]
t4=properties$t4[j]
prop <- list(t1,t2,t3,t4)
comp <- (1:4)
o <- list()
beta <- vector()
for(i in 1:length(comp)){
  o[[i]] <- prop[[comp[i]]]/(1-prop[[comp[i]]])
}
beta <- append(beta,log(o[[1]]))
for(i in 2:length(comp)){
  beta <- append(beta,log(o[[i]]/o[[1]]))
}
#when modelling, factor variable of ttt
sigma2 <- (pi ^ 2) / 3
theta <- sqrt((rho*sigma2)/(1-rho))
names(theta)<-c("site.(Intercept)")
#fitting ss
comp <- expand.grid(1:t,1:t) %>% filter(Var1 != Var2 , Var1 < Var2)
names(comp) <- c("Grp1","Grp2")
expdat$ttt <- factor(expdat$ttt)
expdat$resp <- suppressMessages(simulate( ~ ttt + (1|site), nsim = 1, family = binomial, 
                                            newdata = expdat,newparams = list(beta=beta, theta=theta))) %>% as.matrix()
  
#run stan
sim <- expdat
resp <- as.vector(sim$resp)
N_obs <- dim(sim)[1]
N_site <- length(unique(sim$site))
N_ttt_groups <- length(unique(sim$ttt))
data <- list(N_obs = N_obs, N_site = N_site, N_ttt_groups = N_ttt_groups, 
             site = sim$site, ttt = as.numeric(sim$ttt), resp = resp)
# we recommend running this is a fresh R session or restarting your current session
#install.packages("cmdstanr", repos = c("https://mc-stan.org/r-packages/", getOption("repos")))
library(cmdstanr)

Sys.setenv(Home="C:/Users/Enolan")
path <- "C:/Users/ENolan/Simulations/runbae.stan"
set_cmdstan_path("C:/Users/ENolan/.cmdstan/cmdstan-2.32.2")
mod <- cmdstan_model(path, pedantic = F, compile=T)
fit <<- mod$sample(
  data = data, 
  init=0,
  iter_warmup = 500,
  iter_sampling = 250,
  chains = 4, 
  parallel_chains = 1,
  adapt_delta = 0.8,
  refresh = 50, # print update every 500 iters,
  max_treedepth=12
)
print(data.frame(fit$summary(variables=c("pred_prob_trt"))))
#$summary()[3:5,c(1,2,4,6,7,8,9,10)]

