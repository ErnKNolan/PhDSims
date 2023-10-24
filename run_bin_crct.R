#AUTHOR: Erin Nolan
#TITLE: PROJECT 2 NON-ADAPTIVE SIMULATION
#PURPOSE: RUNS THE SIMULATIONS FOR MULTIARM cRCT UNDER VARIOUS PROPERTIES

#Run the functions defined in make_clusters and fit_bae
pacman::p_load(here,future.apply,tictoc,car,ggforce,rsimsum,dplyr,cmdstanr,rstan)

source(here("Programs","make_clusters.R"))
source(here("Programs","fit_bae.R"))

#PROJECT 2 
#The different trial properties
set.seed(580208819)
properties <- expand.grid(trt_eff_scen = c(1,2,3), ctrl_prop = c(0.1,0.35), icc = c(0.05,0.2), n_per_k = c(5,75), k = c(5,20))
properties <- properties %>%
  mutate(k = ifelse(n_per_k == 75 & k == 20,10,k)) #instead of 20 clusters with 75 per, its 10 clusters

#4OCT Add in n=5 k=10
oct_prop <- expand.grid(trt_eff_scen = c(1,2,3), ctrl_prop = c(0.1,0.35), icc = c(0.05,0.2), n_per_k = 5, k = 10)

#5OCT Add in n=50 and n=100 for ctrl prop = 0.1
oct2_prop <- expand.grid(trt_eff_scen = c(1,2,3), ctrl_prop = c(0.1), icc = c(0.05,0.2), n_per_k = c(50,100), k = 10)

#6OCT Add in n=50 and n=100 for ctrl prop = 0.1 and k=5
oct3_prop <- expand.grid(trt_eff_scen = c(1,2,3), ctrl_prop = c(0.1), icc = c(0.05,0.2), n_per_k = c(50,100), k = 5)

#9OCT Add in n=25 for k=5 and k=10
oct4_prop <- expand.grid(trt_eff_scen = c(1,2,3), ctrl_prop = c(0.1), icc = c(0.05,0.2), n_per_k = c(25), k = c(5,10))

#bind to properties
properties <- rbind(properties,oct_prop,oct2_prop,oct3_prop,oct4_prop) %>%
  mutate(t1 = case_when(trt_eff_scen == 1 ~ ctrl_prop+0.5,
                        trt_eff_scen == 2 ~ ctrl_prop+0.4,
                        trt_eff_scen == 3 ~ ctrl_prop+0),
         t2 = case_when(trt_eff_scen == 1 ~ ctrl_prop+0.3,
                        trt_eff_scen == 2 ~ ctrl_prop+0.3,
                        trt_eff_scen == 3 ~ ctrl_prop+0),
         t3 = case_when(trt_eff_scen == 1 ~ ctrl_prop+0.1,
                        trt_eff_scen == 2 ~ ctrl_prop+0.2,
                        trt_eff_scen == 3 ~ ctrl_prop+0),
         t4 = ctrl_prop)

#loop the data generation
outdat <- list()
for(j in 1:nrow(properties)){
  outdat[[j]] <- makeClusters(t=4,nid=properties$n_per_k[j],t1=properties$k[j],t2=properties$k[j],t3=properties$k[j],t4=properties$k[j])
}
j <- 1
#set how many workers you want to use
plan(multisession,workers=20) 
rstan_options(auto_write = TRUE)
#Work comp (comment out the irrelevant one)
#Sys.setenv(Home="C:/Users/Enolan")
#baepath <- "C:/Users/ENolan/Simulations/runbae.stan"
#set_cmdstan_path("C:/Users/ENolan/.cmdstan/cmdstan-2.33.1")
#outdir <- 
#Home comp (comment out the irrelevant one)
#Sys.setenv(Home="D:/Programs")
baepath <- "D:/Programs/PhDProject2/Programs/runbae.stan"
set_cmdstan_path(path="C:/Users/nolan/Documents/.cmdstan/cmdstan-2.33.1")
#outdir <- "D:/Programs/Simulations"

#run the code
#test <- list()
tic()
mod <- cmdstan_model(baepath, pedantic = F, compile=T)
#test <- readRDS(here("simdata2.RDS"))
#THIS IS LEFT ON WHAT I RAN LAST
for(j in 85:96){
  test[[length(test)+1]] <- future_replicate(2500,testss(expdat=outdat[[j]],t=4,mod=mod,
                                           rho=properties$icc[j],t1=properties$t1[j],t2=properties$t2[j],t3=properties$t3[j],t4=properties$t4[j]),
                                 future.seed = 42L)
}
toc(quiet=FALSE)
#saveRDS(test,here("simdata.RDS"))
#Reframe the output for use
tempd <- test
for(j in c(1:96)){
  for(i in 1:2500){
    tempd[[j]][[i]]$sim <- i
    tempd[[j]][[i]]$property <- j
  }
}

outsim <- bind_rows(tempd)

properties2 <- properties %>% mutate(row = row_number()) 
outsim2 <- merge(outsim,properties2,by.y=c("row"),by.x="property")
#saveRDS(outsim2,file=here("testsim.RDS"))
