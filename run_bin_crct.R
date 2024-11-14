#AUTHOR: Erin Nolan
#TITLE: NON-ADAPTIVE SIMULATIONS FOR OPTIMISING IMPLEMENTATION STRATEGIES
#PURPOSE: RUNS THE SIMULATIONS FOR MULTIARM cRCT UNDER VARIOUS PROPERTIES

#Run the functions defined in make_clusters and fit_bae
pacman::p_load(here,future.apply,tictoc,car,ggforce,rsimsum,dplyr,cmdstanr,rstan)

source("make_clusters.R")
source("fit_bae.R")

#The different trial properties
set.seed(580208819)
#trt_eff_scen = the treatment scenario
#ctrl_prop = the baseline proportion of the event of interest
#icc = intra-class correlation
#n_per_k = number of participants per cluster
#k = number of clusters
properties <- expand.grid(trt_eff_scen = c(1,2,3), ctrl_prop = c(0.1), icc = c(0.05,0.2), n_per_k = c(5,25,50), k = c(5,10))

#bind to properties
properties <- rbind(properties) %>%
  mutate(t4 = case_when(trt_eff_scen == 1 ~ ctrl_prop+0.5,
                        trt_eff_scen == 2 ~ ctrl_prop+0.4,
                        trt_eff_scen == 3 ~ ctrl_prop+0),
         t3 = case_when(trt_eff_scen == 1 ~ ctrl_prop+0.3,
                        trt_eff_scen == 2 ~ ctrl_prop+0.3,
                        trt_eff_scen == 3 ~ ctrl_prop+0),
         t2 = case_when(trt_eff_scen == 1 ~ ctrl_prop+0.1,
                        trt_eff_scen == 2 ~ ctrl_prop+0.2,
                        trt_eff_scen == 3 ~ ctrl_prop+0),
         t1 = ctrl_prop)


#loop the data generation
outdat <- list()
for(j in 1:nrow(properties)){
  outdat[[j]] <- makeClusters(t=4,nid=properties$n_per_k[j],t1=properties$k[j],t2=properties$k[j],t3=properties$k[j],t4=properties$k[j])
}
j <- 1
#set how many workers you want to use
plan(multisession,workers=20) 

#Put in the paths for the simulation
baepath <- "adapt_arm2.stan" #the file path for the Stan model code
set_cmdstan_path(path="/root/.cmdstan/cmdstan-2.33.1") #where cmdstan is located
outdir <- "SimTrash" #where the stan files will be output
mod <- cmdstan_model(baepath, pedantic = F, compile=T)

#run the trial
test <- list() #test is the list of all the output from the simulated trials
for(j in 1:36){ #loops through properties 1 to 60
  test[[length(test)+1]] <- future_replicate(2500,testss(expdat=outdat[[j]],t=4,mod=mod,outdir=outdir,
                                           rho=properties$icc[j],t1=properties$t1[j],t2=properties$t2[j],t3=properties$t3[j],t4=properties$t4[j]),
                                 future.seed = 42L)
  saveRDS(test,"nonadapt.RDS")
}
#save the output for the non-adaptive trial simulations
#saveRDS(test,"nonadapt.RDS")

#Reframe the output for use
tempd <- test
for(j in c(1:36)){
  for(i in 1:2500){
    tempd[[j]][[i]]$sim <- i
    tempd[[j]][[i]]$property <- j
  }
}

outsim <- bind_rows(tempd)
properties2 <- properties %>% mutate(row = row_number()) 
outsim2 <- merge(outsim,properties2,by.y=c("row"),by.x="property")
#saveRDS(outsim2,file=here("Data","outsim_nonadapt.RDS"))

#USE THIS OUTPUT (WITH THE ADAPTIVE DESIGN OUTPUT) IN adapt_sim_convergence and adapt_sim_performance
