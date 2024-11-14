#AUTHOR: Erin Nolan
#TITLE: ADAPTIVE SIMULATION FOR OPTIMISING IMPLEMENTATION STRATEGIES
#PURPOSE: RUNS THE SIMULATIONS FOR MULTIARM cRCT UNDER VARIOUS PROPERTIES

#Run the functions defined and read in required packages
pacman::p_load(future,here,future.apply,tictoc,car,ggforce,rsimsum,dplyr,cmdstanr,rstan)
source("make_clusters.R")
source("testFull.R")
source("testInterim.R")
source("makeDecision.R")
source("runSimTrial.R")

#The different trial properties
set.seed(580208819) #main seed
#set.seed(698426789) #seed sensitivity 1
#set.seed(854187431) #seed sensitivity 2

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
         t1 = ctrl_prop,
         interim = ifelse(k == 5, 3, 5)) #The first interim is after either 3 or 5 clusters (out of 5 and 10)
saveRDS(properties,here("Data","properties.RDS"))

#Put in the paths and options for the trial
plan(multisession,workers=20) #running the model in parallel
baepath <- "adapt_arm.stan" #the file path for the Stan model code
set_cmdstan_path(path="/root/.cmdstan/cmdstan-2.33.1") #where cmdstan is located
outdir <- "SimTrash" #where the stan files will be output
mod <- cmdstan_model(baepath, pedantic = F, compile=T)
adaption <- "both" #this can be early_stopping, arm_dropping, or both
drop_cut <- 0.05 #cutpoint for dropping an arm
stop_cut <- 0.15 #cutpoint for stopping for futility

#this can be ties_prob or ties_opt. Ties opt drops the highest constraint arm if there is a tie for which has the lowest probability of success
#ties_prob drops the arm with the lowest pred prob estimate if there is a tie for which has the lowest probability of success
ties <- "ties_prob" 

#Run the trial
test <- list() #test is the list of all the output from the simulated trials
for(j in c(1:36)){
  test[[length(test)+1]] <- future_replicate(2500,future.seed=42L,runSimTrial(properties,mod,outdir,j,adaption,drop_cut,stop_cut,ties))
  saveRDS(test,"adaptprob.RDS")
}
#save the output for the adaptive trial simulations
#saveRDS(test,"adaptprob.RDS")

#NEXT PROGRAM TO USE IS save_adapt_simulations.R
