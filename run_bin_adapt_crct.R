#AUTHOR: Erin Nolan
#TITLE: PROJECT 2 ADAPTIVE SIMULATION
#PURPOSE: RUNS THE SIMULATIONS FOR MULTIARM cRCT UNDER VARIOUS PROPERTIES

#Run the functions defined in make_clusters and fit_bae
pacman::p_load(future,here,future.apply,tictoc,car,ggforce,rsimsum,dplyr,cmdstanr,rstan)

source(here("Programs","make_clusters.R"))
source(here("Programs","testFull.R"))
source(here("Programs","testInterim.R"))
source(here("Programs","makeDecision.R"))
source(here("Programs","runSimTrial.R"))

#PROJECT 2 
#The different trial properties
set.seed(580208819)
#The first interim is after either 3 or 5 clusters (out of 5 and 10)
properties <- expand.grid(trt_eff_scen = c(1,2,3), ctrl_prop = c(0.1), icc = c(0.05,0.2), n_per_k = c(5,25,50,75,100), k = c(5,10))

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
         interim = ifelse(k == 5, 3, 5))

#Put in the paths and options for the trial
plan(multisession,workers=20)
# point to openCL if you want GPU integration
#path_to_opencl_lib <- "C:/Program Files (x86)/OCL_SDK_Light/lib/x86_64"
#write(paste0("LDFLAGS= -L\"",path_to_opencl_lib,"\" -lOpenCL"), file.path(cmdstan_path(), "make", "local"), append = TRUE)
#
baepath <- "D:/Programs/PhDProject2/Programs/adapt_arm.stan"
set_cmdstan_path(path="C:/Users/nolan/Documents/.cmdstan/cmdstan-2.33.1")
outdir <- "J:/Sims"
#cpp_options = list(stan_opencl = T) add this to cmdstan_model for gpu intergration
mod <- cmdstan_model(baepath, pedantic = F, compile=T)
adaption <- "both" #this can be early_stopping, arm_dropping, or both
drop_cut <- 0.05
stop_cut <- 0.15
#this can be ties_prob or ties_opt. Ties opt drops the highest constraint arm if there is a tie for which has the lowest probability of success
#ties_prob drops the arm with the lowest pred prob estimate if there is a tie for which has the lowest probability of success
ties <- "ties_opt" 
#Run the trial
#test <- list()
for(j in 54){
  test[[length(test)+1]] <- future_replicate(2500,future.seed=42L,runSimTrial(properties,mod,outdir,j,adaption,drop_cut,stop_cut,ties))
}
#saveRDS(test,here("Data","adaptopt15_sim.RDS"))
#Take out the trial properties
trial_props <- list()
for(j in 1:36){
  for(i in seq(3,7500,3)){
    trial_props[[length(trial_props)+1]] <- test[[j]][[i]]
    trial_props[[length(trial_props)]]$sim <- i/3
    trial_props[[length(trial_props)]]$property <- j
  }
}
trial_props <- bind_rows(trial_props)
#saveRDS(trial_props,here("Data","adapt_trial_props.RDS"))

#Take out the interim analyses
interim <- list()
for(j in 1){
  for(i in seq(1,7500,3)){
    interim[[length(interim)+1]] <- test[[j]][[i]]
    interim[[length(interim)]]$sim <- (i+2)/3
    interim[[length(interim)]]$property <- j
  }
}
interim <- bind_rows(interim)
#saveRDS(interim,here("Data","adapt_interim.RDS"))

#Take out the full analyses
tempd <- list()
for(j in c(1:46)){
  for(i in seq(2,7500,3)){
    tempd[[length(tempd)+1]] <- test[[j]][[i]]
    tempd[[length(tempd)]]$sim <- (i+1)/3
    tempd[[length(tempd)]]$property <- j
  }
}
outsim <- bind_rows(tempd)

#merge in the properties of that simulation
properties2 <- properties %>% mutate(row = row_number()) 
outsim2 <- merge(outsim,properties2,by.y=c("row"),by.x="property")
#saveRDS(outsim2,here("Data","adapt_outsim.RDS"))

#Using https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4319656/pdf/nihms657495.pdf page 5 for adaptive early stopping so far