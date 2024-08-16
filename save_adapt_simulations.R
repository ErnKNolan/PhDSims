#AUTHOR: Erin Nolan
#TITLE: ADAPTIVE SIMULATION FOR OPTIMISING IMPLEMENTATION STRATEGIES
#PURPOSE: SAVES THE SIMULATED ADAPTIVE DATASET IN CHUNKS FOR DIFFERENT TESTS

#Run the functions defined and read in required packages
pacman::p_load(future,here,future.apply,tictoc,car,ggforce,rsimsum,dplyr,cmdstanr,rstan)

#read in the simulated dataset and properties
test <- readRDS(here("adaptprob.RDS"))
properties <- readRDS(here("Data","properties.RDS"))

#Take out the trial properties
trial_props <- list()
for(j in 1:60){
  for(i in seq(2,7500,2)){ #the middle number is the total number of items in the list
    trial_props[[length(trial_props)+1]] <- test[[j]][[i]]
    trial_props[[length(trial_props)]]$sim <- i/2 #label each simulation
    trial_props[[length(trial_props)]]$property <- j #label each property combination
  }
}
trial_props <- bind_rows(trial_props)
#saveRDS(trial_props,here("Data","adapt_trial_props.RDS"))

#Take out the interim analyses
interim <- list()
for(j in 1:60){
  for(i in seq(1,7500,3)){ #the middle number is the total number of items in the list
    interim[[length(interim)+1]] <- test[[j]][[i]]
    interim[[length(interim)]]$sim <- (i+2)/3 #label each simulation
    interim[[length(interim)]]$property <- j #label each property combination
  }
}
interim <- bind_rows(interim)
properties2 <- properties %>% mutate(row = row_number()) 
interim2 <- merge(interim,properties2,by.y=c("row"),by.x="property")
#saveRDS(interim2,here("Data","adapt_interim.RDS"))

#Take out the full analyses
tempd <- list()
for(j in c(1:60)){
  for(i in seq(1,7500,2)){ #the middle number is the total number of items in the list
    tempd[[length(tempd)+1]] <- test[[j]][[i]]
    tempd[[length(tempd)]]$sim <- (i+1)/2 #label each simulation
    tempd[[length(tempd)]]$property <- j #label each property combination
  }
}
outsim <- bind_rows(tempd)
#merge in the properties of that simulation
properties2 <- properties %>% mutate(row = row_number()) 
outsim2 <- merge(outsim,properties2,by.y=c("row"),by.x="property")
#saveRDS(properties2,here("Data","properties2.RDS"))
#saveRDS(outsim2,here("Data","adapt_outsim.RDS"))

#USE THIS OUTPUT (WITH THE NON-ADAPTIVE DESIGN OUTPUT) IN adapt_sim_convergence and adapt_sim_performance