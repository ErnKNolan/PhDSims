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
  mutate(t1 = case_when(trt_eff_scen == 1 ~ ctrl_prop+0.5,
                        trt_eff_scen == 2 ~ ctrl_prop+0.4,
                        trt_eff_scen == 3 ~ ctrl_prop+0),
         t2 = case_when(trt_eff_scen == 1 ~ ctrl_prop+0.3,
                        trt_eff_scen == 2 ~ ctrl_prop+0.3,
                        trt_eff_scen == 3 ~ ctrl_prop+0),
         t3 = case_when(trt_eff_scen == 1 ~ ctrl_prop+0.1,
                        trt_eff_scen == 2 ~ ctrl_prop+0.2,
                        trt_eff_scen == 3 ~ ctrl_prop+0),
         t4 = ctrl_prop,
         interim = ifelse(k == 5, 3, 5))

#Put in the paths and options for the trial
plan(multisession,workers=20)
baepath <- "D:/Programs/PhDProject2/Programs/adapt_arm.stan"
set_cmdstan_path(path="C:/Users/nolan/Documents/.cmdstan/cmdstan-2.33.1")
outdir <- "D:/Programs/Simulations"
mod <- cmdstan_model(baepath, pedantic = F, compile=T)
adaption <- "both" #this can be early_stopping, arm_dropping, or both

#Run the trial
test <- list()
for(j in 1:60){
  test[[length(test)+1]] <- future_replicate(100,future.seed=42L,runSimTrial(properties,mod,outdir,j,adaption))
}
#saveRDS(test,here("adapt_100sim.RDS"))
#Take out the trial properties
trial_props <- list()
for(j in 1:60){
  for(i in seq(3,300,3)){
    trial_props[[length(trial_props)+1]] <- test[[j]][[i]]
    trial_props[[length(trial_props)]]$sim <- i/3
    trial_props[[length(trial_props)]]$property <- j
  }
}
trial_props <- bind_rows(trial_props)
#saveRDS(trial_props,"adapt_trial_props.RDS")

stops <- trial_props %>% group_by(property) %>% summarise(stops = (sum(stop)/n()))
trt1drps <- trial_props %>% group_by(property) %>% summarise(trt1drps = (sum(drop=="trt1")/n()))
trt2drps <- trial_props %>% group_by(property) %>% summarise(trt2drps = (sum(drop=="trt2")/n()))
trt3drps <- trial_props %>% group_by(property) %>% summarise(trt3drps = (sum(drop=="trt3")/n()))
trial_drops <- trial_props %>% group_by(property) %>% filter(row_number() == 1)
trial_drops <- merge(trial_drops,stops,by="property") %>%
  merge(trt1drps,by="property") %>%
  merge(trt2drps,by="property") %>%
  merge(trt3drps,by="property") %>%
  select(trt_eff_scen,icc,n_per_k,k,stops,trt1drps,trt2drps,trt3drps)

png(filename=here("Output","trial_drops.png"),width=10,height=6,res=300,units="in")
nested_loop_plot(resdf = trial_drops, 
                 x = "n_per_k", steps = c("icc","k","trt_eff_scen"),
                 steps_y_base = -0.1, steps_y_height = 0.1, steps_y_shift = 0.1,
                 x_name = "Sample size per cluster", y_name = "Proportion",
                 spu_x_shift = 50,
                 steps_values_annotate = TRUE, steps_annotation_size = 2.5, 
                 hline_intercept = c(0,0.05,0.8), 
                 post_processing = list(
                   add_custom_theme = list(
                     axis.text.x = element_text(angle = -90, 
                                                vjust = 0.5, 
                                                size = 5))))+
  scale_colour_manual(values=c("#48157F","#29AF7F","#f7cb48f9","#a65c85ff"))
dev.off()
#Take out the interim analyses
interim <- list()
for(j in 1:60){
  for(i in seq(1,300,3)){
    interim[[length(interim)+1]] <- test[[j]][[i]]
    interim[[length(interim)]]$sim <- (i+2)/3
    interim[[length(interim)]]$property <- j
  }
}
interim <- bind_rows(interim)
#saveRDS(interim,"adapt_interim.RDS")

#Take out the full analyses
tempd <- list()
for(j in c(1:60)){
  for(i in seq(2,300,3)){
    tempd[[length(tempd)+1]] <- test[[j]][[i]]
    tempd[[length(tempd)]]$sim <- (i+1)/3
    tempd[[length(tempd)]]$property <- j
  }
}
outsim <- bind_rows(tempd)

#merge in the properties of that simulation
properties2 <- properties %>% mutate(row = row_number()) 
outsim2 <- merge(outsim,properties2,by.y=c("row"),by.x="property")
#saveRDS(outsim2,"adapt_outsim.RDS")

#Using https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4319656/pdf/nihms657495.pdf page 5 for adaptive early stopping so far