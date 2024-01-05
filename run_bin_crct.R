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


#loop the data generation
outdat <- list()
for(j in 1:nrow(properties)){
  outdat[[j]] <- makeClusters(t=4,nid=properties$n_per_k[j],t1=properties$k[j],t2=properties$k[j],t3=properties$k[j],t4=properties$k[j])
}
j <- 1
#set how many workers you want to use
plan(multisession,workers=20)

if (Sys.getenv('HOME') == '/home/carlos') {
  #Home computer
  baepath <- "D:/Programs/PhDProject2/Programs/adapt_arm.stan"
  set_cmdstan_path(path="C:/Users/nolan/Documents/.cmdstan/cmdstan-2.33.1")
  outdir <- "E:/Simulations"
  
} else {
  #Work computer
  myPgm = '/home/criveros/.local/bin/helloWorld.exe'
}

#run the code
#test <- list()
tic()
mod <- cmdstan_model(baepath, pedantic = F, compile=T)
#test <- readRDS(here("simdata2.RDS"))
#THIS IS LEFT ON WHAT I RAN LAST
for(j in 58:60){
  test[[length(test)+1]] <- future_replicate(2500,testss(expdat=outdat[[j]],t=4,mod=mod,outdir=outdir,
                                           rho=properties$icc[j],t1=properties$t1[j],t2=properties$t2[j],t3=properties$t3[j],t4=properties$t4[j]),
                                 future.seed = 42L)
}
toc(quiet=FALSE)
#saveRDS(test,here("Data","nonadapt_sim.RDS"))
#Reframe the output for use
tempd <- test
for(j in c(1:60)){
  for(i in 1:2500){
    tempd[[j]][[i]]$sim <- i
    tempd[[j]][[i]]$property <- j
  }
}

outsim <- bind_rows(tempd)

properties2 <- properties %>% mutate(row = row_number()) 
outsim2 <- merge(outsim,properties2,by.y=c("row"),by.x="property")
#saveRDS(outsim2,file=here("testsim.RDS"))
