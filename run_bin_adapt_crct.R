
#Run the functions defined in make_clusters and fit_bae
pacman::p_load(here,future.apply,tictoc,car,ggforce,rsimsum,dplyr,cmdstanr,rstan)

source(here("Programs","make_clusters.R"))
source(here("Programs","adapt_fit_bae.R"))

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

#loop the data generation
#NOTE FOR ADAPTIVE, t1 etc read in the n clusters for that INTERIM
outdat <- list()
for(j in 1:nrow(properties)){
  outdat[[j]] <- makeClusters(t=4,nid=properties$n_per_k[j],
                              t1=properties$k[j],t2=properties$k[j],
                              t3=properties$k[j],t4=properties$k[j])
}
j <- 1

baepath <- "D:/Programs/PhDProject2/Programs/adapt_arm.stan"
set_cmdstan_path(path="C:/Users/nolan/Documents/.cmdstan/cmdstan-2.33.1")
outdir <- "D:/Programs/Simulations"

#run the code
test <- list()
tic()
mod <- cmdstan_model(baepath, pedantic = F, compile=T)
#THIS IS LEFT ON WHAT I RAN LAST
for(j in 1:60){
  test[[length(test)+1]] <- future_replicate(3,testss(expdat=outdat[[j]],t=4,mod=mod,outdir=outdir,
                                                         rho=properties$icc[j],t1=properties$t1[j],t2=properties$t2[j],t3=properties$t3[j],t4=properties$t4[j]),
                                             future.seed = 42L)
}
toc(quiet=FALSE)

#save simulated responses
resp <- list()
for(j in c(1:60)){
    resp[[j]] <- test[[j]][seq(2,length(test[[j]]),2)]
}

#save the results
tempd <- list()
for(j in c(1:60)){
  tempd[[j]] <- test[[j]][seq(1,length(test[[j]]),2)]
}

for(j in c(1:60)){
  for(i in 1:3){
    tempd[[j]][[i]]$sim <- i
    tempd[[j]][[i]]$property <- j
  }
}
#take out of a list
outsim <- bind_rows(tempd)
#merge in the properties of that simulation
properties2 <- properties %>% mutate(row = row_number()) 
outsim2 <- merge(outsim,properties2,by.y=c("row"),by.x="property")

#Determining what (if any) treatment group to drop
int_drop <- outsim2 %>% 
  filter(variable %in% c("pp_trt1","pp_trt2","pp_trt3")) %>% 
  dplyr::select(property,sim,variable,mean,ess_bulk,ess_tail) %>%
  pivot_wider(id_cols = c(property,sim),names_from = variable, values_from = mean)

#finding the smallest pred prob
int_drop$drop <- apply(int_drop[,c(3:5)], 1, min, na.rm = TRUE)

#dropping only if pred prob less than 0.05
int_drop <- int_drop %>%
  mutate(drop = ifelse(drop > 0.05, NA, drop),
         droptrt = pmap_chr(list(pp_trt1, pp_trt2, pp_trt3, drop), 
                             ~ifelse(any(c(..1, ..2, ..3) %in% ..4),
                                     c("trt1", "trt2", "trt3")[c(..1, ..2, ..3) %in% ..4],
                                     "none")))

#TODO: need to determine what rule for if multiple trts have the same pred prob
#This occurs when more than 1 is 0 pred prob. Currently takes the first which is always trt1
#given optimisation, could make rule to drop the most expensive
#trt1 > trt2 > trt3 > trt4

#Now outdat with full numbers
#we can do something with the int_drop that inputs something to this function


#this needs to be sim based now
outdat <- list()
for(j in 1:nrow(properties)){
  outdat[[j]] <- makeClusters(t=4,nid=properties$n_per_k[j],
                              t1=properties$k[j],t2=properties$k[j],
                              t3=properties$k[j],t4=properties$k[j])
}
j <- 1

#add the already simulated data to the dataset
#ORRRR - because same seed can we just tell it to sim and the first interim amount will be the same?? check this
#how to tell simulate to not simulate any extra for dropped arms and add for other arms
#if drop 1 arm and therefore 2 clusters left over, give to highest performing arm and control, or the 2 non control arms?
#increased chance of control in an trial might be an ethical no


