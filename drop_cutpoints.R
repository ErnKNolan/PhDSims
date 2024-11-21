
#test multiple cutpoints for trials with 3 and 5 clusters at interim
#e.g. the non-adaptive results at 5 clusters are the same as the adaptive results at their 5 cluster interim
#read in data AFTER compiling the results for non-adaptive designs at k clusters max
pacman::p_load(here,dplyr,tidyr)

#read in the results at 3 and 5 clusters
interim3 <- readRDS("outsim_interim3.RDS")
interim5 <- readRDS("outsim_nonadapt.RDS")
#combine the data and get rid of the 10 cluster interim
outsim2 <- rbind(interim3,interim5) %>% filter(k != 10)

#stopping rule decision for each simulated trial
stops <- outsim2 %>% filter(variable %in% c("pp_trt2","pp_trt3",'pp_trt4')) %>% 
  group_by(property,sim,trt_eff_scen,icc,n_per_k,k) %>% 
  dplyr::select(variable,mean) %>%
  pivot_wider(names_from="variable",values_from="mean") %>% #make the pp_trt data wide
  #trying 5 different stop cutpoints
  mutate(stop15 = ifelse(pp_trt2 < 0.15 & pp_trt3 < 0.15 & pp_trt4 < 0.15,1,0),
         stop2 = ifelse(pp_trt2 < 0.2 & pp_trt3 < 0.2 & pp_trt4 < 0.2,1,0),
         stop25 = ifelse(pp_trt2 < 0.25 & pp_trt3 < 0.25 & pp_trt4 < 0.25,1,0),
         stop35 = ifelse(pp_trt2 < 0.35 & pp_trt3 < 0.35 & pp_trt4 < 0.35,1,0),
         stop45 = ifelse(pp_trt2 < 0.45 & pp_trt3 < 0.45 & pp_trt4 < 0.45,1,0)) %>%
  group_by(property,trt_eff_scen,icc,n_per_k,k) %>% 
  #proportion of times the trials stopped for futility by property
  summarise(stop15 = sum(stop15)/n(),
            stop2 = sum(stop2)/n(),
            stop25 = sum(stop25)/n(),
            stop35 = sum(stop35)/n(),
            stop45 = sum(stop45)/n()) %>%
  mutate(trt_eff_scen = case_when(trt_eff_scen == 1 ~ "Strong effect",
                                  trt_eff_scen == 2 ~ "Moderate effect",
                                  trt_eff_scen == 3 ~ "No effect"))
#save the results
write.csv(stops,"stops.csv")


#Can also look at how different cutpoints for arm dropping would affect decisions
#k=5 aligns with the adaptive decisions specified in paper
armdrop <- outsim2 %>% filter(variable %in% c("pp_trt2","pp_trt3",'pp_trt4',
                                              "pred_prob_trt[2]","pred_prob_trt[3]","pred_prob_trt[4]")) %>% 
  group_by(property,sim,trt_eff_scen,icc,n_per_k,k) %>% 
  dplyr::select(variable,mean) %>%
  pivot_wider(names_from="variable",values_from="mean") %>% #make the pp_trt and pred prob wide
  #this is the cutpoint of 0.05 to drop an arm
  mutate(arm2drop_cut05 = ifelse(pp_trt2 <= 0.05,1,0),
         arm3drop_cut05 = ifelse(pp_trt3 <= 0.05,1,0),
         arm4drop_cut05 = ifelse(pp_trt4 <= 0.05,1,0),
         #determining if each arm has the minimum pp_trt if theres more than 1 pp_trt hitting dropping threshold
         arm_drop2 = ifelse(arm2drop_cut05 + arm3drop_cut05 + arm4drop_cut05 >= 1 & pp_trt2 == min(pp_trt2,pp_trt3,pp_trt4),1,0),
         arm_drop3 = ifelse(arm2drop_cut05 + arm3drop_cut05 + arm4drop_cut05 >= 1 & pp_trt3 == min(pp_trt2,pp_trt3,pp_trt4),1,0),
         arm_drop4 = ifelse(arm2drop_cut05 + arm3drop_cut05 + arm4drop_cut05 >= 1 & pp_trt4 == min(pp_trt2,pp_trt3,pp_trt4),1,0),
         #if there is a tie for pp_trt then judges from pred_prob_trt, if theres no tie then just checking it is the min pp_trt
         armdrop2_cut05 = ifelse((arm_drop2 + arm_drop3 + arm_drop4 > 1 & `pred_prob_trt[2]` == min(`pred_prob_trt[2]`,`pred_prob_trt[3]`,`pred_prob_trt[4]`)) | 
                              (arm_drop2 + arm_drop3 + arm_drop4 == 1 & arm_drop2 == 1),1,0),
         armdrop3_cut05 = ifelse((arm_drop2 + arm_drop3 + arm_drop4 > 1 & `pred_prob_trt[3]` == min(`pred_prob_trt[2]`,`pred_prob_trt[3]`,`pred_prob_trt[4]`)) | 
                              (arm_drop2 + arm_drop3 + arm_drop4 == 1 & arm_drop3 == 1),1,0),
         armdrop4_cut05 = ifelse((arm_drop2 + arm_drop3 + arm_drop4 > 1 & `pred_prob_trt[4]` == min(`pred_prob_trt[2]`,`pred_prob_trt[3]`,`pred_prob_trt[4]`)) | 
                              (arm_drop2 + arm_drop3 + arm_drop4 == 1 & arm_drop4 == 1),1,0)) %>%
  dplyr::select(-arm_drop2,-arm_drop3,-arm_drop4) %>%
  #same rules again but this time with a cutpoint of 0.1
  mutate(arm2drop_cut10 = ifelse(pp_trt2 <= 0.1,1,0),
         arm3drop_cut10 = ifelse(pp_trt3 <= 0.1,1,0),
         arm4drop_cut10 = ifelse(pp_trt4 <= 0.1,1,0),
         arm_drop2 = ifelse(arm2drop_cut05 + arm3drop_cut10 + arm4drop_cut10 >= 1 & pp_trt2 == min(pp_trt2,pp_trt3,pp_trt4),1,0),
         arm_drop3 = ifelse(arm2drop_cut05 + arm3drop_cut10 + arm4drop_cut10 >= 1 & pp_trt3 == min(pp_trt2,pp_trt3,pp_trt4),1,0),
         arm_drop4 = ifelse(arm2drop_cut05 + arm3drop_cut10 + arm4drop_cut10 >= 1 & pp_trt4 == min(pp_trt2,pp_trt3,pp_trt4),1,0),
         armdrop2_cut10 = ifelse((arm_drop2 + arm_drop3 + arm_drop4 > 1 & `pred_prob_trt[2]` == min(`pred_prob_trt[2]`,`pred_prob_trt[3]`,`pred_prob_trt[4]`)) | 
                                   (arm_drop2 + arm_drop3 + arm_drop4 == 1 & arm_drop2 == 1),1,0),
         armdrop3_cut10 = ifelse((arm_drop2 + arm_drop3 + arm_drop4 > 1 & `pred_prob_trt[3]` == min(`pred_prob_trt[2]`,`pred_prob_trt[3]`,`pred_prob_trt[4]`)) | 
                                   (arm_drop2 + arm_drop3 + arm_drop4 == 1 & arm_drop3 == 1),1,0),
         armdrop4_cut10 = ifelse((arm_drop2 + arm_drop3 + arm_drop4 > 1 & `pred_prob_trt[4]` == min(`pred_prob_trt[2]`,`pred_prob_trt[3]`,`pred_prob_trt[4]`)) | 
                                   (arm_drop2 + arm_drop3 + arm_drop4 == 1 & arm_drop4 == 1),1,0)) %>%
  dplyr::select(-arm_drop2,-arm_drop3,-arm_drop4,-arm2drop_cut05,-arm3drop_cut05,-arm4drop_cut05,
                -arm2drop_cut10,-arm3drop_cut10,-arm4drop_cut10)

#summarise the data to see proportion of trials that dropped which arm
droparm <- armdrop %>% group_by(property,trt_eff_scen,icc,n_per_k,k) %>% 
  summarise(armdrop2_cut05 = sum(armdrop2_cut05)/n(),
            armdrop3_cut05 = sum(armdrop3_cut05)/n(),
            armdrop4_cut05 = sum(armdrop4_cut05)/n(),
            armdrop2_cut10 = sum(armdrop2_cut10)/n(),
            armdrop3_cut10 = sum(armdrop3_cut10)/n(),
            armdrop4_cut10 = sum(armdrop4_cut10)/n())


