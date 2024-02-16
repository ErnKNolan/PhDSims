#AUTHOR: Erin Nolan
#TITLE: PROJECT 2 SENSITIVITY ANALYSIS
#PURPOSE: CONVERGENCE AND POWER SENSITIVITY

#Load in packages
pacman::p_load(here,future.apply,tictoc,car,ggforce,dplyr,cmdstanr,looplot,forcats,tidyr,gtsummary)

#identify the interim analyses that had poor convergence
#outsim2 <- readRDS(here("Data","interim_sim.RDS"))
int_sens <- outsim2 %>% 
  group_by(property,sim) %>%
  mutate(converged = max(ifelse(variable %in% c("beta_trt[2]", "beta_trt[3]","beta_trt[4]") & ess_bulk > 400,1,0))) %>%
  filter(converged == 1,variable=="beta_trt[4]")%>%
  dplyr::select(property,sim,converged)

#outsim2 <- readRDS(here("Data","outsim_adapt_opt.RDS"))
#outsim2 <- readRDS(here("Data","outsim_adapt_prob.RDS"))

#Merge interim with final to keep only the converged well interims
out_sens <- merge(outsim2,int_sens,by=c("property","sim"))

#FOR NON-ADAPT
#If non adapt then do this instead of interim
#outsim2 <- readRDS(here("Data","outsim_nonadapt.RDS"))
#out_sens <- outsim2 %>% 
#  group_by(property,sim) %>%
#  mutate(converged = max(ifelse(variable %in% c("beta_trt[2]", "beta_trt[3]","beta_trt[4]") & ess_bulk > 400,1,0))) %>%
#  filter(converged == 1)
  
#summarised values
converge_summ <- out_sens %>%
  group_by(property,sim) %>%
  filter(variable %in% c("beta_trt[2]", "beta_trt[3]","beta_trt[4]")) %>%
  mutate(bess_tail = ifelse(ess_tail < 400 | is.na(ess_tail),1,0),
         brhat = ifelse(rhat > 1.05 | is.na(rhat),1,0),
         sess_tail = ifelse(sum(bess_tail) >= 1,1,0),
         srhat = ifelse(sum(brhat) >= 1,1,0)) %>%
  filter(variable == "beta_trt[4]") %>%
  group_by(property) %>%
  summarise(m_ess_tail = mean(sess_tail == 1),
            mrhat = mean(srhat == 1))

#Calculate if a trial success
trial_success <- out_sens %>% 
  filter(variable %in% c("pp_trt2","pp_trt3","pp_trt4")) %>%
  pivot_wider(id_cols=c(property,sim,trt_eff_scen),names_from=variable,values_from=mean) %>%
  mutate(bayesr = ifelse(trt_eff_scen == 3 & (pp_trt2 >= 0.95 | pp_trt3 >= 0.95 | pp_trt4 >= 0.95),1,
                         ifelse(trt_eff_scen %in% c(1,2) & pp_trt4 >= 0.95,1,0))) 

#Calculate the power
power <- trial_success %>% group_by(property) %>%
  summarise(bayesr = sum(bayesr)/n())

#MC Error for power (ref Morris 2019 table 6)
MCSE_power <- data.frame(MCSE = sqrt((power$bayesr*(1-power$bayesr))/ 2500), property=power$property)

#Merge into ESS and rhat
convergence <- merge(converge_summ,MCSE_power,by="property") %>%
  merge(properties2, by.x="property",by.y="row")
  

#CONVERGENCE MEASURES
icc_perf <- convergence %>% group_by(icc) %>% 
  summarise(m_ess_tail = mean(m_ess_tail),
            mrhat = mean(mrhat),
            mMCSE = mean(MCSE)) %>%
  rename(property = icc) %>%
  mutate(property = case_when(property == 0.05 ~ "ICC = 0.05",
                              property == 0.2 ~ "ICC = 0.2"))

k_perf <- convergence %>% group_by(k) %>% 
  summarise(m_ess_tail = mean(m_ess_tail),
            mrhat = mean(mrhat),
            mMCSE = mean(MCSE)) %>%
  rename(property = k) %>%
  mutate(property = case_when(property == 5 ~ "k = 5",
                              property == 10 ~ "k = 10"))

n_perf <- convergence %>% group_by(n_per_k) %>% 
  summarise(m_ess_tail = mean(m_ess_tail),
            mrhat = mean(mrhat),
            mMCSE = mean(MCSE)) %>%
  rename(property = n_per_k) %>%
  mutate(property = case_when(property == 5 ~ "n = 5",
                              property == 25 ~ "n = 25",
                              property == 50 ~ "n = 50",
                              property == 75 ~ "n = 75",
                              property == 100 ~ "n = 100"))

trt_perf <- convergence %>% group_by(trt_eff_scen) %>% 
  summarise(m_ess_tail = mean(m_ess_tail),
            mrhat = mean(mrhat),
            mMCSE = mean(MCSE)) %>%
  rename(property = trt_eff_scen) %>%
  mutate(property = case_when(property == 1 ~ "Large effect",
                              property == 2 ~ "Mid/small effect",
                              property == 3 ~ "Null scenario"))

ov_perf <- convergence %>% 
  #group_by(trt_eff_scen) %>%
  summarise(m_ess_tail = mean(m_ess_tail),
            mrhat = mean(mrhat),
            mMCSE = mean(MCSE)) %>%
  mutate(property = "Overall") %>%
  dplyr::select(property,m_ess_tail,mrhat,mMCSE)

#COMBINE THE PERFORMANCE MEASURES
performance <- rbind(ov_perf,icc_perf,k_perf,n_perf,trt_perf) %>%
  mutate(m_ess_tail = m_ess_tail*100,
         mrhat = mrhat*100)

tbl_summary(performance,by="property",
            statistic = list(all_continuous() ~ "{mean}"))

#%>% #The remaining code will save the output but needs Chrome. 
#as_gt() %>%
#gt::gtsave(filename = here("Output","performance_tab_prob15.png"))


#PERFORMANCE MEASURES



#Power for these sims-----------------------------------------------
#Power defined as number of sims that the lower CrI of trt 1 greater than all other groups upper CrI
trial_success <- out_sens %>% 
  filter(variable %in% c("pp_trt2","pp_trt3","pp_trt4")) %>%
  pivot_wider(id_cols=c(property,sim,trt_eff_scen),names_from=variable,values_from=mean) %>%
  mutate(bayesr = ifelse(trt_eff_scen == 3 & (pp_trt2 >= 0.95 | pp_trt3 >= 0.95 | pp_trt4 >= 0.95),1,
                         ifelse(trt_eff_scen %in% c(1,2) & pp_trt4 >= 0.95,1,0))) 
power <- trial_success %>% group_by(property) %>%
  summarise(bayesr = sum(bayesr)/n())

outsim3 <- merge(outsim2,power,by="property") %>%
  group_by(property) %>%
  filter(sim == 1, variable == "pred_prob_trt[4]")

#saveRDS(outsim3,here("Data","power_plot_prob15_sens.RDS"))

#plot the power over different features
#setting up the data


power_plot_opt <- readRDS(here("Data","power_plot_opt15_sens.RDS"))
power_plot_prob <- readRDS(here("Data","power_plot_prob15_sens.RDS"))
nonadapt_power <- readRDS(here("Data","nonadapt_power_sens.RDS"))

nonadapt_plot <- nonadapt_power %>% ungroup() %>% filter(k %in% c(5,10), ctrl_prop == 0.1) %>%
  dplyr::select(trt_eff_scen,icc,n_per_k,k,bayesr) 
nest_plot <- power_plot_opt %>% ungroup() %>% dplyr::select(trt_eff_scen,icc,n_per_k,k,bayesr)
nest_plot_prob <- power_plot_prob %>% ungroup() %>% dplyr::select(trt_eff_scen,icc,n_per_k,k,bayesr)

#combine the two datasets
loop_plot <- inner_join(nest_plot,nonadapt_plot,by=c("icc","trt_eff_scen","n_per_k","k")) %>%
  rename(nonadapt = bayesr.y,
         both = bayesr.x) %>%
  inner_join(nest_plot_prob) %>%
  rename(nonadapt_prob = bayesr) %>%
  mutate(k = fct_rev(factor(k))) %>%
  rename(`Early stopping\narm dropping opt` = both,
         `Non-adpative` = nonadapt,
         `Early stopping\narm dropping prob` = nonadapt_prob) 
loop_plot_effs <- loop_plot %>% filter(trt_eff_scen != 3)
loop_plot_null <- loop_plot %>% filter(trt_eff_scen ==3) %>% dplyr::select(-trt_eff_scen)

#The graph
png(filename=here("Output","nest_loop_sens.png"),width=10,height=6,res=300,units="in")
nested_loop_plot(resdf = loop_plot_effs, 
                 x = "n_per_k", steps = c("icc","k", "trt_eff_scen"),
                 steps_y_base = -0.1, steps_y_height = 0.1, steps_y_shift = 0.1,
                 x_name = "Sample size per cluster", y_name = "Power",
                 spu_x_shift = 50,
                 line_alpha = 0.6,
                 point_alpha = 0.6,
                 steps_values_annotate = TRUE, steps_annotation_size = 2.5, 
                 hline_intercept = c(0,0.8,0.9,1), 
                 post_processing = list(
                   add_custom_theme = list(
                     axis.text.x = element_text(angle = -90, 
                                                vjust = 0.5, 
                                                size = 5))))+
  scale_colour_manual(values=c("#8D4585","#003151","orange"))+
  labs(color="Adaptive design",shape="Adaptive design",linetype="Adaptive design",size="Adaptive design")
dev.off()


#Graph for the null scenario
png(filename=here("Output","nest_loop_null_sens.png"),width=10,height=6,res=300,units="in")
nested_loop_plot(resdf = loop_plot_null, 
                 x = "n_per_k", steps = c("icc","k"),
                 steps_y_base = -0.025, steps_y_height = 0.025, steps_y_shift = 0.025,
                 x_name = "Sample size per cluster", y_name = "Type 1 error",
                 spu_x_shift = 50,
                 line_alpha = 0.6,
                 point_alpha = 0.6,
                 steps_values_annotate = TRUE, steps_annotation_size = 2.5, 
                 hline_intercept = c(0), 
                 post_processing = list(
                   add_custom_theme = list(
                     axis.text.x = element_text(angle = -90, 
                                                vjust = 0.5, 
                                                size = 5))))+
  scale_colour_manual(values=c("#8D4585","#003151","orange"))+
  labs(color="Adaptive design",shape="Adaptive design",linetype="Adaptive design",size="Adaptive design")
dev.off()