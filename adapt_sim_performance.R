#AUTHOR: Erin Nolan
#TITLE: PROJECT 2 NON-ADAPTIVE SIMULATION
#PURPOSE: RUNS THE SIMULATIONS FOR MULTIARM cRCT UNDER VARIOUS PROPERTIES
# install.packages("devtools")
#devtools::install_github("matherealize/looplot")

pacman::p_load(here,future.apply,tictoc,car,ggforce,dplyr,cmdstanr,looplot,forcats,tidyr)

#read in the data
#outsim2 <- readRDS(file=here("Data","adapt_outsim.RDS"))

#missings--------------------------------------------
missings <- outsim2 %>% 
  filter(variable != "pp_trt1",
         is.na(ess_tail) | is.na(ess_bulk))

#Convergence diagnostics-----------------------------
# ESS and rhat
ESS <- outsim2 %>% 
  filter(variable %in% c("beta_trt[2]","beta_trt[3]","beta_trt[4]")) %>%
  group_by(property,sim) %>%
  mutate(bess_bulk = ifelse(ess_bulk < 400 | is.na(ess_bulk),1,0), #currently treating missing ess as 'bad' but that isnt necessarily true
         bess_tail = ifelse(ess_tail < 400 | is.na(ess_tail),1,0),
         brhat = ifelse(rhat > 1.05 | is.na(rhat),1,0),
         sess_bulk = ifelse(sum(bess_bulk) >= 1,1,0),
         sess_tail = ifelse(sum(bess_tail) >= 1,1,0),
         srhat = ifelse(sum(brhat) >= 1,1,0)) %>%
  group_by(property) %>%
  filter(variable == "beta_trt[4]") %>%
  summarise(m_ess_bulk = mean(sess_bulk == 1),
            m_ess_tail = mean(sess_tail == 1),
            mrhat = mean(srhat == 1))

# MC Error of each iteration for parameters of interest
#Has to be done after POWER

#Power for these sims-----------------------------------------------
#Power defined as number of sims that the lower CrI of trt 1 greater than all other groups upper CrI
trial_success <- outsim2 %>% 
  filter(variable %in% c("pp_trt2","pp_trt3","pp_trt4")) %>%
  pivot_wider(id_cols=c(property,sim,trt_eff_scen),names_from=variable,values_from=mean) %>%
  mutate(bayesr = ifelse(trt_eff_scen == 3 & (pp_trt2 >= 0.95 | pp_trt3 >= 0.95 | pp_trt4 >= 0.95),1,
                         ifelse(trt_eff_scen %in% c(1,2) & pp_trt4 >= 0.95,1,0))) 
power <- trial_success %>% group_by(property) %>%
  summarise(bayesr = sum(bayesr)/n())

#MC Error for power (ref Morris 2019 table 6)
MCSE_power <- data.frame(MCSE = sqrt((power$bayesr*(1-power$bayesr))/ 2500), property=power$property)
#Merge into ESS and rhat
convergence <- merge(ESS,MCSE_power,by="property") %>%
  merge(properties2,by.x="property",by.y="row")

#BACK INTO POWER
#merge outsim2 back into outsim
outsim3 <- merge(outsim2,power,by="property") %>%
  group_by(property) %>%
  filter(sim == 1, variable == "pred_prob_trt[4]")

#plot the power over different features
#setting up the data
power_plot_prob <- outsim3 %>%
  mutate(join = paste0("n per k=",n_per_k,", k=",k),
         sample_n = n_per_k*k, #total sample size
         deff = 1+icc*((n_per_k)-1), #design effect
         eff_n = sample_n / deff, #effective sample size
         iccf = paste0("ICC = ",factor(icc)),
         ctrlpf = factor(ctrl_prop),
         test = paste0(icc,ctrl_prop),
         kf = factor(k))

#save the power plot
#saveRDS(power_plot_prob,here("Data","power_plot_prob.RDS"))

#graph comparing adaptive and non-adaptive
#cleaning the two datasets
power_plot <- readRDS(here("Data","power_plot.RDS"))
power_plot_prob <- readRDS(here("Data","power_plot_prob.RDS"))
nonadapt_power <- readRDS(here("Data","nonadapt_power.RDS"))
nonadapt_plot <- nonadapt_power %>% ungroup() %>% filter(k %in% c(5,10), ctrl_prop == 0.1) %>%
  dplyr::select(trt_eff_scen,icc,n_per_k,k,bayesr) 
nest_plot <- power_plot %>% ungroup() %>% dplyr::select(trt_eff_scen,icc,n_per_k,k,bayesr)
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
png(filename=here("Output","nest_loop.png"),width=10,height=6,res=300,units="in")
nested_loop_plot(resdf = loop_plot_effs, 
                     x = "n_per_k", steps = c("icc","k", "trt_eff_scen"),
                     steps_y_base = -0.1, steps_y_height = 0.1, steps_y_shift = 0.1,
                     x_name = "Sample size per cluster", y_name = "Power",
                     spu_x_shift = 50,
                     line_alpha = 0.6,
                     point_alpha = 0.6,
                     steps_values_annotate = TRUE, steps_annotation_size = 2.5, 
                     hline_intercept = c(0,0.05,0.8,0.9,1), 
                     post_processing = list(
                       add_custom_theme = list(
                         axis.text.x = element_text(angle = -90, 
                                                    vjust = 0.5, 
                                                    size = 5))))+
  scale_colour_manual(values=c("#8D4585","#003151","orange"))+
  labs(color="Adaptive design",shape="Adaptive design",linetype="Adaptive design",size="Adaptive design")
dev.off()

#Graph for the null scenario
png(filename=here("Output","nest_loop_null.png"),width=10,height=6,res=300,units="in")
nested_loop_plot(resdf = loop_plot_null, 
                 x = "n_per_k", steps = c("icc","k"),
                 steps_y_base = -0.025, steps_y_height = 0.025, steps_y_shift = 0.025,
                 x_name = "Sample size per cluster", y_name = "Power",
                 spu_x_shift = 50,
                 line_alpha = 0.6,
                 point_alpha = 0.6,
                 steps_values_annotate = TRUE, steps_annotation_size = 2.5, 
                 hline_intercept = c(0,0.01,0.05,0.1), 
                 post_processing = list(
                   add_custom_theme = list(
                     axis.text.x = element_text(angle = -90, 
                                                vjust = 0.5, 
                                                size = 5))))+
  scale_colour_manual(values=c("#8D4585","#003151","orange"))+
  labs(color="Adaptive design",shape="Adaptive design",linetype="Adaptive design",size="Adaptive design")
dev.off()

#Trial properties
#trials_props <- readRDS(here("Data","adapt_trial_props.RDS"))

stops <- trial_props %>% group_by(property) %>% summarise(stops = (sum(stop)/n()))
trt1drps <- trial_props %>% group_by(property) %>% summarise(trt1drps = (sum(drop=="trt1")/n()))
trt2drps <- trial_props %>% group_by(property) %>% summarise(trt2drps = (sum(drop=="trt2")/n()))
trt3drps <- trial_props %>% group_by(property) %>% summarise(trt3drps = (sum(drop=="trt3")/n()))
trial_drops <- trial_props %>% group_by(property) %>% filter(row_number() == 1)
trial_drops <- merge(trial_drops,stops,by="property") %>%
  merge(trt1drps,by="property") %>%
  merge(trt2drps,by="property") %>%
  merge(trt3drps,by="property") %>%
  dplyr::select(trt_eff_scen,icc,n_per_k,k,stops,trt1drps,trt2drps,trt3drps) %>%
  mutate(k = fct_rev(factor(k))) %>% 
  rename(`Stop for futility` = stops,
         `Arm 1 dropped` = trt1drps,
         `Arm 2 dropped` = trt2drps,
         `Arm 3 dropped` = trt3drps)

png(filename=here("Output","trialprob_drops.png"),width=10,height=6,res=300,units="in")
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
  scale_colour_manual(values=c("#48157F","#29AF7F","#f7cb48f9","#a65c85ff"))+
  labs(color="Trial outcome",shape="Trial outcome",linetype="Trial outcome",size="Trial outcome")

dev.off()

#PERFORMANCE MEASURES
icc_perf <- convergence %>% group_by(trt_eff_scen,icc) %>% 
  summarise(m_ess_bulk = mean(m_ess_bulk),
            m_ess_tail = mean(m_ess_tail),
            mrhat = mean(mrhat),
            mMCSE = mean(MCSE)) %>%
  rename(property = icc) %>%
  mutate(property = case_when(property == 0.05 ~ "ICC = 0.05",
                              property == 0.2 ~ "ICC = 0.2"))

k_perf <- convergence %>% group_by(trt_eff_scen,k) %>% 
  summarise(m_ess_bulk = mean(m_ess_bulk),
            m_ess_tail = mean(m_ess_tail),
            mrhat = mean(mrhat),
            mMCSE = mean(MCSE)) %>%
  rename(property = k) %>%
  mutate(property = case_when(property == 5 ~ "k = 5",
                              property == 10 ~ "k = 10"))

n_perf <- convergence %>% group_by(trt_eff_scen,n_per_k) %>% 
  summarise(m_ess_bulk = mean(m_ess_bulk),
            m_ess_tail = mean(m_ess_tail),
            mrhat = mean(mrhat),
            mMCSE = mean(MCSE)) %>%
  rename(property = n_per_k) %>%
  mutate(property = case_when(property == 5 ~ "n = 5",
                              property == 25 ~ "n = 25",
                              property == 50 ~ "n = 50",
                              property == 75 ~ "n = 75",
                              property == 100 ~ "n = 100"))

ov_perf <- convergence %>% group_by(trt_eff_scen) %>%
  summarise(m_ess_bulk = mean(m_ess_bulk),
            m_ess_tail = mean(m_ess_tail),
            mrhat = mean(mrhat),
            mMCSE = mean(MCSE)) %>%
  mutate(property = "Overall") %>%
  dplyr::select(trt_eff_scen,property,m_ess_bulk,m_ess_tail,mrhat,mMCSE)

#COMBINE THE PERFORMANCE MEASURES
performance <- rbind(ov_perf,icc_perf,k_perf,n_perf) %>%
  mutate(m_ess_bulk = m_ess_bulk*100,
         m_ess_tail = m_ess_tail*100,
         mrhat = mrhat*100)
