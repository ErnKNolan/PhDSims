#AUTHOR: Erin Nolan
#TITLE: PROJECT 2 SIMULATION PERFORMANCE
#PURPOSE: PLOTTING THE POWER AND ERROR OF THE SIMULATIONS
# install.packages("devtools")
#devtools::install_github("matherealize/looplot")

pacman::p_load(here,future.apply,tictoc,car,ggforce,dplyr,cmdstanr,looplot,forcats,tidyr)

#read in the data
#outsim2 <- readRDS(here("Data","adaptopt_outsim.RDS"))
#outsim2 <- readRDS(here("Data","adaptprob_outsim.RDS"))
#outsim2 <- readRDS(here("Data","outsim_nonadapt.RDS"))
#properties2 <- readRDS(here("Data","properties2.RDS"))

#sensitivity datasets
#outsim2 <- readRDS(here("Data","adaptopt_sens2_outsim.RDS"))
#outsim2 <- readRDS(here("Data","adaptprob_sens1_outsim.RDS"))

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

#BACK INTO POWER
#merge outsim2 back into outsim
outsim3 <- merge(outsim2,power,by="property") %>%
  group_by(property) %>%
  filter(sim == 1, variable == "pred_prob_trt[4]")

#plot the power over different features
#setting up the data
#power_plot_opt <- outsim3 %>%
  mutate(join = paste0("n per k=",n_per_k,", k=",k),
         sample_n = n_per_k*k, #total sample size
         deff = 1+icc*((n_per_k)-1), #design effect
         eff_n = sample_n / deff, #effective sample size
         iccf = paste0("ICC = ",factor(icc)),
         ctrlpf = factor(ctrl_prop),
         test = paste0(icc,ctrl_prop),
         kf = factor(k))

#save the power plot
#saveRDS(power_plot_opt,here("Data","power_plot_opt_sens2.RDS"))

#graph comparing adaptive and non-adaptive
#cleaning the datasets
#sensitivity datasets
#power_plot_opts <- readRDS(here("Data","power_plot_opt_sens2.RDS"))
#power_plot_probs <- readRDS(here("Data","power_plot_prob_sens1.RDS"))
#nest_plot <- power_plot_opts %>% ungroup() %>% dplyr::select(trt_eff_scen,icc,n_per_k,k,bayesr)
#nest_plot_prob <- power_plot_probs %>% ungroup() %>% dplyr::select(trt_eff_scen,icc,n_per_k,k,bayesr)


#main datasets
#power_plot_opt <- readRDS(here("Data","power_plot_opt.RDS"))
power_plot_prob <- readRDS(here("Data","power_plot_prob.RDS"))
nonadapt_power <- readRDS(here("Data","nonadapt_power.RDS"))
nonadapt_plot <- nonadapt_power %>% ungroup() %>% filter(k %in% c(5,10), ctrl_prop == 0.1) %>%
  dplyr::select(trt_eff_scen,icc,n_per_k,k,bayesr) 
#nest_plot <- power_plot_opt %>% ungroup() %>% dplyr::select(trt_eff_scen,icc,n_per_k,k,bayesr)
nest_plot_prob <- power_plot_prob %>% ungroup() %>% dplyr::select(trt_eff_scen,icc,n_per_k,k,bayesr)

#combine the two datasets
loop_plot <- inner_join(nest_plot_prob,nonadapt_plot,by=c("icc","trt_eff_scen","n_per_k","k")) %>%
  rename(nonadapt = bayesr.y,
         prob = bayesr.x) %>%
  mutate(k = fct_rev(factor(k))) %>%
  rename(`Non-adpative` = nonadapt,
         `Adaptive` = prob,
         ICC = icc,
         Scenario = trt_eff_scen) %>%
  mutate(Scenario = case_when(Scenario == 1 ~ "Strong effect",
                            Scenario == 2 ~ "Mid/Weak effect",
                            Scenario == 3 ~ "Null effect"))
loop_plot$Scenario <- factor(loop_plot$Scenario, levels = c("Strong effect","Mid/Weak effect","Null effect"))
loop_plot_effs <- loop_plot %>% filter(Scenario != "Null effect")
loop_plot_null <- loop_plot %>% filter(Scenario =="Null effect") %>% dplyr::select(-Scenario)

#The graph
png(filename=here("Output","nest_loop.png"),width=10,height=6,res=300,units="in")
nested_loop_plot(resdf = loop_plot_effs, 
                     x = "n_per_k", steps = c("ICC","k", "Scenario"),
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
  scale_colour_manual(values=c("#6d6d6d","black","black"))+
  labs(color="Design",shape="Design",linetype="Design",size="Design")
dev.off()

#Graph for the null scenario
png(filename=here("Output","nest_loop_null.png"),width=10,height=6,res=300,units="in")
nested_loop_plot(resdf = loop_plot_null, 
                 x = "n_per_k", steps = c("ICC","k"),
                 steps_y_base = -0.025, steps_y_height = 0.025, steps_y_shift = 0.025,
                 x_name = "Sample size per cluster", y_name = "Type 1 error",
                 spu_x_shift = 50,
                 line_alpha = 0.6,
                 point_alpha = 0.6,
                 steps_values_annotate = TRUE, steps_annotation_size = 2.5, 
                 hline_intercept = c(0), 
                 ylim = c(-0.1,0.25),
                 y_breaks = c(0,0.05,0.1,0.15,0.2,0.25),
                 post_processing = list(
                   add_custom_theme = list(
                     axis.text.x = element_text(angle = -90, 
                                                vjust = 0.5, 
                                                size = 5))))+
  scale_colour_manual(values=c("#6d6d6d","black","black"))+
  labs(color="Design",shape="Design",linetype="Design",size="Design")
dev.off()

#Trial properties
#trial_props <- readRDS(here("Data","adaptopt_trial_props.RDS"))
trial_props <- readRDS(here("Data","adaptprob_trial_props.RDS"))

stops <- trial_props %>% group_by(property) %>% summarise(stops = (sum(stop)/n()))
trt2drps <- trial_props %>% group_by(property) %>% summarise(trt2drps = (sum(drop=="trt2")/n()))
trt3drps <- trial_props %>% group_by(property) %>% summarise(trt3drps = (sum(drop=="trt3")/n()))
trt4drps <- trial_props %>% group_by(property) %>% summarise(trt4drps = (sum(drop=="trt4")/n()))
trial_drops <- trial_props %>% group_by(property) %>% filter(row_number() == 1)
trial_drops <- merge(trial_drops,stops,by="property") %>%
  merge(trt2drps,by="property") %>%
  merge(trt3drps,by="property") %>%
  merge(trt4drps,by="property") %>%
  dplyr::select(trt_eff_scen,icc,n_per_k,k,stops,trt2drps,trt3drps,trt4drps) %>%
  mutate(k = fct_rev(factor(k))) %>% 
  rename(`Stop for futility` = stops,
         `Arm 2 dropped` = trt2drps,
         `Arm 3 dropped` = trt3drps,
         `Arm 4 dropped` = trt4drps,
         ICC = icc,
         Scenario = trt_eff_scen) %>%
  mutate(Scenario = case_when(Scenario == 1 ~ "Strong effect",
                              Scenario == 2 ~ "Mid/Weak effect",
                              Scenario == 3 ~ "Null effect"))
trial_drops$Scenario <- factor(trial_drops$Scenario, levels = c("Strong effect","Mid/Weak effect","Null effect"))


png(filename=here("Output","trialprobprob_drops.png"),width=10,height=6,res=300,units="in")
nested_loop_plot(resdf = trial_drops, 
                 x = "n_per_k", steps = c("ICC","k","Scenario"),
                 steps_y_base = -0.1, steps_y_height = 0.1, steps_y_shift = 0.1,
                 x_name = "Sample size per cluster", y_name = "Proportion",
                 spu_x_shift = 50,
                 steps_values_annotate = TRUE, steps_annotation_size = 2.5, 
                 hline_intercept = c(0), 
                 ylim = c(-0.75,1),
                 y_breaks = c(0,0.25,0.5,0.75,1),
                 post_processing = list(
                   add_custom_theme = list(
                     axis.text.x = element_text(angle = -90, 
                                                vjust = 0.5, 
                                                size = 5))))+
  scale_colour_manual(values=c("#48157F","#29AF7F","#f7cb48f9","#a65c85ff"))+
  labs(color="Trial outcome",shape="Trial outcome",linetype="Trial outcome",size="Trial outcome")

dev.off()
