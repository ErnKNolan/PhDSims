#AUTHOR: Erin Nolan
#TITLE: PROJECT 2 SIMULATION CONVERGENCE
#PURPOSE: PLOTTING THE CONVERGENCE OF THE SIMULATIONS
# install.packages("devtools")
#devtools::install_github("matherealize/looplot")

pacman::p_load(here,future.apply,tictoc,car,ggforce,dplyr,cmdstanr,looplot,forcats,tidyr,gtsummary)

#outsim2 <- readRDS(here("Data","adaptopt_outsim.RDS"))
#outsim2 <- readRDS(here("Data","adaptprob_outsim.RDS"))
#outsim2 <- readRDS(here("Data","outsim_nonadapt.RDS"))
outsim2 <- readRDS(here("Data","outsim_interim.RDS"))
properties2 <- readRDS(here("Data","properties2.RDS"))

#missings--------------------------------------------
missings2 <- outsim2 %>% 
  filter(variable %in% c("beta_trt[2]","beta_trt[3]","beta_trt[4]"),
         is.na(ess_tail) | is.na(ess_bulk))

#Convergence diagnostics-----------------------------
# ESS and rhat
ESS <- outsim2 %>% 
  filter(variable %in% c("beta_trt[2]","beta_trt[3]","beta_trt[4]")) %>%
  group_by(property,sim) %>%
  mutate(bess_bulk = ifelse(ess_bulk < 400 | is.na(ess_bulk),1,0), #currently treating missing ess as 'bad' but that isnt necessarily true
         bess_tail = ifelse(ess_tail < 400 | is.na(ess_tail),1,0),
         bess_bulk300 = ifelse(ess_bulk < 300 | is.na(ess_bulk),1,0),
         brhat = ifelse(rhat > 1.05 | is.na(rhat),1,0),
         sess_bulk = ifelse(sum(bess_bulk) >= 1,1,0),
         sess_tail = ifelse(sum(bess_tail) >= 1,1,0),
         sess_bulk300 = ifelse(sum(bess_bulk300) >= 1,1,0),
         srhat = ifelse(sum(brhat) >= 1,1,0)) %>%
  group_by(property) %>%
  filter(variable == "beta_trt[4]") %>%
  summarise(m_ess_bulk = mean(sess_bulk == 1),
            m_ess_tail = mean(sess_tail == 1),
            m_ess_bulk300 = mean(sess_bulk300 == 1),
            mrhat = mean(srhat == 1))

#by treatment group
ESS_trt <- outsim2 %>% 
  filter(variable %in% c("beta_trt[2]","beta_trt[3]","beta_trt[4]")) %>%
  group_by(property,sim,variable) %>%
  mutate(bess_bulk = ifelse(ess_bulk < 400 | is.na(ess_bulk),1,0), #currently treating missing ess as 'bad' but that isnt necessarily true
         bess_tail = ifelse(ess_tail < 400 | is.na(ess_tail),1,0),
         bess_bulk300 = ifelse(ess_bulk < 300 | is.na(ess_bulk),1,0),
         brhat = ifelse(rhat > 1.05 | is.na(rhat),1,0),
         sess_bulk = ifelse(sum(bess_bulk) >= 1,1,0),
         sess_tail = ifelse(sum(bess_tail) >= 1,1,0),
         sess_bulk300 = ifelse(sum(bess_bulk300) >= 1,1,0),
         srhat = ifelse(sum(brhat) >= 1,1,0)) %>%
  group_by(property,variable) %>%
  summarise(m_ess_bulk = mean(sess_bulk == 1),
            m_ess_tail = mean(sess_tail == 1),
            m_ess_bulk300 = mean(sess_bulk300 == 1),
            mrhat = mean(srhat == 1))

#Plot of ess
outsim2 %>% 
  filter(variable %in% c("beta_trt[2]","beta_trt[3]","beta_trt[4]")) %>% 
  ggplot(aes(x=ess_bulk,group=n_per_k)) + 
  geom_density(aes(fill=n_per_k),alpha=0.3,position="identity") +
  facet_wrap(~k)
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
#saveRDS(convergence,here("Data","converg_interim.RDS"))


#LOOP PLOT
plot_opt <- readRDS(here("Data","converg_opt.RDS"))
plot_prob <- readRDS(here("Data","converg_prob.RDS"))
plot_nonadapt <- readRDS(here("Data","converg_nonadapt.RDS"))
plot_interim <- readRDS(here("Data","converg_interim.RDS"))
plot_nonadapt <- plot_nonadapt %>% filter(k %in% c(5,10), ctrl_prop == 0.1) %>%
  dplyr::select(m_ess_bulk,m_ess_bulk300,m_ess_tail,mrhat,trt_eff_scen,icc,n_per_k,k)
plot_opt <- plot_opt %>% dplyr::select(m_ess_bulk,m_ess_bulk300,m_ess_tail,mrhat,trt_eff_scen,icc,n_per_k,k)
plot_prob <- plot_prob %>% dplyr::select(m_ess_bulk,m_ess_bulk300,m_ess_tail,mrhat,trt_eff_scen,icc,n_per_k,k)
plot_interim <- plot_interim %>% dplyr::select(m_ess_bulk,m_ess_bulk300,m_ess_tail,mrhat,trt_eff_scen,icc,n_per_k,interim) %>%
  rename(k = interim)

#combine the two datasets
loop_plot <- inner_join(plot_prob,plot_opt,by=c("icc","trt_eff_scen","n_per_k","k")) %>%
  rename(bulkopt = m_ess_bulk.y,
         bulkprob = m_ess_bulk.x,
         bulk300opt = m_ess_bulk300.y,
         bulk300prob = m_ess_bulk300.x,
         tailopt = m_ess_tail.y,
         tailprob = m_ess_tail.x,
         rhatopt = mrhat.y,
         rhatprob = mrhat.x) %>%
  inner_join(plot_nonadapt) %>%
  rename(bulknonadapt = m_ess_bulk,
         bulk300nonadapt = m_ess_bulk300,
         tailnonadapt = m_ess_tail,
         rhatnonadapt = mrhat,
         ICC = icc,
         Scenario = trt_eff_scen) %>%
#  inner_join(plot_interim) %>%
#  rename(bulkinterim = m_ess_bulk,
#         bulk300interim = m_ess_bulk300,
#         tailinterim = m_ess_tail,
#         rhatinterim = mrhat) %>%
  mutate(k = fct_rev(factor(k)),
         Scenario = case_when(Scenario == 1 ~ "Strong effect",
                              Scenario == 2 ~ "Mid/Weak effect",
                              Scenario == 3 ~ "Null effect"))
loop_plot$Scenario <- factor(loop_plot$Scenario, levels = c("Strong effect","Mid/Weak effect","Null effect"))

#The graph
loop_plot_bulk <- loop_plot %>% dplyr::select(bulkprob,bulkopt,bulknonadapt,ICC,k,n_per_k,Scenario) %>%
  rename(`Probability ties` = bulkprob,
         `Optimisation\nconstraint ties` = bulkopt,
         `Non-adaptive` = bulknonadapt)
png(filename=here("Output","converge_bulk.png"),width=10,height=6,res=300,units="in")
nested_loop_plot(resdf = loop_plot_bulk, 
                 x = "n_per_k", steps = c("ICC","k", "Scenario"),
                 steps_y_base = -0.1, steps_y_height = 0.1, steps_y_shift = 0.1,
                 x_name = "Sample size per cluster", y_name = "Proportion of trials where ESS bulk < 400",
                 spu_x_shift = 50,
                 line_alpha = 0.6,
                 point_alpha = 0.6,
                 steps_values_annotate = TRUE, steps_annotation_size = 2.5, 
                 hline_intercept = c(0,1), 
                 post_processing = list(
                   add_custom_theme = list(
                     axis.text.x = element_text(angle = -90, 
                                                vjust = 0.5, 
                                                size = 5)))) +
  scale_colour_manual(values=c("#8D4585","#003151","orange","black")) +
  labs(color="Adaptive design",shape="Adaptive design",linetype="Adaptive design",size="Adaptive design")
dev.off()

#bulk 300 graph
loop_plot_bulk300 <- loop_plot %>% dplyr::select(bulk300prob,bulk300opt,bulk300nonadapt,ICC,k,n_per_k,Scenario) %>%
  rename(`Probability ties` = bulk300prob,
         `Optimisation\nconstraint ties` = bulk300opt,
         `Non-adaptive` = bulk300nonadapt)
png(filename=here("Output","converge_bulk300.png"),width=10,height=6,res=300,units="in")
nested_loop_plot(resdf = loop_plot_bulk300, 
                 x = "n_per_k", steps = c("ICC","k","Scenario"),
                 steps_y_base = -0.05, steps_y_height = 0.05, steps_y_shift = 0.05,
                 x_name = "Sample size per cluster", y_name = "Proportion of trials where ESS bulk < 300",
                 spu_x_shift = 50,
                 line_alpha = 0.6,
                 point_alpha = 0.6,
                 steps_values_annotate = TRUE, steps_annotation_size = 2.5, 
                 hline_intercept = c(0,1), 
                 post_processing = list(
                   add_custom_theme = list(
                     axis.text.x = element_text(angle = -90, 
                                                vjust = 0.5, 
                                                size = 5)))) +
  scale_colour_manual(values=c("#8D4585","#003151","orange","black")) +
  labs(color="Adaptive design",shape="Adaptive design",linetype="Adaptive design",size="Adaptive design")
dev.off()

#tail graph
loop_plot_tail <- loop_plot %>% dplyr::select(tailprob,tailopt,tailnonadapt,ICC,k,n_per_k,Scenario) %>%
  rename(`Probability ties` = tailprob,
         `Optimisation\nconstraint ties` = tailopt,
         `Non-adaptive` = tailnonadapt)
png(filename=here("Output","converge_tail.png"),width=10,height=6,res=300,units="in")
nested_loop_plot(resdf = loop_plot_tail, 
                 x = "n_per_k", steps = c("ICC","k", "Scenario"),
                 steps_y_base = -0.025, steps_y_height = 0.025, steps_y_shift = 0.025,
                 x_name = "Sample size per cluster", y_name = "Proportion of trials where ESS tail < 400",
                 spu_x_shift = 50,
                 line_alpha = 0.6,
                 point_alpha = 0.6,
                 steps_values_annotate = TRUE, steps_annotation_size = 2.5, 
                 hline_intercept = c(0,1), 
                 post_processing = list(
                   add_custom_theme = list(
                     axis.text.x = element_text(angle = -90, 
                                                vjust = 0.5, 
                                                size = 5)))) +
  scale_colour_manual(values=c("#8D4585","#003151","orange","black")) +
  labs(color="Adaptive design",shape="Adaptive design",linetype="Adaptive design",size="Adaptive design")
dev.off()

#rhat
loop_plot_rhat <- loop_plot %>% dplyr::select(rhatprob,rhatopt,rhatnonadapt,ICC,k,n_per_k,Scenario) %>%
  rename(`Probability ties` = rhatprob,
         `Optimisation\nconstraint ties` = rhatopt,
         `Non-adaptive` = rhatnonadapt)
png(filename=here("Output","converge_rhat.png"),width=10,height=6,res=300,units="in")
nested_loop_plot(resdf = loop_plot_rhat, 
                 x = "n_per_k", steps = c("ICC","k", "Scenario"),
                 steps_y_base = -0.1, steps_y_height = 0.1, steps_y_shift = 0.1,
                 x_name = "Sample size per cluster", y_name = "Proportion of trials where Rhat > 1.05",
                 spu_x_shift = 50,
                 line_alpha = 0.6,
                 point_alpha = 0.6,
                 steps_values_annotate = TRUE, steps_annotation_size = 2.5, 
                 hline_intercept = c(0,1), 
                 post_processing = list(
                   add_custom_theme = list(
                     axis.text.x = element_text(angle = -90, 
                                                vjust = 0.5, 
                                                size = 5)))) +
  scale_colour_manual(values=c("#8D4585","#003151","orange","black")) +
  labs(color="Adaptive design",shape="Adaptive design",linetype="Adaptive design",size="Adaptive design")
dev.off()

#interim
plot_interim2 <- plot_interim %>%
  rename(`ESS bulk < 400` = m_ess_bulk,
         `ESS bulk < 300` = m_ess_bulk300,
         `ESS tail < 400` = m_ess_tail,
         `Rhat > 1.05` = mrhat,
         ICC = icc,
         Scenario = trt_eff_scen) %>%
  mutate(k = fct_rev(factor(k)),
         Scenario = case_when(Scenario == 1 ~ "Strong effect",
                              Scenario == 2 ~ "Mid/Weak effect",
                              Scenario == 3 ~ "Null effect"))
plot_interim2$Scenario <- factor(plot_interim2$Scenario, levels = c("Strong effect","Mid/Weak effect","Null effect"))

png(filename=here("Output","converge_interim.png"),width=10,height=6,res=300,units="in")
nested_loop_plot(resdf = plot_interim2, 
                 x = "n_per_k", steps = c("ICC","k", "Scenario"),
                 steps_y_base = -0.05, steps_y_height = 0.05, steps_y_shift = 0.05,
                 x_name = "Sample size per cluster", y_name = "Proportion",
                 spu_x_shift = 50,
                 line_alpha = 0.6,
                 point_alpha = 0.6,
                 steps_values_annotate = TRUE, steps_annotation_size = 2.5, 
                 hline_intercept = c(0,1), 
                 post_processing = list(
                   add_custom_theme = list(
                     axis.text.x = element_text(angle = -90, 
                                                vjust = 0.5, 
                                                size = 5)))) +
  scale_colour_manual(values=c("#8D4585","#003151","orange","black")) +
  labs(color="Adaptive design",shape="Adaptive design",linetype="Adaptive design",size="Adaptive design")
dev.off()



#by trt group
ess_trt_rhat <- ESS_trt %>% 
  pivot_wider(id_cols=property,values_from = mrhat, names_from = variable) %>% 
  merge(properties2, by.x="property",by.y = "row") %>%
  dplyr::select(trt_eff_scen,icc,k,n_per_k,`beta_trt[2]`,`beta_trt[3]`,`beta_trt[4]`) %>%
  rename(ICC = icc,
         Scenario = trt_eff_scen)
png(filename=here("Output","rhat_trt_nonadapt.png"),width=10,height=6,res=300,units="in")
nested_loop_plot(resdf = ess_trt_rhat, 
                 x = "n_per_k", steps = c("ICC","k", "Scenario"),
                 steps_y_base = -0.025, steps_y_height = 0.025, steps_y_shift = 0.025,
                 x_name = "Sample size per cluster", y_name = "Proportion of trials where Rhat > 1.05",
                 spu_x_shift = 50,
                 line_alpha = 0.6,
                 point_alpha = 0.6,
                 steps_values_annotate = TRUE, steps_annotation_size = 2.5, 
                 hline_intercept = c(0,1), 
                 ylim = c(-0.2,0.15),
                 y_breaks = c(0,0.05,0.1,0.15),
                 post_processing = list(
                   add_custom_theme = list(
                     axis.text.x = element_text(angle = -90, 
                                                vjust = 0.5, 
                                                size = 5)))) +
  scale_colour_manual(values=c("#8D4585","#003151","orange","black")) +
  labs(color="Treatment group",shape="Treatment group",linetype="Treatment group",size="Treatment group")
dev.off()

ess_trt_bulk <- ESS_trt %>% 
  pivot_wider(id_cols=property,values_from = m_ess_bulk, names_from = variable) %>% 
  merge(properties2, by.x="property",by.y = "row") %>%
  dplyr::select(trt_eff_scen,icc,k,n_per_k,`beta_trt[2]`,`beta_trt[3]`,`beta_trt[4]`)
png(filename=here("Output","essbulk_trt_nonadapt.png"),width=10,height=6,res=300,units="in")
nested_loop_plot(resdf = ess_trt_bulk, 
                 x = "n_per_k", steps = c("icc","k", "trt_eff_scen"),
                 steps_y_base = -0.05, steps_y_height = 0.05, steps_y_shift = 0.05,
                 x_name = "Sample size per cluster", y_name = "Proportion of trials where ESS bulk < 400",
                 spu_x_shift = 50,
                 line_alpha = 0.6,
                 point_alpha = 0.6,
                 steps_values_annotate = TRUE, steps_annotation_size = 2.5, 
                 hline_intercept = c(0,1), 
                 ylim = c(-0.35,0.6),
                 y_breaks = c(0,0.2,0.4,0.6),
                 post_processing = list(
                   add_custom_theme = list(
                     axis.text.x = element_text(angle = -90, 
                                                vjust = 0.5, 
                                                size = 5)))) +
  scale_colour_manual(values=c("#8D4585","#003151","orange","black")) +
  labs(color="Treatment group",shape="Treatment group",linetype="Treatment group",size="Treatment group")
dev.off()

ess_trt_tail <- ESS_trt %>% 
  pivot_wider(id_cols=property,values_from = m_ess_tail, names_from = variable) %>% 
  merge(properties2, by.x="property",by.y = "row") %>%
  dplyr::select(trt_eff_scen,icc,k,n_per_k,`beta_trt[2]`,`beta_trt[3]`,`beta_trt[4]`)
png(filename=here("Output","esstail_trt_nonadapt.png"),width=10,height=6,res=300,units="in")
nested_loop_plot(resdf = ess_trt_tail, 
                 x = "n_per_k", steps = c("icc","k", "trt_eff_scen"),
                 steps_y_base = -0.025, steps_y_height = 0.025, steps_y_shift = 0.025,
                 x_name = "Sample size per cluster", y_name = "Proportion of trials where ESS tail < 400",
                 spu_x_shift = 50,
                 line_alpha = 0.6,
                 point_alpha = 0.6,
                 steps_values_annotate = TRUE, steps_annotation_size = 2.5, 
                 hline_intercept = c(0,1), 
                 ylim = c(-0.2,0.15),
                 y_breaks = c(0,0.05,0.1,0.15),
                 post_processing = list(
                   add_custom_theme = list(
                     axis.text.x = element_text(angle = -90, 
                                                vjust = 0.5, 
                                                size = 5)))) +
  scale_colour_manual(values=c("#8D4585","#003151","orange","black")) +
  labs(color="Treatment group",shape="Treatment group",linetype="Treatment group",size="Treatment group")
dev.off()