#AUTHOR: Erin Nolan
#TITLE: PROJECT 2 SIMULATION CONVERGENCE
#PURPOSE: PLOTTING THE CONVERGENCE OF THE SIMULATIONS
# install.packages("devtools")
#devtools::install_github("matherealize/looplot")

pacman::p_load(here,future.apply,tictoc,car,ggforce,dplyr,cmdstanr,looplot,forcats,tidyr,gtsummary)

#RUN THIS CODE UNTIL NEXT SECTION FOR EACH OUTSIM2

#outsim2 <- readRDS(here("Data","adaptprob_outsim.RDS"))
#outsim2 <- readRDS(here("Data","outsim_nonadapt.RDS"))
#outsim2 <- readRDS(here("Data","outsim_interim.RDS"))
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
  mutate(bess_bulk = ifelse(ess_bulk < 400 | is.na(ess_bulk),1,0), 
         bess_tail = ifelse(ess_tail < 400 | is.na(ess_tail),1,0), #ess in the tail
         brhat = ifelse(rhat > 1.05 | is.na(rhat),1,0), #whether rhat greater than 1.05
         sess_bulk = ifelse(sum(bess_bulk) >= 1,1,0), #proportion of ess bulk < 400
         sess_tail = ifelse(sum(bess_tail) >= 1,1,0), #proportion of ess tail < 400
         srhat = ifelse(sum(brhat) >= 1,1,0)) %>% #proportion of rhat > 1.05
  group_by(property) %>%
  filter(variable == "beta_trt[4]") %>% #keep just 1 row per simulation
  summarise(m_ess_bulk = mean(sess_bulk == 1), #summarise by each property
            m_ess_tail = mean(sess_tail == 1),
            mrhat = mean(srhat == 1))

#by treatment group
ESS_trt <- outsim2 %>% 
  filter(variable %in% c("beta_trt[2]","beta_trt[3]","beta_trt[4]")) %>%
  group_by(property,sim,variable) %>%
  mutate(bess_bulk = ifelse(ess_bulk < 400 | is.na(ess_bulk),1,0), 
         bess_tail = ifelse(ess_tail < 400 | is.na(ess_tail),1,0),
         brhat = ifelse(rhat > 1.05 | is.na(rhat),1,0),
         sess_bulk = ifelse(sum(bess_bulk) >= 1,1,0),
         sess_tail = ifelse(sum(bess_tail) >= 1,1,0),
         srhat = ifelse(sum(brhat) >= 1,1,0)) %>%
  group_by(property,variable) %>%
  summarise(m_ess_bulk = mean(sess_bulk == 1),
            m_ess_tail = mean(sess_tail == 1),
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
#getting the post prob of max(trt eff) taken from Stan code
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
#saveRDS(convergence,here("Data","converg_XXXXX.RDS"))

#THIS IS THE NEXT SECTION
#PLOTTING-----------------------------------------------------------------------------
#LOOP PLOT
plot_prob <- readRDS(here("Data","converg_prob.RDS"))
plot_nonadapt <- readRDS(here("Data","converg_nonadapt.RDS"))
plot_interim <- readRDS(here("Data","converg_interim.RDS"))
plot_nonadapt <- plot_nonadapt %>% filter(k %in% c(5,10), ctrl_prop == 0.1) %>%
  dplyr::select(m_ess_bulk,m_ess_tail,mrhat,trt_eff_scen,icc,n_per_k,k)
plot_prob <- plot_prob %>% dplyr::select(m_ess_bulk,m_ess_tail,mrhat,trt_eff_scen,icc,n_per_k,k)
plot_interim <- plot_interim %>% dplyr::select(m_ess_bulk,m_ess_tail,mrhat,trt_eff_scen,icc,n_per_k,interim) %>%
  rename(k = interim)

#combine the two datasets
loop_plot <- 
  plot_prob %>%
  rename(bulkprob = m_ess_bulk,
         tailprob = m_ess_tail,
         rhatprob = mrhat) %>%
  inner_join(plot_nonadapt) %>%
  rename(bulknonadapt = m_ess_bulk,
         bulk300nonadapt = m_ess_bulk300,
         tailnonadapt = m_ess_tail,
         rhatnonadapt = mrhat,
         ICC = icc,
         Scenario = trt_eff_scen) %>%
  mutate(k = fct_rev(factor(k)),
         Scenario = case_when(Scenario == 1 ~ "Strong effect",
                              Scenario == 2 ~ "Moderate effect",
                              Scenario == 3 ~ "Null effect"))
loop_plot$Scenario <- factor(loop_plot$Scenario, levels = c("Strong effect","Moderate effect","Null effect"))

#The graph
loop_plot_bulk <- loop_plot %>% dplyr::select(bulkprob,bulknonadapt,ICC,k,n_per_k,Scenario) %>%
  rename(`Adaptive` = bulkprob,
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
  scale_colour_manual(values=c("#6d6d6d","black","black")) +
  labs(color="Design",shape="Design",linetype="Design",size="Design")
dev.off()

#tail graph
loop_plot_tail <- loop_plot %>% dplyr::select(tailprob,tailnonadapt,ICC,k,n_per_k,Scenario) %>%
  rename(`Adaptive` = tailprob,
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
  scale_colour_manual(values=c("#6d6d6d","black","black")) +
  labs(color="Design",shape="Design",linetype="Design",size="Design")
dev.off()

#rhat
loop_plot_rhat <- loop_plot %>% dplyr::select(rhatprob,rhatnonadapt,ICC,k,n_per_k,Scenario) %>%
  rename(`Adaptive` = rhatprob,
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
  scale_colour_manual(values=c("#6d6d6d","black","black")) +
  labs(color="Design",shape="Design",linetype="Design",size="Design")
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
                              Scenario == 2 ~ "Moderate effect",
                              Scenario == 3 ~ "Null effect"))
plot_interim2$Scenario <- factor(plot_interim2$Scenario, levels = c("Strong effect","Moderate effect","Null effect"))

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

