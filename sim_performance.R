#AUTHOR: Erin Nolan
#TITLE: PROJECT 2 NON-ADAPTIVE SIMULATION
#PURPOSE: RUNS THE SIMULATIONS FOR MULTIARM cRCT UNDER VARIOUS PROPERTIES
pacman::p_load(here,future.apply,tictoc,car,ggforce,rsimsum,dplyr,cmdstanr)

#read in the data
outsim2 <- readRDS(file=here("testsim.RDS"))

missings <- outsim2 %>% 
  filter(variable != "diff",
         is.na(ess_tail) | is.na(ess_bulk))
#Convergence diagnostics
## ESS
ESS <- outsim2 %>% 
  filter(variable != "diff") %>%
  group_by(property) %>%
  summarise(m_ess_bulk = mean(ess_bulk,na.rm=TRUE),
            m_ess_tail = mean(ess_tail,na.rm=TRUE))
## potential-scale-down-reduction-factor R
scale_downR <- outsim2 %>%
  filter(variable != "diff") %>%
  group_by(property) %>%
  summarise(m_Rhat = mean(rhat,na.rm=TRUE))
## chain number, length, burn-ins
### 4 chains, per chain: 500 total, 250 burn-ins

## MC Error of each iteration for parameters of interest
###simsum
MCerror <- outsim2 %>%
  group_by(property,variable) %>%
  mutate(trt = case_when(variable == "pred_prob_trt[1]" ~ t1,
                         variable == "pred_prob_trt[2]" ~ t2,
                         variable == "pred_prob_trt[3]" ~ t3,
                         variable == "pred_prob_trt[4]" ~ t4)) %>%
  filter(variable == "diff") %>%
  simsum(estvarname="mean",true="trt",se="sd",method=c("property")) 
MCerror[["summ"]] %>% filter(stat == "power") %>% kable(MCerror)

#Power for these sims
#Power defined as number of sims that the lower CrI of trt 1 greater than all other groups upper CrI
power <- outsim2 %>% 
  filter(variable == "diff") %>%
  group_by(property) %>%
  mutate(bayes = ifelse(mean > 0.95,1,0)) %>%
  summarise(bayesr = sum(bayes)/n())

#merge back into outsim
outsim3 <- merge(outsim2,power,by="property") %>%
  group_by(property) %>%
  filter(sim == 1, variable == "pred_prob_trt[1]")

#checking that the model output aligns with the true values
#sums <- outsim2 %>% group_by(property,variable) %>% summarise(mean=mean(mean)) %>% pivot_wider(id_cols=property,names_from=variable,values_from=mean) %>% cbind(properties)

#plot the power over different features
#setting up the data
power_plot <- outsim3 %>%
  mutate(join = paste0("n per k=",n_per_k,", k=",k),
         sample_n = n_per_k*k, #total sample size
         deff = 1+icc*((n_per_k)-1), #design effect
         eff_n = sample_n / deff, #effective sample size
         iccf = paste0("ICC = ",factor(icc)),
         ctrlpf = factor(ctrl_prop),
         test = paste0(icc,ctrl_prop),
         kf = factor(k))
#saveRDS(power_plot,here("nonadapt_power.RDS"))
#plot the power
png(filename=here("Output","power_nonadapt.png"),width=8,height=6,res=300,units="in")
power_plot %>% mutate(trteff = case_when(trt_eff_scen == 1 ~ "Scenario 1",
                                         trt_eff_scen == 2 ~ "Scenario 2",
                                         trt_eff_scen == 3 ~ "Scenario 3")) %>%
  filter(k %in% c(5,10), ctrl_prop == 0.1) %>%
  ggplot(aes(y = bayesr,x = n_per_k)) +
  geom_point(aes(color = kf))+
  geom_line(aes(color = kf)) + 
  facet_wrap(~trteff+iccf) + 
  geom_hline(yintercept = 0.8) + 
  geom_hline(yintercept = 0.05, linetype="dashed") +
  theme_light() + 
  xlab("m") + 
  ylab("Power") + 
  labs(color = "k") +
  scale_color_manual(values=c("#48157F","#29AF7F"))
dev.off()
