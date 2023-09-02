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
### 4 chains, per chain: 750 total, 500 burn-ins

## MC Error of each iteration for parameters of interest
###simsum
MCerror <- outsim2 %>%
  group_by(property,variable) %>%
  mutate(trt = case_when(variable == "pred_prob_trt[1]" ~ t1,
                         variable == "pred_prob_trt[2]" ~ t2,
                         variable == "pred_prob_trt[3]" ~ t3,
                         variable == "pred_prob_trt[4]" ~ t4)) %>%
  simsum(estvarname="mean",true="trt",se="sd") 
kable(MCerror)

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

#plot the power over different features
#setting up the data
power_plot <- outsim3 %>%
  mutate(join = paste0("n per k=",n_per_k,", k=",k),
         sample_n = n_per_k*k,
         iccf = factor(icc),
         ctrlpf = factor(ctrl_prop),
         test = paste0(icc,ctrl_prop))

#the plot code
png(filename=here("Output","power_nonadapt.png"),width=8,height=4,res=300,units="in")
ggplot(data=power_plot,aes(y = bayesr,x = sample_n)) +
  geom_point(aes(color = ctrlpf))+
  geom_line(aes(linetype = iccf,color = ctrlpf)) + 
  facet_wrap(~trt_eff_scen) + 
  geom_hline(yintercept = 0.8) + 
  geom_hline(yintercept = 0.05, linetype="dashed") +
  theme_light() + 
  xlab("Sample size") + 
  ylab("Power") + 
  labs(color = "Control proportion",
       linetype = "ICC") +
  scale_color_manual(values=c("#48157F","#29AF7F"))
dev.off()
  
  
# zip plot illustrating coverage
zip_t1t4 <- outsim2 %>% 
  filter(contrast == "ttt1 - ttt4" & t1 != t4 & warn == "0") %>%
  group_by(property) %>%
  mutate(eff = abs(t1-t4),
         LCL = estimate - (1.96*SE),
         UCL = estimate + (1.96*SE)) %>%
  arrange(p.value) %>%
  mutate(centile = ((1900-row_number())/1900)*100,
         sig = ifelse(`p.value`< 0.05,"Yes","No"))

#https://osf.io/cbr72/
#for true eff = 0, you want a zip. For true eff != 0 then you will get a slide to each side.
p <- ggplot(data=zip_t1t4,aes(y=centile,color=sig)) + 
  geom_errorbar(aes(xmin=LCL,xmax=UCL,width=0)) + 
  theme(legend.position = "none") +
  scale_y_continuous(breaks = c(5,50,95)) +
  scale_color_manual(values=c("#993999", "#6178F0")) + 
  facet_wrap_paginate(~ property, ncol = 3, nrow = 2,scales = "free")

for(i in 1:n_pages(p)){
  p_save <- p + 
    facet_wrap_paginate(~ property, ncol = 3, nrow = 2, page = i,scales="free")
  ggsave(plot = p_save, filename = paste0('Output/zip_page_', i, '.jpg'),
         width=12,height=8,dpi=150,units="in")
}


#protocol for replacing lost simulations
#run the whole thing, count how many failed - run back up to 10k for each.