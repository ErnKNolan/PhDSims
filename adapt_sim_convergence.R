#AUTHOR: Erin Nolan
#TITLE: PROJECT 2 SIMULATION CONVERGENCE
#PURPOSE: PLOTTING THE CONVERGENCE OF THE SIMULATIONS
# install.packages("devtools")
#devtools::install_github("matherealize/looplot")

pacman::p_load(here,future.apply,tictoc,car,ggforce,dplyr,cmdstanr,looplot,forcats,tidyr,gtsummary)

#outsim2 <- readRDS(here("Data","outsim_adapt_opt.RDS"))
#outsim2 <- readRDS(here("Data","outsim_adapt_prob.RDS"))
#outsim2 <- readRDS(here("Data","outsim_nonadapt.RDS"))
#outsim2 <- readRDS(here("Data","interim_sim.RDS))
#properties2 <- readRDS(here("Data","properties2.RDS"))

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

#Plot of ess
outsim2 %>% 
  filter(variable %in% c("beta_trt[2]","beta_trt[3]","beta_trt[4]")) %>% 
  ggplot(aes(x=ess_bulk,group=n_per_k)) + geom_density(aes(fill=n_per_k),alpha=0.3,position="identity") + 
  facet_wrap(~variable)
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

#PERFORMANCE MEASURES
icc_perf <- convergence %>% group_by(icc) %>% 
  summarise(m_ess_bulk = mean(m_ess_bulk),
            m_ess_bulk300 = mean(m_ess_bulk300),
            m_ess_tail = mean(m_ess_tail),
            mrhat = mean(mrhat),
            mMCSE = mean(MCSE)) %>%
  rename(property = icc) %>%
  mutate(property = case_when(property == 0.05 ~ "ICC = 0.05",
                              property == 0.2 ~ "ICC = 0.2"))

k_perf <- convergence %>% group_by(k) %>% 
  summarise(m_ess_bulk = mean(m_ess_bulk),
            m_ess_bulk300 = mean(m_ess_bulk300),
            m_ess_tail = mean(m_ess_tail),
            mrhat = mean(mrhat),
            mMCSE = mean(MCSE)) %>%
  rename(property = k) %>%
  mutate(property = case_when(property == 3 ~ "k = 3",
                              property == 5 ~ "k = 5",
                              property == 10 ~ "k = 10"))

n_perf <- convergence %>% group_by(n_per_k) %>% 
  summarise(m_ess_bulk = mean(m_ess_bulk),
            m_ess_bulk300 = mean(m_ess_bulk300),
            m_ess_tail = mean(m_ess_tail),
            mrhat = mean(mrhat),
            mMCSE = mean(MCSE)) %>%
  rename(property = n_per_k) %>%
  mutate(property = case_when(property == 5 ~ "n = 5",
                              property == 25 ~ "n = 25",
                              property == 50 ~ "n = 50",
                              property == 75 ~ "n = 75",
                              property == 100 ~ "n = 100"))

trt_perf <- convergence %>% group_by(trt_eff_scen) %>% 
  summarise(m_ess_bulk = mean(m_ess_bulk),
            m_ess_bulk300 = mean(m_ess_bulk300),
            m_ess_tail = mean(m_ess_tail),
            mrhat = mean(mrhat),
            mMCSE = mean(MCSE)) %>%
  rename(property = trt_eff_scen) %>%
  mutate(property = case_when(property == 1 ~ "Large effect",
                              property == 2 ~ "Mid/small effect",
                              property == 3 ~ "Null scenario"))

ov_perf <- convergence %>% 
  #group_by(trt_eff_scen) %>%
  summarise(m_ess_bulk = mean(m_ess_bulk),
            m_ess_bulk300 = mean(m_ess_bulk300),
            m_ess_tail = mean(m_ess_tail),
            mrhat = mean(mrhat),
            mMCSE = mean(MCSE)) %>%
  mutate(property = "Overall") %>%
  dplyr::select(property,m_ess_bulk,m_ess_bulk300,m_ess_tail,mrhat,mMCSE)

#COMBINE THE PERFORMANCE MEASURES
performance <- rbind(ov_perf,icc_perf,k_perf,n_perf,trt_perf) %>%
  mutate(m_ess_bulk = m_ess_bulk*100,
         m_ess_bulk300 = m_ess_bulk300*100,
         m_ess_tail = m_ess_tail*100,
         mrhat = mrhat*100)


tbl_summary(performance,by="property",
            statistic = list(all_continuous() ~ "{mean}")) 

converge_head <- convergence %>% 
  mutate(m_ess_bulk = m_ess_bulk*100,
         m_ess_bulk300 = m_ess_bulk300*100,
         m_ess_tail = m_ess_tail*100,
         mrhat = mrhat*100) %>%
  arrange(desc(mrhat)) %>% 
  head(3)
  #%>% #The remaining code will save the output but needs Chrome. 
  #as_gt() %>%
  #gt::gtsave(filename = here("Output","performance_tab_prob15.png"))
