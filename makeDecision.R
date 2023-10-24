

makeDecision <- function(properties = properties, interim_res = interim_res, j=j){}


#Determining what (if any) treatment group to drop
int_drop <- interim_res %>% 
  filter(variable %in% c("pp_trt1","pp_trt2","pp_trt3")) %>% 
  dplyr::select(variable,mean,ess_bulk,ess_tail) %>%
  pivot_wider(names_from = variable, values_from = mean)

#finding the smallest pred prob
int_drop$drop <- apply(int_drop[,c(3:5)], 1, min, na.rm = TRUE)

#dropping only if pred prob less than 0.05
int_drop <- int_drop %>%
  mutate(drop = ifelse(drop > 0.05, NA, drop),
         droptrt = pmap_chr(list(pp_trt1, pp_trt2, pp_trt3, drop), 
                            ~ifelse(any(c(..1, ..2, ..3) %in% ..4),
                                    c("trt1", "trt2", "trt3")[c(..1, ..2, ..3) %in% ..4],
                                    "none")))
properties <- properties[j,]
properties$drop <- int_drop$droptrt
#the different outcomes to change properties
properties_int <- properties %>% 
  mutate(kt1 = case_when(drop == "trt1" ~ interim,
                       drop %in% c("trt2","trt3") ~ k+1,
                       drop == "none" ~ k),
       kt2 = case_when(drop == "trt2" ~ interim,
                       drop %in% c("trt1","trt3") ~ k+1,
                       drop == "none" ~ k),
       kt3 = case_when(drop == "trt3" ~ interim,
                       drop %in% c("trt1","trt2") ~ k+1,
                       drop == "none" ~ k))

return(properties_int)