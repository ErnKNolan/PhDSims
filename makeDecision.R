

makeDecision <- function(properties = properties, interim_res = interim_res, j=j,adaption=adaption){

    #Determining what (if any) treatment group to drop
    int_drop <- interim_res %>% 
      filter(variable %in% c("pp_trt1","pp_trt2","pp_trt3","fu_trt1","fu_trt2","fu_trt3")) %>% 
      dplyr::select(variable,mean) %>%
      pivot_wider(names_from = variable, values_from = mean)
    
    #finding the smallest pred prob
    int_drop$drop <- apply(int_drop[,c(1:3)], 1, min, na.rm = TRUE)
    #determine if stopping rule met
    int_drop$stop <- ifelse(int_drop$fu_trt1 < 0.05 & int_drop$fu_trt2 < 0.05 & int_drop$fu_trt3 < 0.05, 1, 0)
    
    #dropping only if pred prob less than 0.05
    int_drop <- int_drop %>%
      mutate(drop = ifelse(drop > 0.05, NA, drop),
             droptrt = pmap_chr(list(pp_trt1, pp_trt2, pp_trt3, drop), 
                                ~ifelse(any(c(..1, ..2, ..3) %in% ..4),
                                        c("trt1", "trt2", "trt3")[c(..1, ..2, ..3) %in% ..4],
                                        "none")))
    properties <- properties[j,]
    properties$drop <- int_drop$droptrt
    properties$stop <- int_drop$stop
    stop <- int_drop$stop
#This makes decisions for trial properties past the interim based on the adaptions
#Can either be, early stopping, arm dropping, or both
      if (adaption == "both") {
        if (stop == 1) {
          properties_int <- properties %>% 
          mutate(kt1 = interim,
                 kt2 = interim,
                 kt3 = interim)
          
        } else if (stop == 0) {
          properties_int <- properties %>% 
          mutate(kt1 = case_when(
              drop == "trt1" ~ interim,
              drop %in% c("trt2", "trt3") ~ k + ceiling((k - interim) / 2),
              drop == "none" ~ k),
            kt2 = case_when(
              drop == "trt2" ~ interim,
              drop %in% c("trt1", "trt3") ~ k + floor((k - interim) / 2),
              drop == "none" ~ k),
            kt3 = case_when(
              drop == "trt3" ~ interim,
              drop %in% c("trt1", "trt2") ~ k + floor((k - interim) / 2),
              drop == "none" ~ k))
        }
      } else if (adaption == "arm_dropping") {
        properties_int <- properties %>% 
        mutate(kt1 = case_when(
            drop == "trt1" ~ interim,
            drop %in% c("trt2", "trt3") ~ k + ceiling((k - interim) / 2),
            drop == "none" ~ k),
          kt2 = case_when(
            drop == "trt2" ~ interim,
            drop %in% c("trt1", "trt3") ~ k + floor((k - interim) / 2),
            drop == "none" ~ k),
          kt3 = case_when(
            drop == "trt3" ~ interim,
            drop %in% c("trt1", "trt2") ~ k + floor((k - interim) / 2),
            drop == "none" ~ k))
        
      } else if (adaption == "early_stopping") {
        if (stop == 1) {
          properties_int <- properties %>% 
          mutate(kt1 = interim,
                 kt2 = interim,
                 kt3 = interim)
          
        } else {
          properties_int <- properties %>% 
          mutate(kt1 = k,
                 kt2 = k,
                 kt3 = k)
        }
      }
    
return(properties_int)

}