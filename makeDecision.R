

makeDecision <- function(properties = properties, interim_res = interim_res, j=j,adaption=adaption,drop_cut=drop_cut,stop_cut=stop_cut,ties=ties){

  #random number for allocation if it happens
  rand <- as.vector(sample(c(0,1),2,replace=FALSE))
      #Determining what (if any) treatment group to drop
    int_drop <- interim_res %>% 
      filter(variable %in% c("pp_trt2","pp_trt3","pp_trt4","pred_prob_trt[2]","pred_prob_trt[3]","pred_prob_trt[4]")) %>% 
      dplyr::select(variable,mean) %>%
      pivot_wider(names_from = variable, values_from = mean) %>%
      #I want the columns in a specific order to keep the old code the same
      dplyr::select(pp_trt2,pp_trt3,pp_trt4,`pred_prob_trt[2]`,`pred_prob_trt[3]`,`pred_prob_trt[4]`)
    
    #finding the smallest pred prob
    int_drop$drop <- apply(int_drop[,c(1:3)], 1, min, na.rm = TRUE)
    #determine if stopping rule met
    int_drop$stop <- ifelse(int_drop$pp_trt2 < stop_cut & int_drop$pp_trt3 < stop_cut & int_drop$pp_trt4 < stop_cut, 1, 0)
    
    #dropping only if pp less than drop_cut
    #Optimisation option - if tie then drop higher constraint
    if(ties == "ties_opt"){
    int_drop <- int_drop %>%
      mutate(drop = ifelse(drop > drop_cut, NA, drop),
             droptrt = pmap_chr(list(pp_trt2, pp_trt3, pp_trt4, drop), 
                                ~ifelse(any(c(..1, ..2, ..3) %in% ..4),
                                        c("trt2", "trt3", "trt4")[c(..1, ..2, ..3) %in% ..4],
                                        "none")))
    } else{
    #other option when there's a tie - use lower pred prob
    int_drop <- int_drop %>%
      mutate(drop = ifelse(drop > drop_cut, NA, drop),
             total = rowSums(across(pp_trt2:pp_trt4) == drop),
             droptrt = ifelse(total <= 1, 
                              pmap_chr(list(pp_trt2, pp_trt3, pp_trt4, drop), 
                                       ~ifelse(any(c(..1, ..2, ..3) %in% ..4),
                                               c("trt2", "trt3", "trt4")[c(..1, ..2, ..3) %in% ..4],
                                               "none")),
                              case_when(min(`pred_prob_trt[2]`,`pred_prob_trt[3]`,`pred_prob_trt[4]`) == `pred_prob_trt[2]` ~ "trt2",
                                        min(`pred_prob_trt[2]`,`pred_prob_trt[3]`,`pred_prob_trt[4]`) == `pred_prob_trt[3]` ~ "trt3",
                                        min(`pred_prob_trt[2]`,`pred_prob_trt[3]`,`pred_prob_trt[4]`) == `pred_prob_trt[4]` ~ "trt4")),
             #if none to drop then total and droptrt are NA - so assign it 'none'
             droptrt = ifelse(is.na(droptrt),"none",droptrt))
    }
    
    properties <- properties[j,]
    properties$drop <- int_drop$droptrt
    properties$stop <- int_drop$stop
    stop <- int_drop$stop
    
#This makes decisions for trial properties past the interim based on the adaptions
#Can either be, early stopping, arm dropping, or both
      if (adaption == "both") {
        if (stop == 1) {
          properties_int <- properties %>% 
          mutate(kt2 = interim,
                 kt3 = interim,
                 kt4 = interim)
          
        } else if (stop == 0) {
          if(properties$drop == "trt2"){
            rand <- append(rand,0,after=0)
          } else if(properties$drop == "trt3"){
              rand <- append(rand,0,after=1)
            } else if(properties$drop == "trt4"){
              rand <- append(rand,0,after=2)
            } else{rand <- rand}
          if((properties$k - properties$interim) %% 2 == 0){
            rand <- replace(rand,rand==1,0)
          } else(rand <- rand)
          properties_int <- properties %>% 
          mutate(kt2 = case_when(
              drop == "trt2" ~ interim,
              drop %in% c("trt3", "trt4") ~ k + floor((k-interim)/2) + rand[1],
              drop == "none" ~ k),
            kt3 = case_when(
              drop == "trt3" ~ interim,
              drop %in% c("trt2", "trt4") ~ k + floor((k-interim)/2) + rand[2],
              drop == "none" ~ k),
            kt4 = case_when(
              drop == "trt4" ~ interim,
              drop %in% c("trt2", "trt3") ~ k + floor((k-interim)/2) + rand[3],
              drop == "none" ~ k))
        }
      } else if (adaption == "arm_dropping") {
        if(properties$drop == "trt2"){
          rand <- append(rand,0,after=0)
        } else if(properties$drop == "trt3"){
          rand <- append(rand,0,after=1)
        } else if(properties$drop == "trt4"){
          rand <- append(rand,0,after=2)
        } else{rand <- rand}
        if((properties$k - properties$interim) %% 2 == 0){
          rand <- replace(rand,rand==1,0)
        } else(rand <- rand)
        properties_int <- properties %>% 
        mutate(kt2 = case_when(
            drop == "trt2" ~ interim,
            drop %in% c("trt3", "trt4") ~ k + floor((k-interim)/2) + rand[1],
            drop == "none" ~ k),
          kt3 = case_when(
            drop == "trt3" ~ interim,
            drop %in% c("trt2", "trt4") ~ k + floor((k-interim)/2) + rand[2],
            drop == "none" ~ k),
          kt4 = case_when(
            drop == "trt4" ~ interim,
            drop %in% c("trt2", "trt3") ~ k + floor((k-interim)/2) + rand[3],
            drop == "none" ~ k))
        
      } else if (adaption == "early_stopping") {
        if (stop == 1) {
          properties_int <- properties %>% 
          mutate(kt2 = interim,
                 kt3 = interim,
                 kt4 = interim)
          
        } else {
          properties_int <- properties %>% 
          mutate(kt2 = k,
                 kt3 = k,
                 kt4 = k)
        }
      }
    
return(properties_int)

}