#AUTHOR: Erin Nolan
#TITLE: PROJECT 2 FUNCTION MAKE CLUSTERS
#PURPOSE: FUNCTION TO MAKE THE GENERATED DATASET

pacman::p_load(tidyverse,purrr,lme4,MASS,lsmeans)

# Make the generated dataset---------------------------
# makeClusters function
# t = number of arms
# nid = number of participants per arm
# specify the number of clusters for each t
makeClusters <- function(t,nid,...){
  int <- list()
  clusters <- list(...) #the number of clusters in each group
  nclusters <- cumsum(clusters) #cumulative n for the for loop
  for(i in 1:t){
    if (i==1){ #first group has cluster starting at 1
      int[[i]] <- expand.grid(iid = c(1:nid), site = seq(1,nclusters[i],1), trt = i)
    }
    else {
      int[[i]] <- expand.grid(iid = c(1:nid), site = seq(nclusters[i-1]+1,nclusters[i],1), trt = i)
    }
  }
  groups <- map_df(int, ~as.data.frame(.)) #unlist them and make them a dataframe
  return(groups)
}
