
runSimTrial <- function(properties = properties, mod = mod, outdir=outdir, j=j, adaption=adaption,drop_cut=drop_cut,stop_cut=stop_cut) {
  
  #making clusters for the interim dataset
  intclusters <- makeClusters(t=4,nid=properties$n_per_k[j],
                              t1=properties$interim[j],t2=properties$interim[j],
                              t3=properties$interim[j],t4=properties$interim[j])
  #analysing the interim data
  interim <- testInterim(expdat=intclusters,t=4,mod=mod,outdir=outdir,
                    rho=properties$icc[j],t1=properties$t1[j],t2=properties$t2[j],t3=properties$t3[j],t4=properties$t4[j])
  
  #save the interim results and sim
  resp <- unlist(interim[["resp"]]) 
  interim_res <- interim[[1]]
  
  #add the simulated results into the interim dataset
  resp_clus <- cbind(intclusters,resp)
  
  #making the decision and getting the clusters
  properties_int <- makeDecision(properties = properties, interim_res = interim_res, j=j, adaption = adaption,drop_cut = drop_cut,stop_cut=stop_cut)

  fullclusters <- makeClusters(t=4,nid=properties_int$n_per_k,
                               t1=properties_int$kt1,t2=properties_int$kt2,
                               t3=properties_int$kt3,t4=properties_int$k)
  
  
  #need to append the interim data and new data together
  #then run the full trial model
  full <- testFull(expdat=fullclusters,t=4,mod=mod,outdir=outdir,int_dat=resp_clus,
                 rho=properties_int$icc,t1=properties_int$t1,t2=properties_int$t2,t3=properties_int$t3,t4=properties_int$t4)

  res <- list(interim=interim[[1]],full=full[[1]],properties=properties_int)
  return(res)
}