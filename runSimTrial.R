

intclusters <- makeClusters(t=4,nid=properties$n_per_k[j],
                            t1=properties$interim[j],t2=properties$interim[j],
                            t3=properties$interim[j],t4=properties$interim[j])

interim <- testss(expdat=intclusters,t=4,mod=mod,outdir=outdir,
                  rho=properties$icc[j],t1=properties$t1[j],t2=properties$t2[j],t3=properties$t3[j],t4=properties$t4[j])

#save the interim results and sim
resp <- unlist(interim[["resp"]]) 
interim_res <- interim[[1]]

#add the simulated results into the interim dataset
resp_clus <- cbind(intclusters,resp)

#some cleaning to get to decision
properties_int <- makeDecision()

fullclusters <- makeClusters(t=4,nid=properties_int$n_per_k,
                             t1=properties_int$kt1,t2=properties_int$kt2,
                             t3=properties_int$kt3,t4=properties_int$k)


#need to append the interim data and new data together

full <- testss(expdat=fullclusters[[j]],t=4,mod=mod,outdir=outdir,
               rho=properties$icc[j],t1=properties$t1[j],t2=properties$t2[j],t3=properties$t3[j],t4=properties$t4[j])

return(full)