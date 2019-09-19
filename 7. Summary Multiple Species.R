

results <- data.frame(species=0:32,mean.distance=NA,max.distance=NA,Different=NA,Max.Different=NA)

for(i in 0:32){
  
  file <- paste0("../Results/Raw/SP#",i,".MatricesResults.RData")
  load(file)

  data <- connectivity.per.days.matrices[[1]]
                                         
  results[i+1,2] <- mean(data$Distance)
  results[i+1,3] <- max(data$Distance)
  results[i+1,4] <- mean(data$Differantiation)
  results[i+1,5] <- max(data$Differantiation)
  
}

write.csv(results,file=paste0(project.folder,"Results/Summary Distances Bettween Pops.csv"))
