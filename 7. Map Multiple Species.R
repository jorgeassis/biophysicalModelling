
source.sink.xy.i <- source.sink.xy[source.sink.xy$source == 1 ,]

files <- list.files("../Data/Differentiation/",pattern="Coords",full.names = TRUE)

results <- data.frame()
results.i <- data.frame()

for(f in files) {
  
  r.i <- read.table(f,sep=";",header=T)
  colnames(r.i) <- c("Site","Lon","Lat")
  
  dists <- spDists(as.matrix(source.sink.xy.i[,2:3]),as.matrix(r.i[,2:3]),longlat = TRUE)
  
  to.remove <- which(apply(dists,2,min) > 100)
  if(length(to.remove) > 0) { r.i <- r.i[-to.remove,]}
  
  results <- rbind(results,r.i)
  results.i <- rbind(results.i,data.frame(name=f,pops=nrow(read.table(f,sep=";"))))
  
}

plot(results[,2:3])

nrow(results)
results
results.i

write.csv(results,file=paste0(project.folder,"Results/Summary All Sampling Sites.csv"))
write.csv(results.i,file=paste0(project.folder,"Results/Summary All Sampling Pops.csv"))

for(f in files) {
  
  r.i <- read.table(f,sep=";",header=T)
  colnames(r.i) <- c("Site","Lon","Lat")
  print( plot(r.i[,2:3],main=f) )
  readline(prompt="Press [enter] to continue")
  
  
}
