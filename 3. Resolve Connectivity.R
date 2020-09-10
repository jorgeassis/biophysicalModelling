## ------------------------------------------------------------------------------------------------------------------
## PlankTonic
## Assis et al., 2018
## ------------------------------------------------------------------------------------------------------------------
##
## ------------------------------------------------------------------------------------------------------------------

rm(list=(ls()[ls()!="v"]))
gc(reset=TRUE)
source("../Project Config 0.R")

## --------------------------------------------------------------------------------------------------------------
##
##
## 
## --------------------------------------------------------------------------------------------------------------

## Test if connectivity exists

sql <- dbConnect(RSQLite::SQLite(), paste0(sql.directory,"/",project.name,"SimulationResults.sql"))
dbGetQuery(sql, "SELECT * FROM 'Connectivity' LIMIT 5 ")
dbDisconnect(sql)

## ------------------------------------
## Resolve connectivity

sql <- dbConnect(RSQLite::SQLite(), paste0(sql.directory,"/",project.name,"SimulationResults.sql"))
cell.to.process <- 1:nrow(dbReadTable(sql, "SourceSinkSites"))
n.particles.per.cell <- dbReadTable(sql, "Parameters")$n.particles.per.cell[1]
n.new.particles.per.day <- dbReadTable(sql, "Parameters")$n.new.particles.per.day[1]
n.steps.per.day <- dbReadTable(sql, "Parameters")$n.hours.per.day[1]
dbDisconnect(sql)

## ------------------

particles.reference.bm.desc <- dget( paste0(project.folder,"/Results/InternalProc/particles.reference.desc"))

## ------------------

cl.2 <- makeCluster(10 , type="FORK")
registerDoParallel(cl.2)

all.connectivity.pairs.to.sql <- foreach(cell.id.ref.f=cell.to.process, .verbose=FALSE, .combine = rbind ,  .packages=c("gstat","raster","data.table","FNN","bigmemory")) %dopar% { # 
  
  particles.reference.bm.i <- attach.big.matrix(particles.reference.bm.desc)
  
  connectivity.temp.m <- particles.reference.bm.i[ mwhich(particles.reference.bm.i,2,list(cell.id.ref.f), list('eq')) , ]
  connectivity.temp.m <- data.frame(connectivity.temp.m)
  colnames(connectivity.temp.m) <- c("id","start.cell","start.year","start.month","start.day","pos.lon","pos.lat","pos.alt","state","t.start","t.finish","cell.rafted")
  
  connectivity.temp.m <- connectivity.temp.m[connectivity.temp.m$cell.rafted != 0 & connectivity.temp.m$state == 2,]
  connectivity.temp.m <- data.frame(connectivity.temp.m,travel.time= (1 + connectivity.temp.m$t.finish - connectivity.temp.m$t.start) / n.steps.per.day)
  connectivity.temp.m <- as.data.table(connectivity.temp.m)
  connectivity.pairs.to.sql <- data.frame()
  
  for(y in unique(connectivity.temp.m$start.year) ) {
    
    connectivity.temp <- connectivity.temp.m[ start.year == y , ]
    
    for( cell.id.ref.t in unique(connectivity.temp[ , cell.rafted ]) ) {
      
      connectivity.pairs.to.sql <- rbind(connectivity.pairs.to.sql,
                                         
                                         data.frame(  Pair.from = cell.id.ref.f,
                                                      Pair.to = cell.id.ref.t,
                                                      Number.events = nrow(connectivity.temp[ cell.rafted == cell.id.ref.t,]),
                                                      Time.mean = mean(connectivity.temp[ cell.rafted == cell.id.ref.t,]$travel.time),
                                                      Time.min = min(connectivity.temp[ cell.rafted == cell.id.ref.t,]$travel.time),
                                                      Time.max = max(connectivity.temp[ cell.rafted == cell.id.ref.t,]$travel.time),
                                                      Time.sd = sd(connectivity.temp[ cell.rafted == cell.id.ref.t,]$travel.time),
                                                      Probability = nrow(connectivity.temp[ cell.rafted == cell.id.ref.t,]) / round( n.particles.per.cell / length(unique(connectivity.temp.m$start.year)) ),
                                                      Year = y ) )
    }
    
  }
  
  connectivity.pairs.to.sql[is.na(connectivity.pairs.to.sql)] <- 0
  return( connectivity.pairs.to.sql )
  
}

stopCluster(cl.2) ; rm(cl.2) ; gc(reset=TRUE)

# -----------------------------------------

# Save pairs to SQL 

sql <- dbConnect(RSQLite::SQLite(), paste0(sql.directory,"/",project.name,"SimulationResults.sql"))
dbWriteTable(sql, "Connectivity", all.connectivity.pairs.to.sql , overwrite=TRUE, append=FALSE)
dbDisconnect(sql)

rm(all.connectivity.pairs.to.sql) ; gc()

## --------------------------------------------------------------------------------------------------------------
## --------------------------------------------------------------------------------------------------------------

# Direct Overall Connectivity matrix (mean of all years)

resultsFolder <- "Results" # Results

sql <- dbConnect(RSQLite::SQLite(), paste0(sql.directory,"/",project.name,"SimulationResults.sql"))
Connectivity <- data.table(dbReadTable(sql, "Connectivity"))
source.sink.xy <- dbReadTable(sql, "SourceSinkSites")
plot(source.sink.xy[,2:3] )
dbDisconnect(sql)
Connectivity

# Subset Years
# Connectivity <- Connectivity[ Year == 2017, ]

Connectivity <- Connectivity[ , j=list(mean(Probability, na.rm = TRUE) , sd(Probability, na.rm = TRUE) , max(Probability, na.rm = TRUE) , mean(Time.mean, na.rm = TRUE) , sd(Time.mean, na.rm = TRUE) , max(Time.mean, na.rm = TRUE) , mean(Number.events, na.rm = TRUE) , sd(Number.events, na.rm = TRUE) , max(Number.events, na.rm = TRUE) ) , by = list(Pair.from,Pair.to)]
colnames(Connectivity) <- c("Pair.from" , "Pair.to" , "Mean.Probability" , "SD.Probability" , "Max.Probability" , "Mean.Time" , "SD.Time" , "Max.Time" , "Mean.events" , "SD.events" , "Max.events" )
Connectivity[is.na(Connectivity)] <- 0
Connectivity ; gc()

## -------

source.sink.id <- source.sink.xy$cells.id[which(source.sink.xy$source == 1)] 
plot(source.sink.xy[which(source.sink.xy$source == 1),2:3] )

source.sink.xy <- source.sink.xy[source.sink.xy$cells.id %in% source.sink.id,]
source.sink.bm <- as.big.matrix(as.matrix(source.sink.xy))
write.big.matrix(source.sink.bm, paste0(project.folder,"/",resultsFolder,"/source.sink.bm"))

Connectivity <- Connectivity[Connectivity$Pair.from %in% source.sink.id & Connectivity$Pair.to %in% source.sink.id,]
Connectivity.bm <- as.big.matrix(as.matrix(Connectivity))
write.big.matrix(Connectivity.bm, paste0(project.folder,"/",resultsFolder,"/Connectivity.bm")) 

## ------------------------------------------------------------------------------------------------------
## ------------------------------------------------------------------------------------------------------
## Square matrices : Maybe limited in terms of memory

matrix.size <- length(unique(c(Connectivity$Pair.from,Connectivity$Pair.to)))
Connectivity.matrix.probability <- matrix(NA,ncol=matrix.size,nrow=matrix.size)

for(i in 1:matrix.size){
  
  cell.from <- source.sink.xy$cells.id[i]
  cells.to <- Connectivity[Connectivity$Pair.from %in% cell.from,Pair.to]
  if( length(cell.from) == 0 | length(cells.to) == 0 ) { next }
  
  Connectivity.matrix.probability[cell.from,cells.to] <- Connectivity[Connectivity$Pair.from %in% cell.from,Mean.Probability]
  
}

Connectivity.matrix.probability.bm <- as.big.matrix(as.matrix(Connectivity.matrix.probability))
write.big.matrix(Connectivity.matrix.probability.bm, paste0(project.folder,"/Results/Connectivity.matrix.probability.bm"))

rm(Connectivity.matrix.probability) ; rm(Connectivity.matrix.probability.bm) ; gc()

## --------------------------------

Connectivity.matrix.time <- matrix(NA,ncol=matrix.size,nrow=matrix.size)

for(i in 1:matrix.size){
  
  cell.from <- source.sink.xy$cells.id[i]
  cells.to <- Connectivity[Connectivity$Pair.from %in% cell.from,Pair.to]
  
  if( length(cell.from) == 0 | length(cells.to) == 0 ) { next }
  
  Connectivity.matrix.time[cell.from,cells.to] <- Connectivity[Connectivity$Pair.from %in% cell.from,Mean.Time]
  
}

Connectivity.matrix.time.bm <- as.big.matrix(as.matrix(Connectivity.matrix.time))
write.big.matrix(Connectivity.matrix.time.bm, paste0(project.folder,"/Results/Connectivity.matrix.time.bm"))

rm(Connectivity.matrix.time) ; rm(Connectivity.matrix.time.bm) ; gc()

## --------------------------------

Connectivity.matrix.max.time <- matrix(NA,ncol=matrix.size,nrow=matrix.size)

for(i in 1:matrix.size){
  
  cell.from <- source.sink.xy$cells.id[i]
  cells.to <- Connectivity[Connectivity$Pair.from %in% cell.from,Pair.to]
  
  if( length(cell.from) == 0 | length(cells.to) == 0 ) { next }
  
  Connectivity.matrix.max.time[cell.from,cells.to] <- Connectivity[Connectivity$Pair.from %in% cell.from,Max.Time]
  
}

Connectivity.matrix.max.time.bm <- as.big.matrix(as.matrix(Connectivity.matrix.max.time))
write.big.matrix(Connectivity.matrix.max.time.bm, paste0(project.folder,"/Results/Connectivity.matrix.max.time.bm"))

rm(Connectivity.matrix.max.time) ; rm(Connectivity.matrix.max.time.bm) ; gc()

## --------------------------------------------------------------------------------------------------------------
## --------------------------------------------------------------------------------------------------------------
## End of Code [!]
