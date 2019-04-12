## ------------------------------------------------------------------------------------------------------------------
## PlankTonic
## Assis et al., 2018
## ------------------------------------------------------------------------------------------------------------------
##
## ------------------------------------------------------------------------------------------------------------------

rm(list=(ls()[ls()!="v"]))
gc(reset=TRUE)
source("0. Project Config.R")

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

particles.reference.bm.desc <- dget( paste0(project.folder,"/InternalProc/particles.reference.external.desc"))

## ------------------

cl.2 <- makeCluster(number.cores)
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

## --------------------------------------------------------------------------------------------------------------
## --------------------------------------------------------------------------------------------------------------

# Direct Overall Connectivity matrix (mean of all years)

sql <- dbConnect(RSQLite::SQLite(), paste0(sql.directory,"/",project.name,"SimulationResults.sql"))
Connectivity <- data.table(dbReadTable(sql, "Connectivity"))
source.sink.xy <- dbReadTable(sql, "SourceSinkSites")
dbDisconnect(sql)
Connectivity

Connectivity <- Connectivity[ , j=list(mean(Probability, na.rm = TRUE) , sd(Probability, na.rm = TRUE) , max(Probability, na.rm = TRUE) , mean(Time.mean, na.rm = TRUE) , sd(Time.mean, na.rm = TRUE) , max(Time.mean, na.rm = TRUE) , mean(Number.events, na.rm = TRUE) , sd(Number.events, na.rm = TRUE) , max(Number.events, na.rm = TRUE) ) , by = list(Pair.from,Pair.to)]
colnames(Connectivity) <- c("Pair.from" , "Pair.to" , "Mean.Probability" , "SD.Probability" , "Max.Probability" , "Mean.Time" , "SD.Time" , "Max.Time" , "Mean.events" , "SD.events" , "Max.events" )
Connectivity[is.na(Connectivity)] <- 0
Connectivity

Connectivity.bm <- as.big.matrix(as.matrix(Connectivity))
write.big.matrix(Connectivity.bm, "/Volumes/Laminaria/Dropbox/Manuscripts/Transport Simulation in Eastern Asia/Results/Connectivity.bm")

## --------------------------------

source.sink.id <- 1:nrow(source.sink.xy)

source.sink.bm <- as.big.matrix(as.matrix(source.sink.xy))
write.big.matrix(source.sink.bm, "/Volumes/Laminaria/Dropbox/Manuscripts/Transport Simulation in Eastern Asia/Results/source.sink.bm")

## --------------------------------

Connectivity.matrix.probability <- acast(Connectivity[,.(Pair.from , Pair.to ,  Mean.Probability)], Pair.from~Pair.to, value.var="Mean.Probability")
Connectivity.matrix.time <- acast(Connectivity[,.(Pair.from , Pair.to ,  Mean.Time)], Pair.from~Pair.to, value.var="Mean.Time")
Connectivity.matrix.max.time <- acast(Connectivity[,.(Pair.from , Pair.to ,  Max.Time)], Pair.from~Pair.to, value.var="Max.Time")
dim(Connectivity.matrix.probability)

View(Connectivity.matrix.probability)
View(Connectivity.matrix.time)

Connectivity.matrix.probability.bm <- as.big.matrix(as.matrix(Connectivity.matrix.probability))
write.big.matrix(Connectivity.matrix.probability.bm, "/Volumes/Laminaria/Dropbox/Manuscripts/Transport Simulation in Eastern Asia/Results/Connectivity.matrix.probability.bm")

Connectivity.matrix.time.bm <- as.big.matrix(as.matrix(Connectivity.matrix.time))
write.big.matrix(Connectivity.matrix.time.bm, "/Volumes/Laminaria/Dropbox/Manuscripts/Transport Simulation in Eastern Asia/Results/Connectivity.matrix.time.bm")

Connectivity.matrix.max.time.bm <- as.big.matrix(as.matrix(Connectivity.matrix.max.time))
write.big.matrix(Connectivity.matrix.max.time.bm, "/Volumes/Laminaria/Dropbox/Manuscripts/Transport Simulation in Eastern Asia/Results/Connectivity.matrix.max.time.bm")

## --------------------------------------------------------------------------------------------------------------
## --------------------------------------------------------------------------------------------------------------

# Stepping stone Connectivity matrix

Connectivity <- read.big.matrix("/Volumes/Laminaria/Dropbox/Manuscripts/Transport Simulation in Eastern Asia/Results/Connectivity.bm")
Connectivity <- data.table(Connectivity[,])
colnames(Connectivity) <- c("Pair.from" , "Pair.to" , "Probability" , "SD.Probability" , "Max.Probability" , "Mean.Time" , "SD.Time" , "Time.max" , "Mean.events" , "SD.events" , "Max.events" )

source.sink.xy <- read.big.matrix("/Volumes/Laminaria/Dropbox/Manuscripts/Transport Simulation in Eastern Asia/Results/source.sink.bm")
source.sink.xy <- data.table(source.sink.xy[,])
colnames(source.sink.xy) <- c("Pair" , "Lon" , "Lat" , "Source" )
source.sink.xy <- source.sink.xy[Source == 1,]

Connectivity.matrix.ss <- matrix(0,ncol=nrow(source.sink.xy),nrow=nrow(source.sink.xy))
dim(Connectivity.matrix.ss)

network <- produce.network("Prob",Connectivity,60,FALSE,0,source.sink.xy,0)
network.x <- network[[2]]
connectivity.x <- network[[1]]

## ------------------------

for( from in source.sink.xy[,Pair] ) {
  
  ptm <- proc.time()
  
  cl.3 <- makeCluster(10) ; registerDoParallel(cl.3)
  
  connectivity.f <- foreach(to=source.sink.xy[,Pair], .verbose=FALSE, .packages=c("data.table","sp","gdistance","igraph")) %dopar% { 

        possible.paths.y <- get.shortest.paths(network.x,as.character( from ) , as.character( to ),mode="out")$vpath
        stones.t <- as.numeric(names(possible.paths.y[[1]]))
        stones.t.interm <- cbind(stones.t[-length(stones.t)],stones.t[-1])
        path.values <- apply( stones.t.interm , 1 , function(z) { connectivity.x[ connectivity.x[,1] == z[1] & connectivity.x[,2] == z[2] , 3 ][1] }   )
        
        if( length(path.values) > 0 ) { return(apply( t(path.values) , 1 , prod )) }
        if( length(path.values) == 0) { return(0) }
        
  }
  
  stopCluster(cl.3) ; rm(cl.3) ; gc()
  

  proc.time() - ptm
  
  Connectivity.matrix.ss[from,] <- do.call(rbind,connectivity.f)
    
}

## --------------------------------------------------------------------------------------------------------------
## --------------------------------------------------------------------------------------------------------------