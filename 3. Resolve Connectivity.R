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

## ------------------------------------
## Resolve connectivity

sql <- dbConnect(RSQLite::SQLite(), paste0(sql.directory,"/",project.name,"SimulationResults.sql"))
cell.to.process <- 1:nrow(dbReadTable(sql, "SourceSinkSites"))
particles.reference <- as.data.table(dbReadTable(sql, "ReferenceTable"))
n.particles.per.cell <- dbReadTable(sql, "Parameters")$n.particles.per.cell[1]
n.new.particles.per.day <- dbReadTable(sql, "Parameters")$n.new.particles.per.day[1]
n.steps.per.day <- dbReadTable(sql, "Parameters")$n.hours.per.day[1]
dbDisconnect(sql)

## ------------------

particles.reference.bm.desc <- dget( paste0(project.folder,"/InternalProc/particles.reference.desc"))

## ------------------

cl.2 <- makeCluster(number.cores)
registerDoParallel(cl.2)

all.connectivity.pairs.to.sql <- foreach(cell.id.ref.f=cell.to.process, .verbose=FALSE, .combine = rbind ,  .packages=c("gstat","raster","data.table","FNN","bigmemory")) %dopar% { # 
  
  particles.reference.bm.i <- attach.big.matrix(particles.reference.bm.desc)
  connectivity.temp.m <- particles.reference.bm.i[ mwhich(particles.reference.bm.i,2,list(cell.id.ref.f), list('eq')) , ]
  connectivity.temp.m <- as.data.table(connectivity.temp.m)
  connectivity.temp.m <- connectivity.temp.m[connectivity.temp.m$cell.rafted != 0,]
  connectivity.temp.m <- data.frame(connectivity.temp.m,travel.time=(connectivity.temp.m$t.finish - connectivity.temp.m$t.start)/n.steps.per.day)
  connectivity.temp.m <- as.data.table(connectivity.temp.m)
  connectivity.pairs.to.sql <- data.frame()
  
  for(y in unique(connectivity.temp.m$start.year)) {
    
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

# Overall Connectivity matrix (mean of all years)

sql <- dbConnect(RSQLite::SQLite(), paste0(sql.directory,"/",project.name,"SimulationResults.sql"))
Connectivity <- data.table(dbReadTable(sql, "Connectivity"))
source.sink.xy <- dbReadTable(sql, "SourceSinkSites")
dbDisconnect(sql)

Connectivity <- Connectivity[ , j=list(mean(Probability, na.rm = TRUE) , max(Probability, na.rm = TRUE) , mean(Time.mean, na.rm = TRUE) , max(Time.mean, na.rm = TRUE) , mean(Number.events, na.rm = TRUE) ) , by = list(Pair.from,Pair.to)]
colnames(Connectivity) <- c("Pair.from" , "Pair.to" , "Mean.Probability" , "Max.Probability" , "Mean.Time" , "Max.Time" , "N.events" )
Connectivity

source.sink.id <- 1:nrow(source.sink.xy)

Connectivity.matrix.probability <- acast(Connectivity[,.(Pair.from , Pair.to ,  Mean.Probability)], Pair.from~Pair.to, value.var="Mean.Probability")
Connectivity.matrix.time <- acast(Connectivity[,.(Pair.from , Pair.to ,  Mean.Time)], Pair.from~Pair.to, value.var="Mean.Time")
dim(Connectivity.matrix.probability)

View(Connectivity.matrix.probability)
View(Connectivity.matrix.time)

# -----------------------------------------

# save(Connectivity,file=paste0(sql.directory,"/",project.name,"ConnectPairs.RData"))
# save(Connectivity.matrix.probability,file=paste0(sql.directory,"/",project.name,"ConnectMatrixProb.RData"))
# save(Connectivity.matrix.time,file=paste0(sql.directory,"/",project.name,"ConnectMatrixTime.RData"))
# save(source.sink.xy,file=paste0(sql.directory,"/",project.name,"SourceSinkXY.RData"))

## --------------------------------------------------------------------------------------------------------------
## --------------------------------------------------------------------------------------------------------------