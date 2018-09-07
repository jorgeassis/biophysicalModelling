## ------------------------------------------------------------------------------------------------------------------
## PlankTonic
## Assis et al., 2018
## ------------------------------------------------------------------------------------------------------------------

source("Dependences.R")

## ------------------------------------------------------------------------------------------------------------------------------
##
## ------------------------------------------------------------------------------------------------------------------------------

# Resolve connectivity

sql <- dbConnect(RSQLite::SQLite(), paste0(sql.directory,"/",project.name,"SimulationResults.sql"))
cell.to.process <- 1:nrow(dbReadTable(sql, "ReleaseSites"))
particles.reference <- as.data.table(dbReadTable(sql, "ReferenceTable"))
n.particles.per.cell <- dbReadTable(sql, "Parameters")$n.particles.per.cell
dbDisconnect(sql)

cl.2 <- makeCluster(number.cores)
registerDoParallel(cl.2)

all.connectivity.pairs.to.sql <- foreach(cell.id.ref.f=cell.to.process, .verbose=FALSE, .combine = rbind ,  .packages=c("gstat","raster","data.table","FNN")) %dopar% { # 
  
  connectivity.pairs.to.sql <- data.frame()
  connectivity.temp.m <- particles.reference[ start.cell == cell.id.ref.f , ]
  
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
                                                          Probability = nrow(connectivity.temp[ cell.rafted == cell.id.ref.t,]) / n.particles.per.cell,
                                                          Year = simulation.year ) )
        }
    
  }
  
  connectivity.pairs.to.sql[is.na(connectivity.pairs.to.sql)] <- 0
  return( connectivity.pairs.to.sql )
}

stopCluster(cl.2) ; rm(cl.2)

# -----------------------------------------

# Save pairs to SQL 

sql <- dbConnect(RSQLite::SQLite(), paste0(sql.directory,"/",project.name,"SimulationResults.sql"))
dbWriteTable(sql, "Connectivity", all.connectivity.pairs.to.sql , overwrite=TRUE, append=FALSE)
dbDisconnect(sql)

# -----------------------------------------

sql <- dbConnect(RSQLite::SQLite(), paste0(sql.directory,"/",project.name,"SimulationResults.sql"))
Connectivity <- data.table(dbReadTable(sql, "Connectivity")) ; Connectivity
dbDisconnect(sql)

Connectivity <- Connectivity[ , j=list(mean(Probability, na.rm = TRUE) , max(Probability, na.rm = TRUE) , mean(Time.mean, na.rm = TRUE) , max(Time.mean, na.rm = TRUE) , mean(Number.events, na.rm = TRUE) ) , by = list(Pair.from,Pair.to)]
colnames(Connectivity) <- c("Pair.from" , "Pair.to" , "Mean.Probability" , "Max.Probability" , "Mean.Time" , "Max.Time" , "N.events" )
Connectivity



## ------------------------------------------------------------------------------------------------------------------------------
## ------------------------------------------------------------------------------------------------------------------------------
