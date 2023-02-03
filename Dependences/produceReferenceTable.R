## ------------------------------------------------------------------------ ##
## ------------------------------------------------------------------------ ##

## Define reference table
# 0 unborne; 1 living; 2 rafted; 3 out of space; 4 dead by time

file.remove( list.files(paste0(results.folder,"/"), full.names = TRUE, pattern = "particleReferenceTable") ) 
particles.reference.names <- c("id","start.cell","start.year","start.month","start.day","pos.lon","pos.lat","pos.alt","state","t.start","t.finish","cell.rafted","at.sea","travel.time")
save(particles.reference.names,file=paste0(results.folder,"/particleReferenceTableNames.RData"))

particles.reference.bm <- big.matrix(nrow=(n.particles.per.cell * nrow(initial.coords)),ncol=length(particles.reference.names) , backingpath=paste0(results.folder,"/") , backingfile = "particleReferenceTable.bin", descriptorfile = "particleReferenceTable.desc")
particles.reference.bm.desc <- dget( paste0(results.folder,"/particleReferenceTable.desc") )
particles.reference.bm <- attach.big.matrix(particles.reference.bm.desc)

if( nrow(initial.coords) > 999999 ) { stop("Error :: 935")  }
 
for( i in 1:nrow(simulation.parameters.step) ) {
  
  loc.i <- (( sum(source.sink.xy$source == 1) * i + 1 ) - sum(source.sink.xy$source == 1) ):(( sum(source.sink.xy$source == 1) * i - 1 ) + 1)
  
  particles.reference.bm[loc.i,1] <- loc.i
  particles.reference.bm[loc.i,2] <- source.sink.xy[source.sink.xy$source == 1 , "cells.id" ]
  particles.reference.bm[loc.i,6] <- source.sink.xy[source.sink.xy$source == 1 , "x" ]
  particles.reference.bm[loc.i,7] <- source.sink.xy[source.sink.xy$source == 1 , "y" ]
  
}

cat("Reference table defined:", "TRUE","\n")
cat("Dimensions of reference table:", dim(particles.reference.bm)[1],"x",dim(particles.reference.bm)[2],"\n")

particles.reference.template <- data.table(particles.reference.bm[])
names(particles.reference.template) <- particles.reference.names

particles.reference.template.id <- particles.reference.template[ state == 0 , .SD[1] , by=start.cell][,id]
particles.reference.template <- particles.reference.template[particles.reference.template.id, ]
particles.reference.template$id <- 0
rm(particles.reference.bm)