
setwd("/media/nas/One Aquarium V1.0")
library(raster)
library(igraph)

raw.data.dir <- "Data/"
results.directory <- "Results/"
results.file <- "NEAtlantic"

## ------------------------------------------------------------------------------------------------------------------

source(paste0("Digital Aquarium Dependences.R"))
sql.file <- list.files(path=paste0(results.directory,"/SQLite/") , pattern="reference" , full.names = TRUE)

## ------------------------------------------------------------------------------------------------------------------
## Remove data

library(sqldf)

sub.query <- TRUE
sub.query.shapefile <- "Data/Unused2.shp"

if( sub.query ) {
  
  sub.query.shapefile <- shapefile(sub.query.shapefile)
  
  sql.con <- dbConnect(RSQLite::SQLite(), sql.file)
  cells <- dbReadTable(sql.con, "Cells_coordinates")
  dbDisconnect(sql.con)
  
  cells.sp <- cells
  coordinates(cells.sp) <- c("x", "y")
  proj4string(cells.sp) <- proj4string(sub.query.shapefile)
  over.unused <- !is.na(over(cells.sp, as(sub.query.shapefile, "SpatialPolygons")))
  
  final.cells <- which(!over.unused)
}

if( ! sub.query ) {
  
  sql.con <- dbConnect(RSQLite::SQLite(), sql.file)
  cells <- dbReadTable(sql.con, "Cells_coordinates")
  dbDisconnect(sql.con)
  
  final.cells <- cells$cell
}

final.cells

## -----------------------------------

sql.con <- dbConnect(RSQLite::SQLite(), sql.file)
reference.table <- data.table(dbReadTable(sql.con, "Particles_reference_table"))
reference.table <- reference.table[ cell %in% final.cells ,  ] ; nrow(reference.table)
dbDisconnect(sql.con)

sql.con <- dbConnect(RSQLite::SQLite(), sql.file)
Connectivity <- data.table(dbReadTable(sql.con, "Connectivity"))
Connectivity <- Connectivity[ Pair.from %in% final.cells ,  ]
Connectivity <- Connectivity[ Pair.to %in% final.cells ,  ]
dbDisconnect(sql.con)

Connectivity <- Connectivity[ , j=list(mean(Probability, na.rm = TRUE) , mean(Time.mean, na.rm = TRUE)) , by = list(Pair.from,Pair.to)]

# ---------------------------------------------------------
# Define Network for analyses

npops <- length(final.cells)
pop.names <- final.cells

network.type <- "Prob" # Time

if(network.type == "Prob") {
  comb <- Connectivity[,.(Pair.from,Pair.to,V1)]
  comb <- as.data.frame( comb[ sort(comb[,V1] , decreasing = TRUE, index.return =TRUE)$ix , ])
  net.function <- prod
  g2 <- graph.edgelist( cbind(as.character( comb[,1]) , as.character(comb[,2]) ) , directed = TRUE )
  E(g2)$weight = 1 - comb[,3] # The wheight as a negative impact on finding the closest path
}

if(network.type == "Time") {
  comb <- Connectivity[,.(Pair.from,Pair.to,V2)]
  comb <- as.data.frame( comb[ sort(comb[,V2] , decreasing = TRUE, index.return =TRUE)$ix , ])
  net.function <- sum
  g2 <- graph.edgelist( cbind(as.character( comb[,1]) , as.character(comb[,2]) ) , directed = TRUE )
  E(g2)$weight = comb[,3] # The wheight as a negative impact on finding the closest path
}

# ---------------------------------------------------------
# Function of stepping-stones

get.path <- function(i,j) {
  
  if( i != j ) {
    
    stones <- unlist( get.shortest.paths(g2,as.character(i),as.character(j),mode="out")$vpath[[1]] )
    
    if(sum(stones) == 0 ){
      
      try.out <- 1
      
      closest.i <- spDists(as.matrix(cells[,2:3]),as.matrix(cells[i,2:3]),longlat = TRUE)
      closest.j <- spDists(as.matrix(cells[,2:3]),as.matrix(cells[j,2:3]),longlat = TRUE)
      
      seq.try <- expand.grid(x=sort(c(closest.i),index.return=TRUE)$ix[2:6],y=sort(c(closest.j),index.return=TRUE)$ix[2:6])
      
      while( sum(stones) == 0 & try.out <= nrow(seq.try) ) {
        
        if(as.character(seq.try[try.out,1]) %in%  V(g2)$name & as.character(seq.try[try.out,2]) %in%  V(g2)$name) {
          
          stones <- unlist( get.shortest.paths(g2,as.character(seq.try[try.out,1]),as.character(seq.try[try.out,2]))$vpath[[1]] )
          
        }
        
        try.out <- try.out + 1
      }
      
    }
    
    if( sum(stones) == 0 ) {   return( c(i,j,0) ) }
    if( sum(stones) != 0 ) {   stones <- as.numeric(names(stones))
    
    m1 <- stones[-length(stones)]
    m2 <- stones[-1]
    stones.interm <- cbind(m1,m2)
    
    path.values <- apply( stones.interm , 1 , function(x) { comb[ comb[,1] == x[1] & comb[,2] == x[2] , 3 ] }   )
    connectivity.value <- apply( t(path.values) , 1 , net.function )
    
    return(c(i,j,connectivity.value))
    
    }  
    
    
  }
  
  if( i == j ) { return( c(i,j,1) ) }
  
}

# ---------------------------------------------------------
# ---------------------------------------------------------
# ---------------------------------------------------------
# Connectivity Matrix between Pairs

sampling.sites.shp <- "GIS/Shapefiles/by Species/29_Fucus_vesiculosus_msat.shp"
sampling.sites <- as.data.frame(shapefile(sampling.sites.shp))[,3:4]
sampling.sites.names <- as.data.frame(shapefile(sampling.sites.shp))[,5]

sampling.sites <- sampling.sites[-c(1,2),]
sampling.sites.names <- sampling.sites.names[-c(1,2)]

position.matrix <- spDists(as.matrix(cells[,2:3]),as.matrix(sampling.sites),longlat = TRUE)
position.matrix <- apply(position.matrix,2,which.min)
position.matrix <- cells[position.matrix,1]

norm <- expand.grid(from=position.matrix ,to=position.matrix ) 
norm.names <- expand.grid(from=sampling.sites.names,to=sampling.sites.names )

library(parallel)
library(doParallel)
library(reshape2)

cl.2 <- makeCluster(6) ; registerDoParallel(cl.2)
connectivity.stepping.stone <- foreach(x=1:nrow(norm), .verbose=FALSE, .packages=c("igraph","sp")) %dopar% { get.path( norm[x,1] , norm[x,2] )[3] } #
stopCluster(cl.2) ; rm(cl.2)

connectivity.stepping.stone <- data.frame(norm.names,probability=unlist(connectivity.stepping.stone))
connectivity.stepping.stone.matrix <- acast(connectivity.stepping.stone, from~to, value.var="probability")

write.table(connectivity.stepping.stone.matrix, file = "Results/connectivity_29_prob.txt", sep = "\t", row.names = FALSE, col.names = TRUE, na = "NA", dec = ",")
