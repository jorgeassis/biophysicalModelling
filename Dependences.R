## ------------------------------------------------------------------------------------------------------------------
## PlankTonic
## Assis et al., 2018
## ------------------------------------------------------------------------------------------------------------------

if( ! "Results" %in% list.files(project.folder) ) { dir.create(file.path(paste0(project.folder,"/Results"))) }
if( ! "Data" %in% list.files(project.folder) ) { dir.create(file.path(paste0(project.folder,"/Data"))) }
if( ! "SQL" %in% list.files(paste0(project.folder,"/Results")) ) { dir.create(file.path(paste0(project.folder,"/Results/SQL"))) }
if( ! "InternalProc" %in% list.files(project.folder) ) { dir.create(file.path(paste0(project.folder,"/InternalProc"))) }

sql.directory <<- paste0(project.folder,"/Results/SQL/")

## -------------------------

packages.to.use <- c("gdata",
                     "gstat",
                     "compiler",
                     "data.table",
                     "raster",
                     "rgdal",
                     "ncdf4",
                     "parallel",
                     "doParallel",
                     "rgeos",
                     "FNN",
                     "sqldf",
                     "igraph",
                     "reshape2",
                     "gdistance",
                     "ggplot2",
                     "bigmemory",
                     "dismo",
                     "randomcoloR"
                     )

for(package in packages.to.use) {
  sink("/dev/null") 
      if( ! package %in% rownames(installed.packages()) ) { install.packages( package ) }
      if( ! package %in% rownames(installed.packages()) ) { sink() ; stop("Error on package instalation") }
      library(package, character.only = TRUE)
  sink()
}

## ---------------------------------------------------------------------------------------------------------------------
## ---------------------------------------------------------------------------------------------------------------------
## Functions

## Clean Dump

clean.dump.files <- function(clean.dump.files,files,dump.folder) {
  
  if( clean.dump.files ) { 
    file.remove( list.files(dump.folder, full.names = TRUE, pattern = files) ) 
  }
  
}

## ---------------------------------------------------------------------------------------------------------------------

list.memory <- cmpfun( function (pos = 1, pattern, order.by = "Size", decreasing=TRUE, head = TRUE, n = 100) {
  
  napply <- function(names, fn) sapply(names, function(x)
    fn(get(x, pos = pos)))
  names <- ls(pos = pos, pattern = pattern)
  obj.class <- napply(names, function(x) as.character(class(x))[1])
  obj.mode <- napply(names, mode)
  obj.type <- ifelse(is.na(obj.class), obj.mode, obj.class)
  obj.size <- napply(names, object.size) / 10^6 # megabytes
  obj.dim <- t(napply(names, function(x)
    as.numeric(dim(x))[1:2]))
  vec <- is.na(obj.dim)[, 1] & (obj.type != "function")
  obj.dim[vec, 1] <- napply(names, length)[vec]
  out <- data.frame(obj.type, obj.size, obj.dim)
  names(out) <- c("Type", "Size", "Rows", "Columns")
  out <- out[order(out[[order.by]], decreasing=decreasing), ]
  if (head)
    out <- head(out, n)
  out
  
} )

## ---------------------------------------------------------------------------------------------------------------------

errors <- cmpfun( function(t.step.test) {
  

} )

## ---------------------------------------------------------------------------------------------------------------------

monitor.processes <- cmpfun( function (process.name) {  
  if( !exists("ptm") ) { ptm <- proc.time() }
  t.step.time <- proc.time() - ptm
  ptm <- proc.time()
  
  if( !exists("day") ) { day <- 0 }

  cat("P:" ,process.name, " | Day: ",day, " | Time Taken: " , round(t.step.time[3]) , " | Memory: ", round(sum(list.memory()$Size)) , "\n") 
  flush.console()
  ptm <<- proc.time()
} )           


## ---------------------------------------------------------------------------------------------------------------------

produce.network <- function(network.type,connectivity,extract.simulation.days,crop.network,buffer,cells,new.extent) {
  
  Connectivity <- as.data.table(connectivity)
  
  if(crop.network) {  final.cells <- which(   cells[,2] >= (new.extent[1] - buffer) & 
                                              cells[,2] <= (new.extent[2] + buffer) & 
                                              cells[,3] >= (new.extent[3] - buffer) & 
                                              cells[,3] <= (new.extent[4] + buffer) )   
  
  final.cells <- cells[final.cells,1]
  print(plot(cells[final.cells,2:3]))
  
  }
  
  if( ! crop.network ) {  final.cells <- cells[,1]  }
  
  if( network.type == "Prob" ) {
    
    comb <- Connectivity[Time.max <= extract.simulation.days,.(Pair.from,Pair.to,Probability)]
    comb <- comb[Pair.from %in% final.cells,]
    comb <- comb[Pair.to %in% final.cells,]
    comb <- comb[Pair.from != Pair.to,]
    comb <- as.data.frame( comb[ sort(comb[,Probability] , decreasing = TRUE, index.return =TRUE)$ix , ] )
    
    # norm <- t(combn(position.matrix, 2))
    # 
    # for( i in 1:nrow(norm)) {
    #   
    #   t.1 <- which( comb[,1] == norm[i,1] & comb[,2] == norm[i,2] )
    #   
    #   if( length(t.1) == 0 ) { comb <- rbind(comb,data.frame(Pair.from = norm[i,1] , Pair.to = norm[i,2] ,  Mean.Probability = 0)) }
    #   
    #   t.2 <- which( comb[,1] == norm[i,2] & comb[,2] == norm[i,1] )
    #   
    #   if( length(t.2) == 0 ) { comb <- rbind(comb,data.frame(Pair.from = norm[i,2] , Pair.to = norm[i,1] ,  Mean.Probability = 0)) }
    #   
    # }
    
    net.function <<- prod
    graph.obj <- graph.edgelist( cbind( as.character( comb[,1]) , as.character(comb[,2]) ) , directed = TRUE )
    # E(graph.obj)$weight = 1 - comb[,3] # The wheight has a negative impact on finding the closest path
    E(graph.obj)$weight = -log(comb[,3]) # Hock, Karlo Mumby, Peter J 2015
    graph.obj <- simplify(graph.obj, remove.loops = TRUE , remove.multiple = TRUE)
    
  }
  
  if( network.type == "Time" ) {
    
    comb <- Connectivity[Time.max <= extract.simulation.days,.(Pair.from,Pair.to,Time.mean)]
    comb <- comb[Pair.from %in% final.cells,]
    comb <- comb[Pair.to %in% final.cells,]
    comb <- comb[Pair.from != Pair.to,]
    comb <- as.data.frame( comb[ sort(comb[,Time.mean] , decreasing = TRUE, index.return =TRUE)$ix , ] )
    
    # norm <- t(combn(position.matrix, 2))
    # 
    # for( i in 1:nrow(norm)) {
    #   
    #   t.1 <- which( comb[,1] == norm[i,1] & comb[,2] == norm[i,2] )
    #   
    #   if( length(t.1) == 0 ) { comb <- rbind(comb,data.frame(Pair.from = norm[i,1] , Pair.to = norm[i,2] ,  Mean.Time = 9e9999)) }
    #   
    #   t.2 <- which( comb[,1] == norm[i,2] & comb[,2] == norm[i,1] )
    #   
    #   if( length(t.2) == 0 ) { comb <- rbind(comb,data.frame(Pair.from = norm[i,2] , Pair.to = norm[i,1] ,  Mean.Time = 9e9999)) }
    #   
    #   
    # }
    
    net.function <<- sum
    graph.obj <- graph.edgelist( cbind( as.character( comb[,1]) , as.character(comb[,2]) ) , directed = TRUE )
    E(graph.obj)$weight = comb[,3]
    graph.obj <- simplify(graph.obj, remove.loops = TRUE , remove.multiple = TRUE)
    
  }
  
  return(list(comb,graph.obj))
  
}

