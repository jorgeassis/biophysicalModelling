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

packages.to.use <- c("dggridR","gdata","dplyr","sf",
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

list.memory.used <- function(type) {
  
  if(type==1) { return(sum(gc()[,7]) ) }
  if(type==2) { return(sum(list.memory()$Size)) }
  
}

## ---------------------------------------------------------------------------------------------------------------------

list.memory <- cmpfun( function (pos = 1, pattern, order.by = "Size", decreasing=TRUE, head = TRUE, n = 100) {
  
  napply <- function(names, fn) sapply(names, function(x) fn(get(x, pos = pos)))
  names <- ls(pos = pos, pattern = pattern)
  obj.class <- napply(names, function(x) as.character(class(x))[1])
  obj.mode <- napply(names, mode)
  obj.type <- ifelse(is.na(obj.class), obj.mode, obj.class)
  obj.size <- napply(names, object.size) / 10^6 # megabytes
  obj.dim <- t(napply(names, function(x) as.numeric(dim(x))[1:2]))
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

trim.by.distance <- function(xyDF,source.sink.dist,parallel) {
  
  coastline.pts.t <- xyDF

  if(parallel){
    
    seqListing <- round(seq(min(coastline.pts.t$y),max(coastline.pts.t$y),length.out =number.cores))
    parallelChunks <- data.frame(from = c(min(coastline.pts.t$y),seqListing[-c(1,length(seqListing))]), to = c(seqListing[-c(1,length(seqListing))] , max(coastline.pts.t$y) ) )
    
    cl.2 <- makeCluster(number.cores)
    registerDoParallel(cl.2)
    
    source.sink.xy.t <- foreach(section=1:nrow(parallelChunks), .combine=rbind, .verbose=FALSE, .packages=c("dismo","gstat","gdata","raster","data.table","bigmemory","FNN")) %dopar% { 
    
      coastline.pts.t.sec <- coastline.pts.t[coastline.pts.t$y >= parallelChunks[section,1] & coastline.pts.t$y <= parallelChunks[section,2] ,]
    
      source.sink.xy.t <- data.frame(matrix(NA,ncol=2,nrow=nrow(coastline.pts.t.sec)))
      
      iteraction <- 0
      
      while( nrow(coastline.pts.t.sec) > 0 ){
        
        iteraction <- iteraction + 1
        
        pt.i = coastline.pts.t.sec[1,,drop=FALSE]
        
        source.sink.xy.t[iteraction,] <- as.data.frame(pt.i)
        
        coastline.pts.t.i <- coastline.pts.t.sec
        colnames(coastline.pts.t.i) <- c("x","y")
        coordinates(coastline.pts.t.i) <- c("x","y")
        crs(coastline.pts.t.i) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
        circle <- circles(pt.i, lonlat=TRUE, d=source.sink.dist*1000, dissolve=FALSE)
        circle <- geometry(circle)
        crs(circle) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
        
        to.extract <- which(!is.na(over(coastline.pts.t.i,circle)))
        coastline.pts.t.sec <- coastline.pts.t.sec[-to.extract,]
        
      }
      
      source.sink.xy.t <- source.sink.xy.t[complete.cases(source.sink.xy.t),]
      return(source.sink.xy.t)

    }
      
    stopCluster(cl.2) ; rm(cl.2)
    
    to.remove <- which(duplicated(source.sink.xy.t))
    if(length(to.remove) > 0) { source.sink.xy.t <- source.sink.xy.t[-to.remove,]}
    
  }
  
  if( ! parallel){
    
  source.sink.xy.t <- data.frame(matrix(NA,ncol=2,nrow=nrow(coastline.pts.t)))
    
  iteractions <- nrow(coastline.pts.t)
  iteraction <- 0
  
  while( nrow(coastline.pts.t) > 0 ){
    
    iteraction <- iteraction + 1
    progress.percent <- 100 - round((nrow(coastline.pts.t) / iteractions) * 100)
    
    cat('\014')
    cat('\n')
    
    cat('\n',paste0(rep("-",100),collapse = ""))
    cat('\n',paste0(rep("-",progress.percent),collapse = ""),"||",progress.percent,"% (",iteraction,")")
    cat('\n',paste0(rep("-",100),collapse = ""))
    
    pt.i = coastline.pts.t[1,,drop=FALSE]
    
    source.sink.xy.t[iteraction,] <- as.data.frame(pt.i)
    
    coastline.pts.t.i <- coastline.pts.t
    colnames(coastline.pts.t.i) <- c("x","y")
    coordinates(coastline.pts.t.i) <- c("x","y")
    crs(coastline.pts.t.i) <- dt.projection
    circle <- circles(pt.i, lonlat=TRUE, d=source.sink.dist*1000, dissolve=FALSE)
    circle <- geometry(circle)
    crs(circle) <- dt.projection
    
    to.extract <- which(!is.na(over(coastline.pts.t.i,circle)))
    coastline.pts.t <- coastline.pts.t[-to.extract,]
    
  }
  
  source.sink.xy.t <- source.sink.xy.t[complete.cases(source.sink.xy.t),]
  
  }
  
  return(source.sink.xy.t)
  
}

## ---------------------------------------------------------------------------------------------------------------------

produce.network <- function(network.type,comb,n.days,crop.network,buffer,cells,new.extent) {
  
  if(crop.network) {  final.cells <- which(   cells[,2] >= (new.extent[1] - buffer) & 
                                              cells[,2] <= (new.extent[2] + buffer) & 
                                              cells[,3] >= (new.extent[3] - buffer) & 
                                              cells[,3] <= (new.extent[4] + buffer) )   
  
  final.cells <- cells[final.cells,1]
  plot(cells[final.cells,2:3])
  
  }
  
  if( ! crop.network ) {  final.cells <- cells[,1]  }
  
  comb <- comb[Time.max <= n.days,]
  comb <- comb[,.(Pair.from,Pair.to,Probability)]
  comb <- comb[Pair.from %in% as.vector(unlist(final.cells)) & Pair.to %in% as.vector(unlist(final.cells)) ,]
  comb <- comb[Pair.from != Pair.to,]
  comb <- as.data.frame( comb[ sort(comb[,Probability] , decreasing = TRUE, index.return =TRUE)$ix , ] )
  
  if( network.type == "Prob" ) {
    
    net.function <<- prod
    graph.obj <- graph.edgelist( cbind( as.character( comb[,1]) , as.character(comb[,2]) ) , directed = TRUE )
    
    # E(graph.obj)$weight = 1 - comb[,3] # The wheight has a negative impact on finding the closest path
    E(graph.obj)$weight = -log(comb[,3]) # Hock, Karlo Mumby, Peter J 2015
    
    graph.obj <- simplify(graph.obj, remove.loops = TRUE , remove.multiple = TRUE)
    
  }
  
  if( network.type == "Time" ) {
    
    net.function <<- sum
    graph.obj <- graph.edgelist( cbind( as.character( comb[,1]) , as.character(comb[,2]) ) , directed = TRUE )
    E(graph.obj)$weight = comb[,3]
    graph.obj <- simplify(graph.obj, remove.loops = TRUE , remove.multiple = TRUE)
    
  }
  
  return(list(comb,graph.obj))
  
}

## ---------------------------------------------------------------------------------------------------------------------

padlock <- function(dir,fun,id) {
  
  if( missing(id) ) { id <- NULL } 
  if( missing(dir) ) { dir <- NULL } 
  if( missing(fun) ) { fun <- NULL } 
  
  options(warn=-1)
  
  file.t <- paste0("Padlocker#",id,".Lk")
  
  ## ---------------------------
  
  if( fun == "LockAndHaltOn" ) {
  
    int <- 0
    
    repeat {
      
      if( ! padlock(dir,"isLocked") ) { 
        
        int <- int + 1
        
        padlock(dir,"Lock",id=id)
        
        Sys.sleep(sample(seq(1,2,length.out = 10),1))

        if( ! padlock(dir,"uniqueLocker",id=id) ) { padlock(dir,"Unlock",id=id) ; next }
        
        if( padlock(dir,"uniqueLocker",id=id) ) { break }
        
      }
      
      if( int > 9999 ) { stop("More than 9999 tries") }
      
  }
    
  }
  
  ## ---------------------------
    
  if( fun == "Lock" ) {
    
    write("Locked",file=paste0(dir,"/",file.t),append=FALSE)
    
  }
  
  ## ---------------------------
  
  if( fun == "Unlock" ) {
    
    file.remove( paste0(dir,"/",file.t) )
    
    if( id == -1 ) {
      
      file.remove( list.files(dir,"Padlocker#",full.names = TRUE) )
      
    }
  }
  
  ## ---------------------------
  
  if( fun == "isLocked" ) {
    
    Locked <- TRUE
    
    fs <- list.files(dir,"Padlocker#",full.names = TRUE)
    
    if( length(fs) != 0 ) {
      Locked <- TRUE
    } 
    if( length(fs) == 0 ) {
      Locked <- FALSE
    } 
    
    return( Locked  )
    
  }
  
  ## ---------------------------
  
  if( fun == "uniqueLocker" ) {
    
    u.locker <- FALSE
    
    fs <- list.files(dir,"Padlocker#",full.names = TRUE)
    
    if( length(fs) == 0) { u.locker <- FALSE  }
    if( length(fs) > 1) { u.locker <- FALSE  }
    if( length(fs) == 1 ) { 
      
      if( grepl(file.t,fs) ) { u.locker <- TRUE  }  
      
    }
    
    return( u.locker  )
    
  }
  
  ## ---------------------------
  
  if( fun == "countLockersAge" ) {
    
    fs <- list.files(dir,"Padlocker#",full.names = TRUE)
    age <- as.numeric( difftime(Sys.time(), file.info(fs )$atime, units ="mins") )
    return( data.frame(fs,age)  )
    
  }
  
  ## ---------------------------
  
  if( fun == "Age" ) {
    
    fs <- list.files(dir,"Padlocker#",full.names = TRUE)
    age.i <- numeric(0)
    
    if( length(fs) > 0 ) {
      
      for( i in 1:length(fs)) {
        
        age <- as.numeric( difftime(Sys.time(), file.info(fs[i] )$atime, units ="mins") )
        age.i <- c(age.i,age)
        
      }
    }
    if( length(fs) == 0 ) {
      age.i <- 0
    }
    
    return(min(age.i))
    
  }
  
  ## ---------------------------
  
  options(warn=0)
  
}

## ---------------------------------------------------------------------------------------------------------------------
## ---------------------------------------------------------------------------------------------------------------------
