## ------------------------------------------------------------------------------------------------------------------
##
## Digital Aquarium 2015
## Simulation of Dispersal using ocean fields
##
## Assis, et al. 2015
##
## ------------------------------------------------------------------------------------------------------------------
## ------------------------------------------------------------------------------------------------------------------

packages.to.use <- c("gdata",
                     "compiler",
                     "data.table",
                     "raster",
                     "ncdf4",
                     "parallel",
                     "doParallel",
                     "rgeos",
                     "FNN",
                     "sqldf",
                     "igraph",
                     "reshape2",
                     "gdistance",
                     "ggplot2")

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



