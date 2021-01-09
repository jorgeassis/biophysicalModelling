## ------------------------------------------------------------------------------------------------------------------
## PlankTonic
## Assis et al., 2018
## ------------------------------------------------------------------------------------------------------------------
##
## ------------------------------------------------------------------------------------------------------------------

rm(list=(ls()[ls()!="v"]))
gc(reset=TRUE)
source("../Project Config 0.R")
source("Dependences.R")

## --------------------------------------------------------------------------------------------------------------
##
##
## 
## --------------------------------------------------------------------------------------------------------------

pipeLiner <- TRUE
type <- "polygons" # coordinates polygons

load(paste0(project.folder,"/Results/",project.name,"/InternalProc/","Parameters.RData"))
n.days.max <- global.simulation.parameters$particle.max.duration
n.days.max

## ----------------------------------------------------------------

n.seasons <- c("","Spring","Summer","Autumn","Winter")
n.repetitions <- 1:120

for( n.days in n.repetitions ) {
  
  gc(reset=TRUE)
  n.season <- "Summer"
  
  if( type == "polygons" ) { source("5.1. Connectivity Estimates [Polygons].R") }
  if( type == "coordinates" ) { stop("Revise") ; source("5.1. Connectivity Estimates [Coord Sites].R") }
  
}

## --------------------------------------------------------------------------------------------------------------
## --------------------------------------------------------------------------------------------------------------