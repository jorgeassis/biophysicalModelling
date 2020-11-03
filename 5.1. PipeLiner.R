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
n.days.max <- global.simulation.parameters
n.days.max

## ----------------------------------------------------------------


if( type == "polygons") {
  
  isolatedAgg <- data.frame()
  betweennessAgg <- data.frame()
  higherBetweennessAgg <- data.frame()
  eighenCentralityAgg <- data.frame()
  highereighenCentralityAgg <- data.frame()
  closenessAgg <- data.frame()
  higherclosenessAgg <- data.frame()
  clusterAssignmentAgg <- data.frame()
  resistanceAgg <- data.frame()
  higherresistanceAgg <- data.frame()
  
}

n.repetitions <- c(1,5,30,60,120)

for( n.days in n.repetitions ) {
  
  if( type == "polygons" ) { source("5.1. Connectivity Estimates [Polygons].R") }
  if( type == "coordinates" ) { stop("Revise") ; source("5.1. Connectivity Estimates [Coord Sites].R") }
  
}

