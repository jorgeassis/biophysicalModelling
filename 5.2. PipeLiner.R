## ------------------------------------------------------------------------------------------------------------------
## PlankTonic
## Assis et al., 2018
## ------------------------------------------------------------------------------------------------------------------
##
## ------------------------------------------------------------------------------------------------------------------

rm( list=(ls()[ls()!="v"]) )
gc(reset=TRUE)

source("0. Project Config.R")
source("Dependences.R")

## --------------------------------------------------------------------------------------------------------------
##
##
## 
## --------------------------------------------------------------------------------------------------------------

pipeLiner <- TRUE
type <- "points" # points polygons

load(paste0(project.folder,"/Results/",project.name,"/InternalProc/","Parameters.RData"))
n.days.max <- global.simulation.parameters$particle.max.duration
n.days.max

## ----------------------------------------------------------------

list.dirs(path = paste0("../Results"), recursive = FALSE)

popCoordinates <- read.csv("../Data/DiffRecords.csv", sep=";")
popFST <- read.csv("../Data/DiffFST.csv", sep=";")

pld.period <- 1:60 # 1:120 c(10 , 30 , 90 , 120 , 200)
n.seasons <- "YearRound" # c("YearRound","Spring","Summer","Autumn","Winter")
combinations <- expand.grid(season=n.seasons,pld.period=pld.period,stringsAsFactors = F)

for( c in 1:nrow(combinations) ) {
  
  gc(reset=TRUE)
  
  season <- combinations[c,1]
  pld.period <- combinations[c,2]

  if( season == "Spring" ) { spawn.p <- c(3,4,5) }
  if( season == "Summer" ) { spawn.p <- c(6,7,8) }
  if( season == "Autumn" ) { spawn.p <- c(9,10,11) }
  if( season == "Winter" ) { spawn.p <- c(12,1,2) }
  if( season == "YearRound" ) { spawn.p <- 1:12 }
  
  if( type == "polygons" ) { source("5.2. Connectivity vs Differentiaiton.R") }
  if( type == "points" ) { source("5.2. Connectivity vs Differentiaiton.R") }
  
}

## --------------------------------------------------------------------------------------------------------------
## --------------------------------------------------------------------------------------------------------------