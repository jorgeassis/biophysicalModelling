## ------------------------------------------------------------------------------------------------------------------
## PlankTonic
## Assis et al., 2018
## ------------------------------------------------------------------------------------------------------------------
##
## ------------------------------------------------------------------------------------------------------------------

closeAllConnections()
rm(list=(ls()[ls()!="v"]))
gc(reset=TRUE)
source("0. Config.R")
source("Dependences/mainFunctions.R")

## --------------------------------------------------------------------------------------------------------------
##
##
## 
## --------------------------------------------------------------------------------------------------------------

pipeLiner <- TRUE
doParallelCalculations <- TRUE # repeat all parallel computations
type <- "points" # points polygons

load(paste0(results.folder,"/modelParameters.RData"))
n.days.max <- global.simulation.parameters$particle.max.duration
n.days.max

## ----------------------------------------------------------------

pld.period <- 1:n.days.max # c(10 , 30 , 90 , 120 , 200)
n.seasons <- "" # c("","Spring","Summer","Autumn","Winter")
combinations <- expand.grid(season=n.seasons,pld.period=pld.period,stringsAsFactors = F)

for( c in 1:nrow(combinations) ) {
  
  gc(reset=TRUE)
  cat(c,"\n")
  dev.off()
  
  season <- combinations[c,1]
  pld.period <- combinations[c,2]

  if( type == "polygons" ) { source("5.1. Connectivity Estimates [Polygons].R") }
  if( type == "points" ) { source("5.1. Connectivity Estimates [Points].R") }
  
}

## ----------------------------------------------------------------

mainTheme <- theme(panel.grid.major = element_blank() ,
                   text = element_text(size=12) ,
                   axis.title.y = element_text(margin = margin(t = 0, r = 12, b = 0, l = 0)) ,
                   axis.title.x = element_text(margin = margin(t = 12, r = 0, b = 0, l = 0)) )

names(combResults)

y.lab <- "Isolation degree (number of sites)"
y.lab <- "Degree centrality (average number of connections)"
y.lab <- "Degree centrality (maximum number of connections)"
y.lab <- "Number of clusters"

p3 <- ggplot() +
  geom_point(data = combResults, aes(x=pld, y=numberClusters), shape = 21,colour = "black", fill = "black", size = 2, stroke = 0.75, alpha = 0.5) +
  theme_minimal() + mainTheme + xlab("Dispersal period (day)") + ylab(y.lab) # + geom_smooth(data = combResults, aes(x=pld, y=numberClustersConsensus),span = 0.8)
p3

pdf( file=paste0(results.folder,"/",y.lab,".pdf"), width = 10, height = 8 )
print(p3)
dev.off()

## ----------------------------------------------------------------------------------------------------------
## ----------------------------------------------------------------------------------------------------------