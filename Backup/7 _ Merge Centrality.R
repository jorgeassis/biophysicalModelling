
library(igraph)


gist.directory <- "/Volumes/Jellyfish/Dropbox/Gist/One Aquarium V2.0/" # Laminaria Jellyfish
raw.data.dir <- "/Volumes/Jellyfish/Dropbox/Manuscripts/Genetic diversity drivers of Sargassum thunbergii/Data/"
results.directory <- "/Volumes/Jellyfish/Dropbox/Manuscripts/Genetic diversity drivers of Sargassum thunbergii/Results"
results.file <- "Sargassum"

study.region.file <- "Data/_ Asia HR 0.01/coast_line.tif"
ocean.region.file <- "Data/_ Asia HR 0.01/ocean.tif"

## ------------------------------------------------------------------------------------------------------------------------------
##
## ------------------------------------------------------------------------------------------------------------------------------

setwd(gist.directory)
source(paste0(gist.directory,"Dependences.R"))

sql.file <- list.files(path=paste0(results.directory,"/SQLite/") , pattern="reference" , full.names = TRUE)

study.region <- raster(paste0(gist.directory,study.region.file))
ocean.region <- raster(paste0(gist.directory,ocean.region.file))

## -------------------

sql <- dbConnect(RSQLite::SQLite(), sql.file)
cells <- dbReadTable(sql, "Cells_coordinates")
cells.id <- cells$cell
global.simulation.parameters <- dbReadTable(sql, "Simulation_parameters")
dbDisconnect(sql)

plot(cells[,2:3])

## -----------------------------------

sql.con <- dbConnect(RSQLite::SQLite(), sql.file)
reference.table <- data.table(dbReadTable(sql.con, "Particles_reference_table")) ; reference.table

Connectivity <- data.table(dbReadTable(sql.con, "Connectivity")) ; Connectivity
Connectivity <- Connectivity[ Pair.from %in% cells.id ,  ]
Connectivity <- Connectivity[ Pair.to %in% cells.id ,  ]
dbDisconnect(sql.con)

Connectivity <- Connectivity[ , j=list(mean(Probability, na.rm = TRUE) , max(Probability, na.rm = TRUE) , mean(Time.mean, na.rm = TRUE) , max(Time.mean, na.rm = TRUE) , mean(Number.events, na.rm = TRUE) ) , by = list(Pair.from,Pair.to)]
colnames(Connectivity) <- c("Pair.from" , "Pair.to" , "Mean.Probability" , "Max.Probability" , "Mean.Time" , "Max.Time" , "N.events" )
Connectivity

## -----------------------------------

comb <- Connectivity[Max.Time <= 30,.(Pair.from,Pair.to,Mean.Probability)]
comb <- comb[Pair.from != Pair.to,]
comb <- as.data.frame( comb[ sort(comb[,Mean.Probability] , decreasing = TRUE, index.return =TRUE)$ix , ] )

graph.obj <- graph.edgelist( cbind( as.character( comb[,1]) , as.character(comb[,2]) ) , directed = TRUE )
E(graph.obj)$weight = -log(comb[,3]) # Hock, Karlo Mumby, Peter J 2015
graph.obj <- simplify(graph.obj, remove.loops = TRUE , remove.multiple = TRUE)

sampling.sites <- "/Volumes/Jellyfish/Dropbox/Manuscripts/Genetic diversity drivers of Sargassum thunbergii/Data/sampling.sites.csv"
sampling.sites <- read.csv(sampling.sites,header = FALSE)[,2:3]

ocean.region <- raster(paste0(gist.directory,ocean.region.file))
plot(ocean.region,box=FALSE,legend=FALSE,col=c("black"))
points(sampling.sites,col="red")

position.matrix <- spDists(as.matrix(cells[,2:3]),as.matrix(sampling.sites),longlat = TRUE)
position.matrix <- apply(position.matrix,2,which.min)
position.matrix <- cells[position.matrix,1]
plot(cells[,2:3])
points(cells[position.matrix,2:3],col="green")

## -----------------------------------

centr <- alpha_centrality(graph.obj, nodes = as.character(position.matrix))
between <- betweenness(graph.obj, v = as.character(position.matrix), directed = TRUE, nobigint = TRUE)
close <- closeness(graph.obj, vids = as.character(position.matrix), mode = "in")
deg <- degree(graph.obj, v = as.character(position.matrix),loops = FALSE, normalized = FALSE)


library(gdata)

diversity <- read.xls("/Volumes/Jellyfish/Dropbox/Manuscripts/Genetic diversity drivers of Sargassum thunbergii/Results/Diversity.xlsx")

plot(centr,diversity$A)


