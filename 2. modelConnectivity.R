## ------------------------------------------------------------------------ ##
## ------------------------------------------------------------------------ ##
##                                                                          ##
##  Jorge Assis [ jorgemfa@gmail.com ]                                      ##
##  biodiversityDS                                                          ##
##  biodiversityDataScience.com                                             ##
##                                                                          ##
##  Biophysical modelling framework to predict oceanographic connectivity   ##
##  version: 2.00                                                           ##
##                                                                          ##
## ------------------------------------------------------------------------ ##
## ------------------------------------------------------------------------ ##

## Notes:
## 1. Numerical particles should be released along the actual path of the animal at regular, relatively short time intervals, i.e. between 6 h and 2 d, as that might give a better estimate of the current speed and direction. 

## ------------------------------------------------------------------ ##

closeAllConnections()
rm(list=(ls()))
gc(reset=TRUE)

# library(credentials)
# set_github_pat()

# exportRequirements()
# exportPackages()

source("0. Config.R")
source("Dependences/mainFunctions.R")

## ------------------------------------------------------------------ ##
## File structure

if( ! dir.exists(results.folder) ) { dir.create(file.path(results.folder), showWarnings = FALSE) } 
if( ! dir.exists(paste0(results.folder,"/Figures/")) ) { dir.create(file.path(paste0(results.folder,"/Figures/")), showWarnings = FALSE) } 

## ------------------------------------------------------------------ ##
## Define region of interest and hexagon list

source("Dependences/produceEnvironment.R")

plot1
pdf(file=paste0(results.folder,"/Figures/sourceSinkSitesHexagons.pdf"),width=12)
print(plot1)
dev.off()

## ------------------------------------------------------------------ ##

source("Dependences/listData.R")

# ------------------

source("Dependences/defineParticleConditions.R")

# ------------------

source("Dependences/produceReferenceTable.R")

# ------------------

source("Dependences/prepareAnimation.R")

plot2
pdf(file=paste0(results.folder,"/Figures/animationSites.pdf"),width=12)
print(plot2)
dev.off()

## ------------------------------------------------------------------ ##

source("Dependences/cleanEnvironment.R")

# save.image(paste0(results.folder,"/","environmentState.RData"))
# load(paste0(results.folder,"/","environmentState.RData"))

## ------------------------------------------------------------------ ##

## Run biophysical modelling

source("Dependences/biophysicalModel.R")

## ------------------------------------------------------------------ ##

## Export simulation parameters

particles.reference.bm <- attach.big.matrix(particles.reference.bm.desc)
if( length(mwhich(particles.reference.bm,c(9),list(0), list('eq'))) > 0 ) { print("An error may have occurred [!]")}

nonRafters <- length(mwhich(particles.reference.bm,c(12),list(0), list('eq')))
rafters <- length(mwhich(particles.reference.bm,c(12),list(0), list('neq')))

global.simulation.parameters <- data.frame(   project.name = project.name,
                                              sim.years = paste(c(from.year,to.year),collapse="-"),
                                              sim.months = paste(months.all,collapse=","),
                                              n.hours.per.day = n.hours.per.day , 
                                              n.new.particles.per.day = n.new.particles.per.day , 
                                              remove.new.particles.last.days = remove.new.particles.last.days , 
                                              longevity = longevity , 
                                              particle.max.duration = particle.max.duration , 
                                              n.particles.per.cell = n.particles.per.cell,
                                              movie.year = movie.year, 
                                              movie.sites.id = paste(movie.sites.id,collapse=",") , 
                                              sim.resolution = sim.resolution,
                                              extent = paste(c(min.lon,max.lon,min.lat,max.lat),collapse=","),
                                              particlesLost=nonRafters,
                                              particlesConnectEvents=rafters )       

global.simulation.parameters
save(global.simulation.parameters, file = paste0(results.folder,"/","modelParameters.RData"))

## -----------------
## -----------------

particles.reference.bm <- particles.reference.bm[]
particles.reference.bm <- particles.reference.bm[particles.reference.bm[,12] != 0, ]

file.remove( list.files(results.folder, full.names = TRUE, pattern = "particleReferenceTable.bin") )
file.remove( list.files(results.folder, full.names = TRUE, pattern = "particleReferenceTable.desc") )

particles.reference.bm <- as.big.matrix(particles.reference.bm , backingpath=paste0(results.folder,"/") , backingfile = "particleReferenceTable.bin", descriptorfile = "particleReferenceTable.desc")

## -----------------
## -----------------

if( file.exists(paste0(results.folder,"/particles.video.location.csv")) ) {
  
  particles.video.location.dt <- read.csv(paste0(results.folder,"/particles.video.location.csv"), sep=";")
  particles.video.location.dt <- as.data.table(particles.video.location.dt)
  save(particles.video.location.dt, file = paste0(results.folder,"/","particles.video.location.RData"))
  file.remove(paste0(results.folder,"/particles.video.location.csv"))
  
}

## ------------------------------------------------------------------------------------------------------------------
# Time taken

seq.t <- numeric(length(time.i)) ; seq.t[1] <- 0
for( t.seq in 2:length(seq.t)) {
  seq.t[t.seq] <- as.numeric(difftime(time.i[t.seq], time.i[t.seq-1], units = "mins"))
}
plot(1:length(seq.t),seq.t)

## ---------------------------------------------------------------------------------------------------------------------------------
## ---------------------------------------------------------------------------------------------------------------------------------
## End of Code