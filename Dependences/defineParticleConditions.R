## ------------------------------------------------------------------------ ##
## ------------------------------------------------------------------------ ##

## Define conditions

norm.time <- data.frame()
for(n in 1:nrow(simulation.parameters.step)) {
  norm.time <- rbind(norm.time,data.frame(do.call("rbind", replicate(n.hours.per.day, simulation.parameters.step[n,c("year" , "month" , "day")] , simplify = FALSE)) , hour=1:n.hours.per.day))
}
norm.time <- data.table(norm.time)

# ------------------

release.particles.condition <- round(seq(1,n.hours.per.day,n.hours.per.day/n.new.particles.per.day))

norm.time[ hour == release.particles.condition, release.particles := TRUE ]
norm.time[ hour != release.particles.condition, release.particles := FALSE ]

if(remove.new.particles.last.days) {
  norm.time[ ( nrow(norm.time)-(remove.new.particles.last.days.n.days*n.hours.per.day)):nrow(norm.time) , release.particles := FALSE ]
}

n.particles.per.cell <- (nrow(simulation.parameters.step)) * n.new.particles.per.day

cat("Particle conditions defined:", "TRUE","\n")
cat("Number of time steps:", nrow(simulation.parameters.step),"\n")
cat("Number of particles per day:", n.new.particles.per.day,"\n")
cat("Number particles per source site:", n.particles.per.cell,"\n")
