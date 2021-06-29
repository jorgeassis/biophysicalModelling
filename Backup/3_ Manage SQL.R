## ------------------------------------------------------------------------------------------------------------------
## ------------------------------------------------------------------------------------------------------------------
##
## Digital Aquarium 2016 Versio 3.0
## Simulation of Dispersal using ocean fields
##

## -------------------------------------------------------------------
## Remove Years from Database (if needed) - Question for last year or other!

sql <- dbConnect(RSQLite::SQLite(), paste0(results.directory,"/SQLite/",results.files,".reference.particles.sql"))
dbListFields(sql, "Connectivity")
dbListFields(sql, "Particles_reference_table")
all.connectivity <- dbReadTable(sql, "Connectivity") ; unique(all.connectivity$Year)
all.particles.reference <- dbReadTable(sql, "Particles_reference_table") ; unique(all.particles.reference$Year)
dbExecute(sql, 'DELETE FROM Connectivity WHERE "Year" == 2006')
dbExecute(sql, 'DELETE FROM Particles_reference_table WHERE "Year" == 2006')
dbDisconnect(sql)