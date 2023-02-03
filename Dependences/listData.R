## ------------------------------------------------------------------------ ##
## ------------------------------------------------------------------------ ##

raw.data.files <- list.files(data.folder,full.names = TRUE,pattern="nc")
raw.data.files <- raw.data.files[unlist(sapply(from.year:to.year,function(x){ which(grepl(x,raw.data.files)) }))]
simulation.parameters.step <- data.frame()
simulation.dimension.test <- data.frame()

cat("\n","Available raw data files:", length(raw.data.files),"\n")

for( file in 1:length(raw.data.files)) {
  nc <- nc_open( raw.data.files[file] , verbose=FALSE )
  nc.date <- as.Date(ncvar_get( nc, "Date"), origin = as.Date("1970-01-01"))
  simulation.parameters.step <- rbind( simulation.parameters.step, data.frame( simulation=file, file=raw.data.files[file],  year=substr(nc.date, 1, 4) , month=substr(nc.date, 6, 7) , day=substr(nc.date, 9, 10) , stringsAsFactors = FALSE) )
  simulation.dimension.test <- rbind(simulation.dimension.test,data.frame(dim.i=length(ncvar_get( nc, "X")),dim.j=length(ncvar_get( nc, "Y"))))
  nc_close( nc )
}
simulation.parameters.step <- simulation.parameters.step[simulation.parameters.step$year %in% as.character(from.year:to.year) & simulation.parameters.step$day %in% sapply(from.day:to.day,function(x){ ifelse(nchar(x) > 1,x,paste0("0",x))}) & simulation.parameters.step$month %in% sapply(months.all,function(x){ ifelse(nchar(x) > 1,x,paste0("0",x))}) ,  ]

if(sum( ! from.year:to.year %in% as.numeric(simulation.parameters.step[,"year"])) + sum( ! months.all %in% as.numeric(simulation.parameters.step[,"month"])) > 0 ) { stop("Data is not available for time window") }
if( length(unique(simulation.dimension.test$dim.i)) > 1 | length(unique(simulation.dimension.test$dim.j)) > 1 ) { stop("Error :: 992")}

