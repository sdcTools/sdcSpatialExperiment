#######################################################################
# Script to apply sdcSpatial on 500/1000/2000 m grid count tables
# 
# Script applies the following steps (in that order)
# 
# - load needed libraries
# - load dummy population data
# - define "raster"-object for grid cells to be used with sdcSpatial
# - define sdc_raster object using population data and raster-object from before
# - check sensitive cells
# - apply protection method; check sensitive cells again and suppress if necessary
# - get protected table and apply information loss function
# - check if the original value of protected cells can be re-estimated
#

##################################
# load libraries
library(data.table)
library(raster)
library(sf)
library(sdcSpatial)


##################################
# load population data
load("test_data.RData")
pop_data
# data.table should be one row per person and contain variables
# coordinates for each person (x, y)
# IDs for gridcells
#           - L000500 (500m)
#           - L001000 (1000m)
#           - L002000 (2000m)
# PID (optional) 

raster_cells
# data.table containing all raster cells for a country even unpopulated ones
# should contain all grid cell IDs
#           - L000500 (500m)
#           - L001000 (1000m)
#           - L002000 (2000m)

# define cell size
m <- 500 # dimension of rasters in meters
geo <-  paste(c("L",paste(rep(0,6-nchar(m))),m),collapse="")

# make table original cell counts
grid_data <- pop_data[,.(count_original=.N),by=c(geo)]
grid_data <- merge(grid_data, raster_cells, by=c(geo),all=TRUE)

# get center coordinates of grid cells - needed later
grid_data[,y:=as.numeric(paste0(substr(get(geo),6,10),"00"))]
grid_data[,x:=as.numeric(paste0(substr(get(geo),12,16),"00"))]
grid_data[is.na(count_original),count_original:=0]

##################################
# define "raster"-object for grid cells to be used with sdcSpatial

# download and read in shapefiles
path <- paste0("https://data.statistik.gv.at/data/OGDEXT_RASTER_1_STATISTIK_AUSTRIA_",
               geo,"_LAEA.zip")
layer <- paste0("STATISTIK_AUSTRIA_",
                geo,"_LAEA")

zipname <- basename(path)
# download .zip
my_shapefile <- paste0(getwd(),"/",zipname)
download.file(path,my_shapefile)
# unzip file
unzipfolder <- paste0(sub(".zip","",my_shapefile))
unzip(my_shapefile,exdir=unzipfolder)
# read shape file
area_shape <- sf::st_read(dsn = unzipfolder, layer = layer)
# delete zip file
file.remove(my_shapefile)
# delete folder optional
# unlink(unzipfolder, recursive = TRUE)

# create raster object to be used by sdcSpatial
bounds <- sf::st_bbox(area_shape)
ncol <- (bounds$xmax-bounds$xmin)/m+1
nrow <- (bounds$ymax-bounds$ymin)/m+1
projection <-  sf::st_crs(area_shape)$proj4string
r <- raster::raster(nrows=nrow,ncols=ncol,
                    xmn = bounds$xmin-m/2,
                    xmx = bounds$xmax+m/2,
                    ymn = bounds$ymin-m/2,
                    ymx = bounds$ymax+m/2,
                    crs = sp::CRS(projection))

# for overview on projections see:
# https://www.nceas.ucsb.edu/sites/default/files/2020-04/OverviewCoordinateReferenceSystems.pdf

##################################
# define sdc_raster object using population data and raster-object from before

# convert populatin data to SpatialPointsDataFrame
pop_data[,value:=1] # set value needed for sdcSpatial
sp_data <- sp::SpatialPointsDataFrame(coords = as.data.frame(pop_data[,.(x,y)]),
                                proj4string = sp::CRS(projection),
                                data = as.data.frame(pop_data))

# define sdc_raster object
pop_raster <- sdc_raster(x = sp_data, # SpatialPointsDataFrame holding the data
                         variable = "value", # variable cannot be missing
                         r = r, # set specific raster
                         min_count = 5 #  ~ cells should have at least min_count persons
                         )

# object summary 
pop_raster
# mean sensitivity score: share of sensitive cells

# plot object
plot(pop_raster,
     value="count", # what value should be plotted
     col=RColorBrewer::brewer.pal(8,"Blues")[-c(1:2)])



##################################
# check sensitive cells
sensitive_cells <- is_sensitive(pop_raster) # get sensitive cells
data_sens <- as.data.frame(sensitive_cells,xy=TRUE) # convert to data frame
setDT(data_sens)

# compare with input data (sanity check)
grid_data[,sensitive_direct:=count_original<5 & count_original>0]
grid_data <- merge(grid_data,data_sens,by=c("x","y"),all=TRUE)
# filter by all grid cells which are part of country
# use shape file for this
grid_data[,part_of_country:=get(geo) %in% area_shape$ID]
grid_data[,.N,by=.(sensitive,sensitive_direct,part_of_country)]



##################################
# apply protection method; check sensitive cells again and suppress if necessary

# 3 different methods implemented in sdcSpatial
# protect_quadtree() -> specifically for count data
# protect_neighborhood() -> sum over neighborhood in circle?
# protect_smooth() -> seems to be more applicable for numerical data

###
# apply quadtree -> example from ?protect_quadtree()
fined <- sdc_raster(enterprises, enterprises$fined, r=50)
plot(fined)
fined_qt <- protect_quadtree(fined)
plot(fined_qt)
###

# apply on own data
protected_raster <- protect_quadtree(pop_raster,
                                     max_zoom = 2 # number of zoom-steps - 1
                                     )
protected_raster # mean sensitivity for protected cells
plot(protected_raster,
     value="count", # what value should be plotted
     col=RColorBrewer::brewer.pal(8,"Blues")[-c(1:2)])

# get protected data in data.frame format for further analysis
protected_data <- as.data.frame(protected_raster$value,xy=TRUE)
setDT(protected_data)
grid_data[protected_data,count_quadtree:=count,on=.(x,y)]
grid_data

# apply new sensitivity flag 
# cells which are still sensitive will be supressed afterwards
sensitive_cells <- is_sensitive(protected_raster)
data_sens <- as.data.frame(sensitive_cells,xy=TRUE)
setDT(data_sens)
grid_data[data_sens,sensitive:=i.sensitive,on=.(x,y)]
grid_data[,.N,by=.(sensitive,part_of_country)]
grid_data[!is.na(sensitive) & part_of_country==FALSE]
# some cells have values which are not party of country

##################################
# get protected table and apply information loss function
grid_data[,suppress_value:=FALSE]
grid_data[sensitive==TRUE,suppress_value:=TRUE] # supress all cells which are still sensitive

grid_data[is.na(count_original),count_original:=0]
grid_data[is.na(count_quadtree),count_quadtree:=0]
infoLoss <- cellKey::ck_cnt_measures(grid_data[suppress_value==FALSE]$count_original,
                                     grid_data[suppress_value==FALSE]$count_quadtree)
infoLoss
# infoloss/noise is quite high in some cases
# appears if neighbouring cells have very different cell count and are averaged
# infoloss should probably be scaled with some distance

# share of supressed values
grid_data[part_of_country==TRUE & count_original>0,mean(suppress_value)]
# this is less than mean sensitivity score [0,1]:  0.1247862
# additional 500m raster which were not populated before
# some are not part of country

# get share of cells which have been aggregated
cells_aggregated <- grid_data[count_original!=count_quadtree,uniqueN(L000500)]
cells_total <- grid_data[,uniqueN(L000500)]
cells_aggregated/cells_total


##################################
# use higher level grid data and try to get an estimate on grid cell values

# higher grid level
geo2 <- paste0("L00",m*2)
count_geo <- paste0("count",geo2)
grid_data[part_of_country==TRUE,c(count_geo):=sum(count_original),by=c(geo2)]

# aggregated grid is the same as predefined higher level grid used by protect_quadtree
grid_data[part_of_country==TRUE,.(sum(count_original),sum(count_quadtree)),by=c(geo2)][V1==V2]

# example of raster plot
r <- raster(nrows=4, ncols=4)
r <- setValues(r, sample(1:ncell(r)))
plot(r,col=RColorBrewer::brewer.pal(9,"Blues"))

# check how close protected values are to original ones
# filter on cells which where used during protection and are now no longer sensitive
# all cells with values after the decimal point were surely used during protection
cells_protected <- grid_data[count_quadtree%%1!=0 & part_of_country==TRUE & sensitive == FALSE,unique(get(geo2))]
grid_data[get(geo2) %in% cells_protected,quantile(count_quadtree-count_original,probs=seq(0,1,.1))]

# check distance between original and protected value of previously sensitive cells
tab1 <- grid_data[get(geo2) %in% cells_protected & sensitive_direct==TRUE & sensitive == FALSE]
tab1 <- tab1[,.N,by=.(diff=abs(count_quadtree-count_original))] # group by abs difference
tab1 <- tab1[order(diff)]
tab1[,cumulative_p:=round(cumsum(N/sum(N))*100,2)]
tab1
tab1[diff<1] # % of protected cells which are less than 1 away from original
tab1[1:10]

# are there sensitive cells which are exactly the same after applying quadtree
cells_protected <- grid_data[sensitive_direct == TRUE,unique(get(geo2))]
tab1 <- grid_data[get(geo2) %in% cells_protected & sensitive_direct==TRUE & sensitive == FALSE]
tab1 <- tab1[,.N,by=.(diff=abs(count_quadtree-count_original))] # group by abs difference
tab1 <- tab1[order(diff)]
tab1[,cumulative_p:=round(cumsum(N/sum(N))*100,2)]
tab1
tab1[diff==0]  # % of sensitive cells which stay the same
tab1[1:10]
