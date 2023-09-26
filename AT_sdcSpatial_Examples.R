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
# - define focus areas and calculate information loss
#


## -----------------------------------------------------------------------------------------
# load libraries
library(data.table)   # r package for data wrangling with large(r) data sets
library(raster)       # r package for defining raster-objects -> used by sdcSpatial
library(terra)        # used to read in shape file
library(sdcSpatial)   # r package to apply sdc methods on spatial data
library(viridis)      # r package for color shema (used in plots)

source("helpers.R")

## -----------------------------------------------------------------------------------------
# load population data
load("AT_test_data.RData")
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
gridsize <- 500 # dimension of rasters in meters
geo <-  paste(c("L",paste(rep(0,6-nchar(gridsize))),gridsize),collapse="")

# make table original cell counts
grid_data <- pop_data[,.(count_original=.N),by=c(geo)]
grid_data <- merge(grid_data, raster_cells, by=c(geo),all=TRUE)

# get center coordinates of grid cells - needed later
grid_data[,y:=as.numeric(paste0(substr(get(geo),6,10),"00"))]
grid_data[,x:=as.numeric(paste0(substr(get(geo),12,16),"00"))]
grid_data[is.na(count_original),count_original:=0]


## -----------------------------------------------------------------------------------------
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
area_shape <- terra::vect(file.path(unzipfolder, paste0(layer,".shp")))
# delete zip file
file.remove(my_shapefile)
# delete folder optional
# unlink(unzipfolder, recursive = TRUE)


# create raster object to be used by sdcSpatial
bounds <- ext(area_shape)
xmax <- bounds[2] + 3*gridsize # 3*gridsize -> include extra boundary on the map
xmin <- bounds[1] - 3*gridsize
ymax <- bounds[4] + 3*gridsize
ymin <- bounds[3] - 3*gridsize

ncol <- (xmax - xmin) / gridsize + 1
nrow <- (ymax - ymin) / gridsize + 1
projection <-  sf::st_crs(area_shape)$proj4string
r <- raster::raster(
  nrows = nrow,
  ncols = ncol,
  xmn = xmin - gridsize / 2,
  xmx = xmax + gridsize / 2,
  ymn = ymin - gridsize / 2,
  ymx = ymax + gridsize / 2,
  crs = sp::CRS(projection)
)
print(r)


# for overview on projections see:
# https://www.nceas.ucsb.edu/sites/default/files/2020-04/OverviewCoordinateReferenceSystems.pdf


## -----------------------------------------------------------------------------------------
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

pal <- rev(viridis(10)) 
# plot object
plot(pop_raster,
     value="count", # what value should be plotted
     col = pal)


##
# check sensitive cells
sensitive_cells <- is_sensitive(pop_raster) # get sensitive cells
data_sens <- as.data.frame(sensitive_cells,xy=TRUE) # convert to data frame
setDT(data_sens)

# compare with input data (sanity check)
grid_data[,sensitive_direct:=count_original<5 & count_original>0]
grid_data_check <- merge(grid_data,data_sens,by=c("x","y"),all=TRUE)
# filter by all grid cells which are part of country
# use shape file for this
grid_data_check[,part_of_country:=get(geo) %in% area_shape$ID]
grid_data_check[,.N,by=.(sensitive,sensitive_direct,part_of_country)]



## -----------------------------------------------------------------------------------------
# apply protection method; check sensitive cells again and suppress if necessary

# lets compare the following 3 different methods 
# protect_quadtree()
# protect_remove()
# protect_smooth()

# apply removal of sensitive cells
protected_rm <- remove_sensitive(x = pop_raster, min_count = 5)
print(protected_rm)

# apply quadtree
protected_qt <- protect_quadtree(x = pop_raster, min_count = 5, max_zoom = 2)
print(protected_qt)

# apply smoothing
protected_smooth_500 <- protect_smooth(x = pop_raster, bw = 500)
print(protected_smooth_500)


## -----------------------------------------------------------------------------------------
# define focus areas and calculate info loss for them

# set coordinates of center and radius for each focus area
focus_areas <- list(
  Area1 = list(center = c(x = 4790500, y = 2809500),
               radius = 26000),
  Area2 = list(center = c(x = 4547500, y = 2747000),
               radius = 15000),
  Area3 = list(center = c(x = 4298000, y = 2707500),
               radius = 14500) 
)

protected_grids <- list(remove=protected_rm,
                        quadtree = protected_qt,
                        smooth500 = protected_smooth_500)

# focus areas for unprotected maps
sdcObj_orig <- lapply(focus_areas,function(z,pop_ras){
  focus_sdcSpatialObj(pop_ras,focus=z)
},pop_ras=pop_raster)

# focus areas for protected maps
sdcObj_focus1 <- lapply(protected_grids,focus_sdcSpatialObj, focus=focus_areas$Area1)
sdcObj_focus2 <- lapply(protected_grids,focus_sdcSpatialObj, focus=focus_areas$Area2)
sdcObj_focus3 <- lapply(protected_grids,focus_sdcSpatialObj, focus=focus_areas$Area3)

sdcObj_focus1
sdcObj_focus2
sdcObj_focus3

plot(sdcObj_orig[[1]], value = "count", col = pal)
plot(sdcObj_focus1$remove, value = "count", col = pal)
plot(sdcObj_focus1$quadtree, value = "count", col = pal)
plot(sdcObj_focus1$smooth, value = "count", col = pal)

plot(sdcObj_orig[[2]], value = "count", col = pal)
plot(sdcObj_focus2$remove, value = "count", col = pal)
plot(sdcObj_focus2$quadtree, value = "count", col = pal)
plot(sdcObj_focus2$smooth, value = "count", col = pal)

plot(sdcObj_orig[[3]], value = "count", col = pal)
plot(sdcObj_focus3$remove, value = "count", col = pal)
plot(sdcObj_focus3$quadtree, value = "count", col = pal)
plot(sdcObj_focus3$smooth, value = "count", col = pal)


## -----------------------------------------------------------------------------------------
# Calculate HD and KWD for focus areas
sdcObj_focus <- list(Area1 = sdcObj_focus1,
                     Area2 = sdcObj_focus2,
                     Area3 = sdcObj_focus3)

# get raster from original input
original_raster <- pop_raster$value[["count"]]
v_o <- raster::getValues(original_raster)
v_o[is.na(v_o)] <- 0
xy <- raster::xyFromCell(original_raster, 1:raster::ncell(original_raster))
# transform coordinates to integer
xy_int <- copy(xy)
xy_int[,"x"] <- as.numeric(factor(xy_int[,"x"]))
xy_int[,"y"] <- as.numeric(factor(xy_int[,"y"]))

out <- list()
# iterate over focus areas
for(f in names(sdcObj_focus)){

  x_center <- focus_areas[[f]]$center[["x"]]
  y_center <- focus_areas[[f]]$center[["y"]]
  r_focus <- focus_areas[[f]]$radius
  
  # get position of center and radius for integer transformed coordinates
  int_center <- xy_int[xy[,"x"]==x_center & xy[,"y"]==y_center]
  int_radius <- r_focus/gridsize
  
  # get filter for focus area
  xy_filter <- xy[,"x"] %between% c(x_center-r_focus,x_center+r_focus) & 
    xy[,"y"] %between% c(y_center-r_focus,y_center+r_focus)
  
  # filter for focus area + 5*gridsize borders
  r_focus_ext <- r_focus+5*gridsize
  xy_filter_kwd <- xy[,"x"] %between% c(x_center-r_focus_ext,x_center+r_focus_ext) & 
    xy[,"y"] %between% c(y_center-r_focus_ext,y_center+r_focus_ext)
  
  # iterarte over protection methods
  for(m in names(sdcObj_focus[[f]])){
    
    protected_raster <- protected_grids[[m]]$value[["count"]]
    v_x <- raster::getValues(protected_raster)
    v_x[is.na(v_x)] <- 0
    # v_x[xy_filter==FALSE] <- 0 # needed for KWDSpatial::focusArea()
    
    sum(v_o[xy_filter_kwd==TRUE & xy_filter==TRUE])
    sum(v_x[xy_filter_kwd==TRUE & xy_filter==TRUE])
    
    # calculate HD
    v_o_hd <- v_o[xy_filter==TRUE]
    v_x_hd <- v_x[xy_filter==TRUE]
    hd <- (1/sqrt(2)) * sqrt((sum(sqrt(v_o_hd/sum(v_o_hd)) - sqrt(v_x_hd/sum(v_x_hd)))^2))
    
    
    # calculate KWD
    kwd <- SpatialKWD::focusArea(Coordinates = xy_int[xy_filter_kwd,],
                                 Weights = cbind(v_o[xy_filter_kwd],v_x[xy_filter_kwd]),
                                 x = int_center[1],
                                 y = int_center[2],
                                 radius = int_radius,
                                 method = "exact", area = "Inf")
    
    out_i <- data.table(Area=f,Method=m,hd=hd,kwd=kwd$distance/sum(v_x_hd))
    out <- c(out,list(out_i))
  }
}

out <- rbindlist(out)
out

# for 0 distance for KWD see current issue also raised in Example_spatialKWD_focusAreas.Rmd





