#################################################
# script to generate dummy data for Austrian population
#


## ---------------------------------------------------------------
# load libraries
library(data.table)
library(sf)


## ---------------------------------------------------------------
# download data from open portal
url <- "https://data.statistik.gv.at/data/OGDEXT_RASTER_1_STATISTIK_AUSTRIA_LAEA_GPKG.zip"

# download .zip
zipfile <- file.path(getwd(),basename(url))
download.file(url,zipfile)

# unzip file
unzipfolder <- paste0(sub(".zip","",zipfile))
unzip(zipfile,exdir=unzipfolder)
# delete zip file
file.remove(zipfile)

# read geopackage file
st_layers(file.path(unzipfolder,"STATISTIK_AUSTRIA_RASTER_LAEA_3035.gpkg"))

l500 <- st_read(file.path(unzipfolder,"STATISTIK_AUSTRIA_RASTER_LAEA_3035.gpkg"), layer = "l000500")

# delete folder
unlink(unzipfolder, recursive = TRUE)


## ---------------------------------------------------------------
# define lookup for 500m, 1km and 2km gridcells
raster_cells <- data.table(L000500=l500$cellcode)
raster_cells[,y:=as.numeric(paste0(substr(L000500,6,10),"00"))]
raster_cells[,x:=as.numeric(paste0(substr(L000500,12,16),"00"))]
raster_cells[,c("y_1k","x_1k"):=.(floor(y/1000),floor(x/1000))]
raster_cells[,L001000:=paste0("1kmN",y_1k,"E",x_1k)]
raster_cells[,c("y_2k","x_2k"):=.(y_1k-y_1k%%2,x_1k-x_1k%%2)]
raster_cells[,L002000:=paste0("2kmN",y_2k,"E",x_2k)]
raster_cells[,c("x_1k","y_1k","x_2k","y_2k","x","y"):=NULL]
setorder(raster_cells,L002000,L001000,L000500)


## ---------------------------------------------------------------
# randomly populate raster_cells to build dummy population
set.seed(202305)
share_unpop_cells <- runif(1,min=0.65,0.75)

pop_data <- copy(raster_cells)
# populate some grid cells with 0
pop_data[sample(.N,round(share_unpop_cells*.N)),pop_count:=0]
# remaining gridcells populate with some values between 1 and 10000
pop_data[is.na(pop_count),pop_count:=sample(1:10000,.N,replace=TRUE,prob=1/(abs(1:10000 - 10)+1)^(1.455))]
pop_data[,sum(pop_count)]

pop_data <- pop_data[pop_count>0,.(1:pop_count),by=.(L000500,L001000,L002000)]
pop_data[,V1:=NULL]
pop_data[,PID:=1:.N] # dummy PID
# coordinates of each person = 500m grid cell coordinates
pop_data[,y:=as.numeric(paste0(substr(L000500,6,10),"00"))]
pop_data[,x:=as.numeric(paste0(substr(L000500,12,16),"00"))]


## ---------------------------------------------------------------
# save data
save(pop_data,raster_cells,file="AT_test_data.RData",compress=TRUE)

