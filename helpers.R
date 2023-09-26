# variouse help functions

# get "focused" sdcSpatialObj 
focus_sdcSpatialObj <- function(pop_ras,focus,gridsize=500){#
  
  # find cell index of centroid in larger raster
  cntr <- c(raster::colFromX(pop_ras$value, x = focus$center["x"]),
            raster::rowFromY(pop_ras$value, y = focus$center["y"]))
  
  wcells <- focus$radius/gridsize
  
  # create extent objects for chosen positions & sizes
  extt <- extent(pop_ras$value, cntr[2] - wcells , 
                 cntr[2] + wcells ,
                 cntr[1] - wcells ,
                 cntr[1] + wcells)
  
  pop_focus <- copy(pop_ras)
  pop_focus$value <- crop(pop_focus$value,extt)
  
  return(pop_focus)
}