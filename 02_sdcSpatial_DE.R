##### Apply sdcSpatial on 100m grid counts tables [DE] #####
# 
# Script does the following:
# (1) load required tools + pre-computed pop. and grid data [see: 01_read_data_DE.R]
# (2) create + inspect sdc_raster object from pop. data
# (3) apply protection algorithms (removal, quadtree, kernel smoothing)
# (4) inspect protected raster maps, check & compare various utility measures
# (5) illustrate two aspects of KWD as info loss measure


# ----- (1a) load libraries -----

library(raster)       # grid data handling
library(sdcSpatial)   # map protection algorithms
library(spatstat)     # geospatial functions
library(SpatialKWD)   # fast KWD approximation
library(ggplot2)      # plotting
library(viridis)      # plotting
library(dplyr)        # convenience
library(tidyr)        # convenience


# ----- (1b) load data  -----

# customize to local setting
path_data <- "sdcSpatialExperiment_DE_ext"

load(file.path(path_data, "pop_data_DE.RData"))

head(pop_100m)
## data.table is one row per person ############################## ###
# cell_id:        INSPIRE name of grid cell
# person_id_cell: unique ID for person within a cell
# x_mp:           x-coordinate of person (cell centroid)
# y_mp:           y-coordinate of person (cell centroid)
# person_id:      unique person ID
# [all coordinates in ETRS89-LAEA]

head(grid_100m)
## data.table is one row per grid cell (incl. uninhabited) ####### ###
# cell_id: INSPIRE name of grid cell
# x_sw:    x-coordinate of southwest-corner
# y_sw:    y-coordinate of southwest-corner
# x_mp:    x-coordinate of centroid
# y_mp:    y-coordinate of centroid
# ags:     municipality identifier
# [all coordinates in ETRS89-LAEA]


# ----- (2a) set up sdc_raster  -----

# prepare input raster (incl. uninhabited)
r <- raster(xmn = min(grid_100m$x_sw), xmx = max(grid_100m$x_sw + 100),
            ymn = min(grid_100m$y_sw), ymx = max(grid_100m$y_sw + 100),
            res = 100,
            crs = CRS("+init=epsg:3035"))

# to sdc_raster
pop100m_raster <- sdc_raster(as.matrix(pop_100m[, c("x_mp", "y_mp")]),
                             variable  = 1,
                             r         = r,
                             min_count = 10,  # require 10-anonymity
                             max_risk  = 1.0) # no magnitude / freq. criterion

pal <- rev(viridis(10)) # color palette for plotting
plot(pop100m_raster, value = "count", col = pal)


# ----- (2b) inspect sdc_raster -----

# subset sensitive cells
raster_sens <- is_sensitive(pop100m_raster)

# no. of sensitive cells
(ncells_sens <- sum(raster_sens@data@values, na.rm = TRUE))
# share of sensitive cells (in all inhabited)
sensitivity_score(pop100m_raster)
# verify:
ncells_sens / length(unique(pop_100m$cell_id))

# How many individuals are at risk?
sens_vals <- pop100m_raster$value$count[raster_sens]
sum(sens_vals)                  # no. of indiv. at risk
sum(sens_vals) / nrow(pop_100m) # share of indiv. at risk

# How are cell counts < k distributed?
ggplot(data = data.frame(sens_vals), aes(sens_vals)) +
  geom_histogram(binwidth = 1, color = "black", fill = "lightgrey") +
  scale_x_continuous(n.breaks = 7) +
  xlab("no. inhabitants in sensitive cells")


# ----- (3) apply map protection algorithms -----

# simply remove sensitive cells
pop100m_raster_rm <- remove_sensitive(pop100m_raster)

# quadtree method
pop100m_raster_qt1 <- protect_quadtree(pop100m_raster, max_zoom = 2)
pop100m_raster_qt2 <- protect_quadtree(pop100m_raster, max_zoom = 3)

# kernel smoothing (~ 12min)
pop100m_raster_sm <- protect_smooth(pop100m_raster, bw = 100)


sensitivity_score(pop100m_raster_rm)
sensitivity_score(pop100m_raster_qt1)
sensitivity_score(pop100m_raster_qt2)
sensitivity_score(pop100m_raster_sm) 

# remove remaining sensitive cells (if any)
pop100m_raster_qt1 <- remove_sensitive(pop100m_raster_qt1)
pop100m_raster_qt2 <- remove_sensitive(pop100m_raster_qt2)
pop100m_raster_sm  <- remove_sensitive(pop100m_raster_sm)


# ----- (4a) inspect protected rasters -----

# visual comparison
plot(pop100m_raster_rm,  "count", col = pal, sensitive = FALSE)
plot(pop100m_raster_qt1, "count", col = pal, sensitive = FALSE)
plot(pop100m_raster_qt2, "count", col = pal, sensitive = FALSE)
plot(pop100m_raster_sm,  "count", col = pal, sensitive = FALSE)


# What's the cell-wise abs. difference between protected vs. unprotected?

difference_rasters <- function(x, y) {
  
  # function is only needed when rasters contain NAs
  na <- is.na(x) & is.na(y)
  x[is.na(x)] <- 0
  y[is.na(y)] <- 0
  z <- x - y
  z[na] <- NA
  z
}

# calculate differences
pop100m_raster_rm$value$diff  <- difference_rasters(pop100m_raster_rm$value$count, 
                                                    pop100m_raster$value$count)
pop100m_raster_qt1$value$diff <- difference_rasters(pop100m_raster_qt1$value$count, 
                                                    pop100m_raster$value$count)
pop100m_raster_qt2$value$diff <- difference_rasters(pop100m_raster_qt2$value$count,
                                                    pop100m_raster$value$count)
pop100m_raster_sm$value$diff  <- difference_rasters(pop100m_raster_sm$value$count,
                                                    pop100m_raster$value$count)

# plot differences
plot(pop100m_raster_rm$value$diff,  col = pal, main = "removal")
plot(pop100m_raster_qt1$value$diff, col = pal, main = "quadtree zoom 2")
plot(pop100m_raster_qt2$value$diff, col = pal, main = "quadtree zoom 3")
plot(pop100m_raster_sm$value$diff,  col = pal, main = "kernel smoothing")

# How much have the values of sensitive cells been changed?
hist(pop100m_raster_rm$value$diff[raster_sens],  main = "removal")
hist(pop100m_raster_qt1$value$diff[raster_sens], main = "quadtree zoom 2")
hist(pop100m_raster_qt2$value$diff[raster_sens], main = "quadtree zoom 3")
hist(pop100m_raster_sm$value$diff[raster_sens],  main = "kernel smoothing")

# How many sensitive cells are unchanged?
sum(pop100m_raster_rm$value$diff[raster_sens]  == 0)
sum(pop100m_raster_qt1$value$diff[raster_sens] == 0)
sum(pop100m_raster_qt2$value$diff[raster_sens] == 0)
sum(pop100m_raster_sm$value$diff[raster_sens]  == 0)
# What share of sensitive cells is unchanged?
sum(pop100m_raster_rm$value$diff[raster_sens]  == 0) / ncells_sens
sum(pop100m_raster_qt1$value$diff[raster_sens] == 0) / ncells_sens
sum(pop100m_raster_qt2$value$diff[raster_sens] == 0) / ncells_sens
sum(pop100m_raster_sm$value$diff[raster_sens]  == 0) / ncells_sens


# Which cells were NA in orig. and now aren't?

pop100m_raster_rm$value$outside <- 
  is.na(pop100m_raster$value$count) & !is.na(pop100m_raster_rm$value$count)
pop100m_raster_qt1$value$outside <- 
  is.na(pop100m_raster$value$count) & !is.na(pop100m_raster_qt1$value$count)
pop100m_raster_qt2$value$outside <- 
  is.na(pop100m_raster$value$count) & !is.na(pop100m_raster_qt2$value$count)
pop100m_raster_sm$value$outside <- 
  is.na(pop100m_raster$value$count) & !is.na(pop100m_raster_sm$value$count)

sum(getValues(pop100m_raster_rm$value$outside))  # removal
sum(getValues(pop100m_raster_qt1$value$outside)) # quadtree zoom 2
sum(getValues(pop100m_raster_qt2$value$outside)) # quadtree zoom 3
sum(getValues(pop100m_raster_sm$value$outside))  # kernel smoothing


# ----- (4b) utility measures - general -----

###### fun: calculate utility values for protected maps ##### ###
# x:       sdc_raster object - after map protection
# orig:    sdc_raster object - original before protection
# value:   value to compare ("count", "sum", ...)
# measure: one or several of the following
#           - bin-by-bin: -
#            RMSE: Root Mean Squared Error
#            HD:   Hellinger's Distance
#           - cross-bin: -
#            Moran: change in (global) Moran's I
#            KWD: Kantorovic-Wasserstein Distance
# ...  :   additional parameters for SpatialKWD::compareOneToOne()

get_utility <- function(x, orig, 
                        value = "count", 
                        measure = c("RMSE", "HD", "Moran", "KWD"), 
                        ...) {
  
  u <- vector("numeric", length = length(measure))
  names(u) <- measure
  
  # extract value raster of interest
  r_x <- x$value[[value]]
  r_o <- orig$value[[value]]
  
  # extract cell values
  v_x <- raster::getValues(r_x)
  v_o <- raster::getValues(r_o)
  
  # set NAs to 0
  v_x[is.na(v_x)] <- 0
  v_o[is.na(v_o)] <- 0
  
  for (i in seq(u)) {
    
    if (measure[i] == "RMSE") { 
      
      # compute root mean squared error
      u[i]  <- sqrt(mean((v_x - v_o)^2))
    }
    
    if (measure[i] == "HD") {
      
      # compute Hellinger's distance
      u[i] <- sqrt(sum((sqrt(v_o) - sqrt(v_x))^2)) * (1/(sqrt(2)))
    }
    
    if (measure[i] == "Moran") {
      
      # change in (global) Moran's I
      u[i] <- Moran(r_x) - Moran(r_o)
    }
    
    if (measure[i] == "KWD") {
      
      # rescale to balance mass
      v_x_sc <- v_x * (sum(v_o) / sum(v_x))
      
      # approximate Kantorovic-Wasserstein distance
      xy <- raster::xyFromCell(r_x, 1:raster::ncell(r_x))
      u[i] <- SpatialKWD::compareOneToOne(Coordinates = xy,
                                          Weights = cbind(v_o, v_x_sc),
                                          ...)$distance
    }
  }
  
  u
} 


# assess protected maps with utility measures

u_measures <- c("RMSE", "HD")

u_rm  <- get_utility(pop100m_raster_rm, 
                     orig = pop100m_raster, value = "count",
                     measure = u_measures)
u_qt1 <- get_utility(pop100m_raster_qt1, 
                     orig = pop100m_raster, value = "count",
                     measure = u_measures)
u_qt2 <- get_utility(pop100m_raster_qt2, 
                     orig = pop100m_raster, value = "count",
                     measure = u_measures)
u_sm  <- get_utility(pop100m_raster_sm,
                     orig = pop100m_raster, value = "count",
                     measure = u_measures)

results <- as.data.frame(rbind(u_rm, u_qt1, u_qt2, u_sm))
rownames(results) <- c("removal", "qt_zoom2", "qt_zoom3", "smoothing")
results


# ----- (4c) utility measures - sub-area (incl. KWD) -----

# Even with the fast approximation algorithm, computing KWD for a large map with
# small cell size can take a long time!
# Therefore we only demonstrate on a smaller sub-area

# select a rectangular sub-area as bounding box of regional ID
ags_sub <- "01002000" # city of Kiel
grid_sub <- grid_100m[grid_100m$ags == ags_sub, ]
ext <- extent(min(grid_sub$x_sw), max(grid_sub$x_sw) + 100,
              min(grid_sub$y_sw), max(grid_sub$y_sw) + 100)

# crop rasters to sub-area
sub_raster     <- pop100m_raster
sub_raster_rm  <- pop100m_raster_rm
sub_raster_qt1 <- pop100m_raster_qt1
sub_raster_qt2 <- pop100m_raster_qt2
sub_raster_sm  <- pop100m_raster_sm
sub_raster$value     <- crop(sub_raster$value,     ext)
sub_raster_rm$value  <- crop(sub_raster_rm$value,  ext)
sub_raster_qt1$value <- crop(sub_raster_qt1$value, ext)
sub_raster_qt2$value <- crop(sub_raster_qt2$value, ext)
sub_raster_sm$value  <- crop(sub_raster_sm$value,  ext)

# visually inspect sub-area
par(mfrow = c(2, 3))
plot(sub_raster,     "count", col = pal, sensitive = FALSE, main = "unprotected")
plot(sub_raster_rm,  "count", col = pal, sensitive = FALSE, main = "removal")
plot(sub_raster_qt1, "count", col = pal, sensitive = FALSE, main = "quadtree zoom 2")
plot(sub_raster_qt2, "count", col = pal, sensitive = FALSE, main = "quadtree zoom 3")
plot(sub_raster_sm,  "count", col = pal, sensitive = FALSE, main = "kernel smoothing (100m bw)")
par(mfrow = c(1, 1))

# assess sub-area with full set of utility measures

u_sub_rm  <- get_utility(sub_raster_rm,  orig = sub_raster, value = "count")
u_sub_qt1 <- get_utility(sub_raster_qt1, orig = sub_raster, value = "count")
u_sub_qt2 <- get_utility(sub_raster_qt2, orig = sub_raster, value = "count")
u_sub_sm  <- get_utility(sub_raster_sm,  orig = sub_raster, value = "count")

results_sub <- as.data.frame(rbind(u_sub_rm, u_sub_qt1, u_sub_qt2, u_sub_sm))
rownames(results_sub) <- c("removal", "qt_zoom2", "qt_zoom3", "smoothing")
results_sub

# inspect results
ggplot(results_sub %>% mutate(method = rownames(results_sub))) +
  geom_point(aes(HD, KWD, color = method))


# ----- (5a) KWD illustration 1 -----

# KWD can yield a different rank-order of protection mechanisms than bin-by-bin
# utility metrics, if the mechanism shifts distribution mass without concern for
# geographic neighborhood (e.g. removal, random swapping).

# If, however, the protection mechanism respects geographic neighborhood by design
# (quadtree, smoothing, small random perturbation), then bin-by-bin metrics and
# KWD give the same answer (same rank order w.r.t. utility).

data("dwellings") # pseudo household locations (from sdcSpatial)
dw_or <- sdc_raster(dwellings[, c("x", "y")], variable = 1)

# do both methods with increasingly strong parameters
dw_qt1 <- protect_quadtree(dw_or, max_zoom = 2)
dw_qt2 <- protect_quadtree(dw_or, max_zoom = 3)
dw_qt3 <- protect_quadtree(dw_or, max_zoom = 4)
dw_sm1 <- protect_smooth(dw_or, bw = 200)
dw_sm2 <- protect_smooth(dw_or, bw = 300)
dw_sm3 <- protect_smooth(dw_or, bw = 400)
  
# get utility
u_qt1 <- round(get_utility(dw_qt1, dw_or, value = "count"), 3)
u_qt2 <- round(get_utility(dw_qt2, dw_or, value = "count"), 3)
u_qt3 <- round(get_utility(dw_qt3, dw_or, value = "count"), 3)
u_sm1 <- round(get_utility(dw_sm1, dw_or, value = "count"), 3)
u_sm2 <- round(get_utility(dw_sm2, dw_or, value = "count"), 3)
u_sm3 <- round(get_utility(dw_sm3, dw_or, value = "count"), 3)

(results_dw <- as.data.frame(rbind(u_qt1, u_qt2, u_qt3, u_sm1, u_sm2, u_sm3)))

# plot results
ggplot(results_dw, aes(RMSE, KWD)) +
  geom_point() +
  geom_text(aes(label = rownames(results_dw), x = RMSE + 0.5, y = KWD + 0.03)) +
  geom_abline(intercept = coef(lm(KWD ~ RMSE, data = results_dw))[1],
              slope     = coef(lm(KWD ~ RMSE, data = results_dw))[2],
              color = "blue", lty = "dashed")

# inspect visual impressions
par(mfrow = c(2, 3))
plot(dw_qt1, "count", sensitive = FALSE, col = pal,
     main = paste("KWD", u_qt1[4], " RMSE", u_qt1[1]))
plot(dw_qt2, "count", sensitive = FALSE, col = pal,
     main = paste("KWD", u_qt2[4], " RMSE", u_qt2[1]))
plot(dw_qt3, "count", sensitive = FALSE, col = pal,
     main = paste("KWD", u_qt3[4], " RMSE", u_qt3[1]))
plot(dw_sm1, "count", sensitive = FALSE, col = pal,
     main = paste("KWD", u_sm1[4], " RMSE", u_sm1[1]))
plot(dw_sm2, "count", sensitive = FALSE, col = pal,
     main = paste("KWD", u_sm2[4], " RMSE", u_sm2[1]))
plot(dw_sm3, "count", sensitive = FALSE, col = pal,
     main = paste("KWD", u_sm3[4], " RMSE", u_sm3[1]))
par(mfrow = c(1, 1))


# ----- (5b) KWD illustration 2 -----

# KWD seems to give preference to the quadtree method over smoothing;
# this does not necessarily correspond to human visual impression

# make small toy data set
set.seed(3000)
pts_exa <- as.data.frame(spatstat.random::rMatClust(20, 0.1, 100))
pts_exa <- pts_exa * 1000

test_raster <- sdc_raster(pts_exa, variable = 1, r = 50,
                          min_count = 5, max_risk = 1.0)

## compare: quadtree vs. smoothing

test_raster_qt <- protect_quadtree(test_raster)
test_raster_sm <- protect_smooth(test_raster, bw = 50)

u_test_qt <- round(get_utility(test_raster_qt, orig = test_raster,
                               value = "count", measure = c("RMSE", "KWD")), 4)
u_test_sm <- round(get_utility(test_raster_sm, orig = test_raster,
                               value = "count", measure = c("RMSE", "KWD")), 4)

par(mfrow = c(1, 3))
plot(test_raster, "count", sensitive = FALSE,
     main = "original")
plot(test_raster_qt, "count", sensitive = FALSE,
     main = "quadtree", 
     sub = paste("KWD:", u_test_qt[2], "  RMSE:", u_test_qt[1]))
plot(test_raster_sm, "count", sensitive = FALSE,
     main = "smoothing", 
     sub = paste("KWD:", u_test_sm[2], "  RMSE:", u_test_sm[1]))
par(mfrow = c(1, 1))

# The KWD is, in principle, well suited to compare 2 raster maps, but does it 
# capture 'similarity of visual impression'? Some doubt is in order.
# KWD seems to penalize small distribution mass, transported over medium distances, 
# more than the eyeball.

# -> Idea: remove small distribution mass; compare only 'hot spots'


# remove the lower [tsh]% of distribution mass
tsh <- 0.2
hs    <- tsh * max(getValues(test_raster$value$count),    na.rm = TRUE)
hs_qt <- tsh * max(getValues(test_raster_qt$value$count), na.rm = TRUE)
hs_sm <- tsh * max(getValues(test_raster_sm$value$count), na.rm = TRUE)

# raster layer of 'hot spots' - original
hsr <- test_raster$value$count
hsr@data@values <- ifelse(hsr@data@values > hs, hsr@data@values, NA)
test_raster$value$hs <- hsr
# raster layer of 'hot spots' - quadtree
hsr <- test_raster_qt$value$count
hsr@data@values <- ifelse(hsr@data@values > hs_qt, hsr@data@values, NA)
test_raster_qt$value$hs <- hsr
# raster layer of 'hot spots' - smoothing
hsr <- test_raster_sm$value$count
hsr@data@values <- ifelse(hsr@data@values > hs_sm, hsr@data@values, NA)
test_raster_sm$value$hs <- hsr

# re-compute KWD & RMSE for 'hot spots' only
u_test_qt_hs <- round(get_utility(test_raster_qt, orig = test_raster,
                                  value = "hs", measure = c("RMSE", "KWD")), 4)
u_test_sm_hs <- round(get_utility(test_raster_sm, orig = test_raster,
                                  value = "hs", measure = c("RMSE", "KWD")), 4)

par(mfrow = c(1, 3))
plot(test_raster, "hs", sensitive = FALSE,
     main = "original")
plot(test_raster_qt, "hs", sensitive = FALSE,
     main = "quadtree", 
     sub = paste("KWD:", u_test_qt_hs[2], "  RMSE:", u_test_qt[1]))
plot(test_raster_sm, "hs", sensitive = FALSE,
     main = "smoothing", 
     sub = paste("KWD:", u_test_sm_hs[2], "  RMSE:", u_test_sm[1]))
par(mfrow = c(1, 1))

# KWD for 'hot spot'-maps, rather than full maps, seems to be more in line with 
# 'visual experience'. How to finally assess the map?
# Maybe some (weighted) average of KWD for full map vs. for 'hot spot' map?

