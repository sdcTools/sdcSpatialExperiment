###### Prepare German pop. data for sdcSpatial #####

# ----- load libraries -----

library(sf)         # handling shape data
library(data.table) # fast import and data manipulation
library(dplyr)      # convenience
library(ggplot2)    # (optional) visualization


# ----- get required data -----

data_folder <- file.path(getwd(), "sdcSpatialExperiment_DE_ext") # customize if needed

web_pop   <- paste0("https://www.zensus2011.de/SharedDocs/Downloads/DE/Pressemitteilung/",
                    "DemografischeGrunddaten/csv_Bevoelkerung_100m_Gitter.zip",
                    "?__blob=publicationFile&v=2")
web_shp   <- paste0("https://daten.gdz.bkg.bund.de/produkte/vg/vg250_ebenen_0101/aktuell/",
                    "vg250_01-01.utm32s.shape.ebenen.zip")
web_100km <- paste0("https://daten.gdz.bkg.bund.de/produkte/sonstige/geogitter/aktuell/",
                    "DE_Grid_ETRS89-LAEA_100km.gpkg.zip")
web_100m  <- paste0("https://daten.gdz.bkg.bund.de/produkte/sonstige/geogitter/aktuell/",
                    "DE_Grid_ETRS89-LAEA_100m.csv.zip")

# download files (may need to raise timeout in R options for the larger files)
download.file(web_pop, paste0(data_folder, "/csv_Bevoelkerung_100m_Gitter.zip"),
              mode = "wb") # ~ 105 MB
download.file(web_shp, paste0(data_folder, "/vg250_01-01.utm32s.shape.ebenen.zip"),
              mode = "wb") # ~ 66 MB
download.file(web_100km, paste0(data_folder, "/DE_Grid_ETRS89-LAEA_100km.gpkg.zip"),
              mode = "wb") # ~ 1 MB
download.file(web_100m, paste0(data_folder, "/DE_Grid_ETRS89-LAEA_100m.zip"),
              mode = "wb") # ~ 291 MB

# unzip files
unzip(paste0(data_folder, "/csv_Bevoelkerung_100m_Gitter.zip"),    exdir = data_folder)
unzip(paste0(data_folder, "/vg250_01-01.utm32s.shape.ebenen.zip"), exdir = data_folder)
unzip(paste0(data_folder, "/DE_Grid_ETRS89-LAEA_100km.gpkg.zip"),  exdir = data_folder)
unzip(paste0(data_folder, "/DE_Grid_ETRS89-LAEA_100m.zip"),        exdir = data_folder)

path_grid_100m <- "DE_Grid_ETRS89-LAEA_100m/geogitter"
shp_grid_100km <- "DE_Grid_ETRS89-LAEA_100km.gpkg/geogitter/DE_Grid_ETRS89-LAEA_100km.gpkg"
shp_admi <- "vg250_01-01.utm32s.shape.ebenen/vg250_ebenen_0101/VG250_LAN.shp"
data_pop <- "Zensus_Bevoelkerung_100m-Gitter.csv"


# ----- select data -----

# The 100m x 100m grid is delivered in 100km x 100km chunks for the sake of file size.
# We start with the large grid in order to decide which chunks to load.

# read large grid
grid_100km <- st_read(file.path(data_folder, shp_grid_100km))

# In order to do the analysis for a single federal state, load state borders
admi_data <- st_read(file.path(data_folder, shp_admi)) %>%
  filter(GF == 4) %>%               # use only land area
  st_transform(st_crs(grid_100km))  # UTM to LAEA

# select a federal state to work with
admi_data$GEN
(fedstate <- admi_data$GEN[1])
fedstate_ags <- admi_data[admi_data$GEN == fedstate, ]$AGS # state identifier

# which 100km x 100km must be loaded for that state?
grid_load <- grid_100km[st_intersects(admi_data[admi_data$GEN == fedstate, ], 
                                           grid_100km)[[1]], ]$id

# (optional) visualize the selected chunks
ggplot() +
  geom_sf(data = grid_100km) +
  geom_sf(data = admi_data, fill = NA) +
  geom_sf(data = admi_data[admi_data$AGS == fedstate_ags, ], fill = "cyan") +
  geom_sf(data = grid_100km[grid_100km$id %in% grid_load, ], 
          fill = NA, color = "red")


# ----- read data -----

# read in selected chunks
grd <- vector("list", length = length(grid_load))
for(i in seq(grid_load)) {
  grd[[i]] <- fread(file.path(data_folder, path_grid_100m, 
                              paste0(grid_load[i], "_DE_Grid_ETRS89-LAEA_100m.csv")),
                    colClasses = c("character", rep("numeric", 10), "character"))
}

# combining selected chunks
grid_100m <- dplyr::bind_rows(grd)[, c(1:5, 12)]
names(grid_100m) <- c("cell_id", "x_sw", "y_sw", "x_mp", "y_mp", "ags")

# get rid of overlaps with other fed. states
grid_100m <- filter(grid_100m, substr(ags, 1, 2) == fedstate_ags)


# ----- add pop. numbers -----

# read in population data
pop <- fread(file.path(data_folder, data_pop)) %>%
  dplyr::select(Gitter_ID_100m, Einwohner) %>%
  dplyr::filter(Einwohner > 0)

# merge population to grid cells
pop_100m <- merge(grid_100m, pop, by.x = "cell_id", by.y = "Gitter_ID_100m")

sum(pop_100m$Einwohner) # inhabitants of selected state (w/o suppressed cells)

# expand to record-level data
pop_100m <- pop_100m[, .(1:Einwohner, x_mp, y_mp), by = .(cell_id)] %>%
  rename(person_id_cell = V1)
pop_100m$person_id <- 1:nrow(pop_100m)


# ----- save data -----

save(pop_100m, grid_100m, file = paste0(data_folder, "/pop_data_DE.RData"), compress = TRUE)

rm(list = ls())

