# Scripts for testing R-package [sdcSpatial](https://github.com/edwindj/sdcSpatial) on population grid cells

Scripts to apply R-Package [sdcSpatial](https://github.com/edwindj/sdcSpatial) on population grids. First install needed packages:

```
install.packages(c("data.table","terra","sf","raster","sdcSpatial", "dplyr", "tidyr", "ggplot2", "viridis", "spatstat", "SpatialKWD"))
```

## Example Austrian dummy data

Run script `AT_generate_data.R` which generates a random population based on the cell grids from Austria (see https://data.statistik.gv.at/web/meta.jsp?dataset=OGDEXT_RASTER_1).
This data set is the input for the script `AT_sdcSpatial_Examples.R` which goes through the following steps:

- Load needed libraries
- Load dummy population data
- Define "raster"-object for grid cells to be used with sdcSpatial
- Define sdc_raster object using population data and raster-object from before
- Check sensitive cells
- Apply protection method; check sensitive cells again and suppress if necessary
- Get protected table and apply information loss function, look at HD and KWD

## Example German data

2 scripts are supplied:

`DE_read_data.R` can be run to recreate pop_data_DE.RData
`DE_sdcSpatial_Examples.R` runs sdcSptial experiments on pop_data_DE.RData


`DE_read_data.R`

- creates geocoded microdata from published 2011 census results for test purposes
- test area is one of the German federal states (can be chosen in the script)
- outputs: pop_data_DE.RData (stored default: Schleswig-Holstein)
- needs the following data:
  a) German population aggregates on 100m grid
source: 
https://www.zensus2011.de/SharedDocs/Downloads/DE/Pressemitteilung/DemografischeGrunddaten/csv_Bevoelkerung_100m_Gitter.zip
  b) 100m grid definitions for Germany
source:
https://gdz.bkg.bund.de/index.php/default/open-data/geographische-gitter-fur-deutschland-in-lambert-projektion-geogitter-inspire.html
  c) 100km grid definitions for Germany (to select chunks of the 100m grid)
source:
https://gdz.bkg.bund.de/index.php/default/open-data/geographische-gitter-fur-deutschland-in-lambert-projektion-geogitter-inspire.html
  d) shapefile for German federal states (in UTM32)
https://daten.gdz.bkg.bund.de/produkte/vg/vg250_ebenen_0101/aktuell/


`DE_sdcSpatial_Examples.R`

- run the actual sdcSpatial experiments and utility analyses
- inputs: pop_data_DE.RData

