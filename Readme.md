# Scripts for testing R-package [sdcSpatial](https://github.com/edwindj/sdcSpatial) on German grid data

Contact:
Martin Möhler
martin.moehler@destatis.de


2 scripts are supplied:

`01_read_data_DE.R` can be run to recreate pop_data_DE.RData; it can also
be skipped and the pre-computed file used to run `02_sdcSpatial_DE.R`.

Packages needed:

```
install.packages(c("data.table","sf","raster","sdcSpatial", "dplyr", "tidyr", "ggplot2", "viridis", "spatstat", "SpatialKWD"))
```

`01_read_data_DE.R`

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


`02_sdcSpatial_DE.R`

- run the actual sdcSpatial experiments and utility analyses
- inputs: pop_data_DE.RData

