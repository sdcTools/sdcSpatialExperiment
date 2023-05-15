# Scripts for testing R-package [sdcSpatial](https://github.com/edwindj/sdcSpatial) on population grid cells

Scripts to apply R-Package [sdcSpatial](https://github.com/edwindj/sdcSpatial) on population grids. First install need packages:

```
install.packages(c("data.table","sf","raster","sdcSpatial"))
```

Then run script `generate_data.R` which generates a random population based on the cell grids from Austria (see https://data.statistik.gv.at/web/meta.jsp?dataset=OGDEXT_RASTER_1).
This data set is the input for the script `sdcSpatial_template.R` which goes through the following steps:

- Load needed libraries
- Load dummy population data
- Define "raster"-object for grid cells to be used with sdcSpatial
- Define sdc_raster object using population data and raster-object from before
- Check sensitive cells
- Apply protection method; check sensitive cells again and suppress if necessary
- Get protected table and apply information loss function
- Check if the original value of protected cells can be re-estimated


