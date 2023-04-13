Official download from HydroSheds website (version 1.1) is in Tiff format:
hyd_glo_dir_5m.zip

Converted into HydroSheds_v1.1.asc in R by:
```
library(raster)
hyd <- raster("hyd_glo_dir_5m.tif")
# DDM does not cover all land cells; need to expand to global extent and fill with missing value
global_extent <- extent(-180, 180, -90, 90)
hyd_expanded <- extend(hyd, global_extent)
# Write expanded raster as ASCII grid
writeRaster(hyd_expanded, filename = "HydroSheds_v1.1.asc", datatype = "INT2S", NAflag = -9)
```
