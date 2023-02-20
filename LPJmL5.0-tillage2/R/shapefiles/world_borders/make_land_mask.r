rm(list=ls());gc()

library(maptools)
library(raster)
library(rgeos)

wdir <- "./"


world_wgs <- readShapePoly(paste(wdir,"world_country_admin_boundary_shapefile_with_fips_codes.shp",sep=""),proj4string=CRS("+proj=longlat +datum=WGS84"))
world_wgs <- as(world_wgs,"SpatialPolygons")
mask_land <- as(extent(world_wgs),"SpatialPolygons")
proj4string(mask_land) <- proj4string(world_wgs)
mask_land <- gDifference(mask_land,world_wgs)

save(mask_land,file=paste0(wdir,"mask_land.RData"))
