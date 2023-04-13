# GADM administrative units
-----
This directory contains scripts to derive gridded land masks and 
country/region/district datasets from GADM shapes.

Unless you specify a predefined list of coordinates the scripts in this
directory will create a list of coordinates (the grid) which forms the basis
of all other input datasets.

## Input
GADM data needs to be downloaded from https://gadm.org/. 
  - Download the entire world as six separate layers (one for each level of 
    subdivision/aggregation). The scripts support either GeoPackage format or
    ESRI Shapefiles. 
  - File names: `gadm36_levels_gpkg.zip` or `gadm36_levels_shp.zip`.
  - Unzip in place.

If you use a different GADM version check file names in `gadm_load()` in
gadm_helper.R.

Optional: A CSV file providing the coordinates of a precribed grid.

This sample application is conducted for two resolutions, 5 arc-minute and 30
arc-minute, using GADM version 3.6. The grid at 30 arc-minutes uses a predefined
list of coordinates as specified in the included `gridlist_CRU.csv`.

## Software requirements
- `R`
- R packages `foreach`, `lwgeom`, `raster`, `rgdal`, `stringi`, `sf`, `units`,
  to use MPI backend for parallelization: `doMPI`, `Rmpi`, 
  or use `doParallel`, `parallel` to run on multiple CPUs of a local machine.


## Files included in this directory
  - 1_map_countries_to_grid_preparation.R: Script to prepare GADM (level 0
    and 1) data for grid intersection.
  - 2_map_countries_to_grid_intersection.R: Script to perform grid intersection
    between GADM data and grid.
  - 2_map_countries_to_grid_intersection_SLURM.sh: Example batch script for
    SLURM to run step 2 on the PIK cluster.
  - 3_map_countries_to_grid_collection.R: Script to collect data from grid
    intersection and create model input files.
  - 4_map_districts_to_grid_preparation.R: Script to prepare GADM (level 0-2)
    data for grid intersection.
  - 5_map_districts_to_grid_intersection.R: Script to perform grid intersection
    between GADM level 0-2 data and grid.
  - 5_map_districts_to_grid_intersection_SLURM.sh: Example batch script for
    SLURM to run step 5 on the PIK cluster.
  - 6_map_districts_to_grid_collection.R: Script to collect data from grid
    intersection and create files for further processing.
  - gadm_helper.R: Script defines a number of utility functions.
  - gadm_setup.R: Main setup script.
  - gridlist_CRU.csv: Example list of coordinates at 0.5° spatial resolution,
    corresponds to land grid used by the Climate Research Unit's time-series
    datasets of variations in climate with variations in other phenomena.
  - paper_analysis.R: Script to analyse generated input files for paper.
  - README.md: This file.
  - r_with_spatial_libs.sh: Bash script that loads all required modules to run
    R scripts on the PIK cluster.

## Set up for sample application
  - gadm_setup.R: 
    - Set `gadm_dir`, the base directory where scripts and unzipped GADM data
      are stored. Set to current working directory for sample application.
    - Set `lpj_res`, the spatial resolution of gridded data to be produced. Set
      to either 1/12 degree (5 arc-minutes) or 1/2 degree (30 arc-minutes)
    - Set `gadm_format` to either GeoPackage or ESRI Shapefile. Support may
      depend on your location installation. "ESRI Shapefile" used for sample
      application.
    - Setting `force_grid`: This setting allows to use a predefined grid
      (list of coordinates). Set to FALSE for 5 arc-minute dataset and set to
      TRUE for 30 arc-minute dataset.
    - Variable `griddata` set to empty data.frame for 5 arc-minute dataset,
      loaded from `gridlist_CRU.csv` for 30 arc-minute dataset.
    - All other settings left at default values in sample application.

By default the following files are created by step 3:
  - Grid: list of center coordinates of all cells that contain land area (in
    LPJmL input format).
  - Country code: largest country/region in each cell, list of two columns: 1)
    country index, 2) region index for some large countries defined in setting
    `include_regions` (in LPJmL input format).
  - Tables with meta information mapping indices used in country code file to
    country names and ISO codes (CSV tables).
  - Raster files of information contained in country code file.
  - Number of countries in each cell: Information used by land use data 
    processing scripts in this toolbox (LPJmL input format and raster file).
  - Land fraction: fraction of each cell covered by land polygons according to
    GADM (LPJmL input format and raster file).

Results of scripts for step 4-6 are required by land use data processing
scripts in this toolbox. By default the following files are created by step 6:
  - GADM level 0-2 code: largest country/region/district in each cell, list
    of three columns: 1) country index, 2) region index for all countries,
    not just countries defined in setting `include_regions`, 3) district index
    (LPJmL input format and raster file).
  - Table with meta information mapping indices used in GADM level 0-2 code
    file to names and ISO codes (CSV table).
