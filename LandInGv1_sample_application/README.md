# LandInG: Land Input Generator
-----
Main contact: Sebastian Ostberg, ostberg@pik-potsdam.de

Affiliation: Potsdam Institute for Climate Impact Research (PIK), Potsdam,
Germany

This project contains a collection of scripts to derive basic input datasets for
terrestrial ecosystem models from diverse and partially conflicting data
sources. It has been developed with a focus on the open-source dynamic global
vegetation, hydrology and crop growth model LPJmL (Lund-Potsdam-Jena with
managed Land) and complements the
[LPJmL model code release](https://github.com/PIK-LPJmL/LPJmL). This toolbox
does *not* cover climate inputs.

This specific code collection covers a sample application to produce sets of
LPJmL inputs at 5 arc-minute and 30 arc-minute spatial resolution, as presented
in the paper: Ostberg, S., Müller, C., Heinke, J., and Schaphoff, S.:
LandInG 1.0: A toolbox to derive input datasets for terrestrial ecosystem
modelling at variable resolutions from heterogeneous sources, submitted to
Geosci. Model Dev.

The generic version of this toolbox is distributed through
https://github.com/PIK-LPJmL/LandInG

All source code, configuration and parameter files are subject to 
Copyright (C) by the Potsdam Institute for Climate Impact Research (PIK), see
the file COPYRIGHT, and are licensed under the GNU AFFERO GENERAL PUBLIC LICENSE
Version 3, see the file LICENSE.

## Structure:
  - elevation_paper: subdirectory containing scripts to derive an elevation input
  - fertilizer_paper: subdirectory containing scripts to derive fertilizer and
    manure inputs from various sources
  - gadm_paper: subdirectory containing scripts to process GADM data (grid,
    country/region/district codes, land fraction in each cell etc.)
  - lakes_rivers_paper: subdirectory containing scripts to process GLWD data and
    extract cell fractions covered by lakes and rivers
  - landuse_paper: subdirectory containing scripts to generate a landuse dataset
    from various sources
  - river_routing_paper: subdirectory containing scripts to derive inputs
    related to river routing
  - reservoirs_paper: subdirectory containing scripts to derive dam/reservoir
    input
  - soil_paper: subdirectory containing scripts to derive soil inputs
  - COPYRIGHT: copyright information
  - LICENSE: License information
  - README.md: this file
  - lpjml_format_helper_functions.R: script definining a number of utility
    functions to work with the LPJmL input format

## Steps used for sample application
  1) First run the scripts in the `gadm_paper` directory for 5 arc-minutes
    (no forced grid) and 30 arc-minutes (forced grid), which generate a grid
    file that is used by all the other input generating scripts. Also generate
    several administrative masks for both spatial resolutions.
  2) Run scripts in `soil_paper` and `lakes_rivers_paper` to create a soil input
    and a lake input based on the grid files created in step 1. Lake inputs
    created twice per spatial resolution, once based on GLWD raster data and
    second based on GLWD polygon data.
  3) Run scripts in `river_routing_paper` to create river routing for two
    different drainage direction maps per spatial resolution (4 in total).
    Create neighbour irrigation inputs using 4 different neighbour assigment
    options for each of the 4 drainage direction maps (16 in total).
  4) Run scripts in `elevation_paper` to create elevation inputs for each grid
    created in step 1. Run scripts in `reservoirs_paper` to generate reservoir
    inputs: 4 different assigment strategies per drainage direction map
    processed in step 3 (16 in total).
  5) Run scripts in `landuse_paper` to generate land use inputs. Data processing
    from various source datasets to gridded time series of crop-specific rainfed
    and irrigated harvested areas at 5 arc-minute spatial resolution using grid
    and administrative masks created for 5 arc-minute resolution in step 1.
    30 arc-minute inputs aggregated in last processing step.
  6) Run scripts in `fertilizer_paper` to generate fertilizer and manure inputs.
    Data processing at 5-arc minute spatial resolution using grid and
    administrative masks created for 5 arc-minute resolution in step 1. Uses
    crop-specific rainfed and irrigated harvested areas at 5 arc-minute spatial
    resolution created in step 5. 30 arc-minute inputs aggregated in last
    processing step. Maximum application thresholds enforced at the final
    resolution (5 arc-minute, 30 arc-minute)

See README files in subdirectories for further detail.
