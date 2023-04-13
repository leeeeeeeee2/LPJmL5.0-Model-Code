# AQUASTAT data

This file describes the AQUASTAT data used by the toolbox.

AQUASTAT is released by the Food and Agriculture Organization of the United
Nations (FAO).

## Input
AQUASTAT data can be downloaded from: 
http://www.fao.org/aquastat/statistics/query/index.html

## How to use
1. In the top-left field select the variables you would like to downloaded.
  The variables of interest are found under:
  `All variables` -> `Irrigation and drainage development` -> `Irrigated crop
  area and cropping intensity`.
  There select all crop-specific `Harvested irrigated temporary crop area: 
  [[crop name]]` variables and all crop-specific `Harvested irrigated permanent
  crop area: [[crop name]]` variables
2. In the top-right field select `All Countries`. If you don't plan to generate
  a global landuse dataset you may select individual countries of interest.
3. In the middle field select the period for which to download data. We suggest
  to select the full period.
4. Click on `Submit` button to retrieve data
5. Once the data view opens click on the `CSV (Flat)` download button in the
  top-right corner.

## Data for sample application
Since AQUASTAT data are updated regularly and previous versions are no longer
available a copy of the exact AQUASTAT data used for the sample application is
included in this archive (`aquastat_irrigated_harvested_areas.csv`).

Data credit: FAO.AQUASTAT Core Database.
License: [CC BY-NC-SA 3.0 IGO](https://creativecommons.org/licenses/by-nc-sa/3.0/igo/).
Extracted from: http://www.fao.org/aquastat/statistics/query/index.html.
Date of Access: 17-04-2020
