/**************************************************************************************/
/**                                                                                \n**/
/**       i  n  p  u  t  _  c  r  u  m  o  n  t  h  l  y  .  j  s                  \n**/
/**                                                                                \n**/
/** Configuration file for input dataset for LPJ C Version 4.0.002                 \n**/
/**                                                                                \n**/
/** (C) Potsdam Institute for Climate Impact Research (PIK), see COPYRIGHT file    \n**/
/** authors, and contributors see AUTHORS file                                     \n**/
/** This file is part of LPJmL and licensed under GNU AGPL Version 3               \n**/
/** or later. See LICENSE file or go to http://www.gnu.org/licenses/               \n**/
/** Contact: https://github.com/PIK-LPJmL/LPJmL                                    \n**/
/**                                                                                \n**/
/**************************************************************************************/

#include "include/conf.h" /* include constant definitions */


"inpath" : "/p/projects/lpjml/input",

"input" :
{
//  "soil" :         { "fmt" : META, "name" : "ISIMIP3/soil_ISISMIP3.descr"},
  "soil" :         { "fmt" : META, "name" : "historical/input_VERSION2/soil.descr"},
//  "soil" :         { "fmt" : META, "name" : "historical/input_VERSION2/soil_climasteppe.descr"},
 // "soil" :         { "fmt" : CLM, "name" : "/home/rolinski/LPJ/Climasteppe/Preprocessing/soil_climasteppe.bin"},
  "coord" :        { "fmt" : CLM,  "name" : "ISIMIP3/grid.bin"},
  "countrycode" :  { "fmt" : CLM,  "name" : "historical/input_VERSION2/cow_full_2018.bin"},
  "no3deposition" : { "fmt" : CLM,  "name" : "ISIMIP3/no3_deposition_2015soc_1850-2100.clm"},
  "nh4deposition" : { "fmt" : CLM,  "name" : "ISIMIP3/nh4_deposition_2015soc_1850-2100.clm"},
  "soilpH" :        { "fmt" : CLM,  "name" : "historical/input_VERSION2/soil_ph.clm"},
  "lakes" :        { "fmt" : META, "name" : "historical/input_VERSION2/glwd_lakes_and_rivers.descr"},
  "drainage" :     { "fmt" : CLM,  "name" : "historical/input_VERSION2/drainagestn.bin"},
  "neighb_irrig" : { "fmt" : CLM,  "name" : "historical/input_VERSION2/neighb_irrig_stn.bin"},
  "elevation" :    { "fmt" : CLM,  "name" : "historical/input_VERSION2/elevation.bin"},
  "reservoir" :    { "fmt" : CLM,  "name" : "historical/input_VERSION2/reservoir_info_grand5.bin"},
  "landuse" :      { "fmt" : CLM,  "name" : "historical/input_VERSION2/cft1700_2005_irrigation_systems_64bands.bin"},
  "fertilizer_nr" : { "fmt" : CLM,  "name" : "historical/input_VERSION3/fertilizer_luh2v2_1900-2015_32bands.clm"},
  "manure_nr" :    { "fmt" : CLM, "name" : "historical/input_VERSION3/manure_zhang17_1860-2014_32bands_clm2.clm"},
  "sdate" :        {"fmt" : CLM, "name" : "crop_calendar/sdates_ggcmi_phase3_v1.01_67420_24bands.clm"},  /* insert prescribed sdate file name here */
  "crop_phu" :     {"fmt" : CLM, "name" : "crop_calendar/phu_agmerra_ggcmi_phase3_v1.01_67420_24bands.clm"},  /* insert prescribed phu file name here */
  "with_tillage" : { "fmt" : CLM, "name" : "historical/input_VERSION3/lpj_tillage_CA_1973-2010.clm"},
  "residue_on_field" : { "fmt" : CLM, "name" : "/p/projects/lpjml/input/MADRAT/residues_madrat_1850-2015_16bands.clm"},
  "temp" :         { "fmt" : CLM,  "name" : "historical/ISIMIP3a/obsclim/GSWP3-W5E5/tas_gswp3-w5e5_obsclim_1901-2016.clm"},
  "tmax" :         { "fmt" : CLM,  "name" : "historical/ISIMIP3a/obsclim/GSWP3-W5E5/tasmax_gswp3-w5e5_obsclim_1901-2016.clm"},
  "tmin" :         { "fmt" : CLM,  "name" : "historical/ISIMIP3a/obsclim/GSWP3-W5E5/tasmin_gswp3-w5e5_obsclim_1901-2016.clm"},
  "prec" :         { "fmt" : CLM,  "name" : "historical/ISIMIP3a/obsclim/GSWP3-W5E5/pr_gswp3-w5e5_obsclim_1901-2016.clm"},
  "lwnet" :        { "fmt" : CLM,  "name" : "historical/ISIMIP3a/obsclim/GSWP3-W5E5/lwnet_gswp3-w5e5_obsclim_1901-2016.clm"},
  "swdown" :       { "fmt" : CLM,  "name" : "historical/ISIMIP3a/obsclim/GSWP3-W5E5/rsds_gswp3-w5e5_obsclim_1901-2016.clm"},
  "wind":          { "fmt" : CLM,  "name" : "historical/ISIMIP3a/obsclim/GSWP3-W5E5/sfcwind_gswp3-w5e5_obsclim_1901-2016.clm"},
#ifdef CONSTCO2
  "co2" :          { "fmt" : TXT,  "name" : "/p/projects/landuse/users/heinke/Grazing/lpjml_internal/co2_const278ppm_1850_2018.txt"},
#else
  "co2" :          { "fmt" : TXT,  "name" : "historical/ISIMIP3a/co2/co2_obsclim_annual_1850_2018.txt"},
#endif
},


