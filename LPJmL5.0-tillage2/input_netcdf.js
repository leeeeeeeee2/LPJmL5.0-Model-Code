/**************************************************************************************/
/**                                                                                \n**/
/**        i  n  p  u  t  _  n  e  t  c  d  f  .  j  s                             \n**/
/**                                                                                \n**/
/** Configuration file for NetCDF input dataset for LPJmL                          \n**/
/** Version 5.1.001                                                                \n**/
/**                                                                                \n**/
/** (C) Potsdam Institute for Climate Impact Research (PIK), see COPYRIGHT file    \n**/
/** authors, and contributors see AUTHORS file                                     \n**/
/** This file is part of LPJmL and licensed under GNU AGPL Version 3               \n**/
/** or later. See LICENSE file or go to http://www.gnu.org/licenses/               \n**/
/** Contact: https://github.com/PIK-LPJmL/LPJmL                                    \n**/
/**                                                                                \n**/
/**************************************************************************************/

#include "include/conf.h" /* include constant definitions */

"inpath" : "/p/projects/biodiversity",

"input" :
{
  "soil" : { "fmt" : CDF, "var" : "soilcode", "name" : "cru_netcdf/soil_new_67420.nc"},
  "no3deposition" : { "fmt" : CDF, "var" : "no3dep",  "name" : "cru_netcdf/no3_deposition_rcp8p5.nc"},
  "nh4deposition" : { "fmt" : CDF,  "var" : "nh4dep", "name" : "cru_netcdf/nh4_deposition_rcp8p5.nc"},
  "soilpH" :        { "fmt" : CDF,  "var" : "soilph", "name" : "cru_netcdf/soil_ph.nc"},
  "countrycode" : { "fmt" : CDF, "var" : "country", "name" : "cru_netcdf/cow_full_2018.nc"},
  "regioncode" : { "fmt" : CDF, "var" : "region", "name" : "cru_netcdf/reg_full_2018.nc"},
  "landuse" : { "fmt" : CDF, "var" :  "landfrac", "name" : "cru_netcdf/cft1700_2005_irrigation_systems_64bands.nc"},
  "fertilizer_nr" : { "fmt" : CDF,  "var" : "nfert", "name" : "cru_netcdf/fertilizer_ggcmi.nc"},
  "temp" : { "fmt" : CDF, "var" : "temp", "name" : "cru_netcdf/cru_ts_3_10.1901.2009.tmp.nc"},
  "prec" : { "fmt" : CDF, "var" : "prec", "name" : "cru_netcdf/gpcc_cru09_prec_monthly_1901_2009.nc"},
  "cloud" : { "fmt" : CDF, "var" : "cloud", "name" : "cru_netcdf/cru_ts_3_10.1901.2009.cld.nc"},
  "wind" : { "fmt" : CDF, "var" : "windspeed", "name" : "cru_netcdf/mwindspeed_1860-2100_67420.nc"},
  "tamp" : { "fmt" : CDF, "var" : "tamp", "name" : "cru_netcdf/cru_ts_3_10.1901.2009.dtr.nc"},           /* diurnal temp. range */
  "lightning" : { "fmt" : CDF, "var" : "lightning", "name" : "cru_netcdf/lightning.nc"},
  "human_ignition" : { "fmt" : CDF, "var" : "human_ignition", "name" : "cru_netcdf/human_ignition.nc"},
  "popdens" : { "fmt" : CDF, "var" : "popdens", "name" : "cru_netcdf/popdens_HYDE_1901_2010_bi.nc "},
  "co2" : { "fmt" : TXT, "name" : "input_VERSION2/co2_1841-2010.dat"},

  "wetdays" : { "fmt" : CDF, "var" : "wet", "name" : "cru_netcdf/gpcc_cru09_wet_monthly_1901_2009.nc"},
  "wateruse" : { "fmt" : CDF, "var" : "wateruse", "name" : "cru_netcdf/wateruse_1900_2000.nc"} /* water consumption for industry, household and livestock */
},
