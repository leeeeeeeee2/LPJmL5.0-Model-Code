/**************************************************************************************/
/**                                                                                \n**/
/**              s  o  i  l  .  j  s                                               \n**/
/**                                                                                \n**/
/**  Soil parameter data for LPJmL version 5.2.001                                 \n**/
/**                                                                                \n**/
/** (C) Potsdam Institute for Climate Impact Research (PIK), see COPYRIGHT file    \n**/
/** authors, and contributors see AUTHORS file                                     \n**/
/** This file is part of LPJmL and licensed under GNU AGPL Version 3               \n**/
/** or later. See LICENSE file or go to http://www.gnu.org/licenses/               \n**/
/** Contact: https://github.com/PIK-LPJmL/LPJmL                                    \n**/
/**                                                                                \n**/
/**************************************************************************************/

#include "../include/soilpar.h"

"soildepth" :
[
  200.,       /* soil depth layer (mm) */
  300.,       /* soil depth layer (mm) */
  500.,       /* soil depth layer (mm) */
  1000.,      /* soil depth layer (mm) */
  1000.,      /* soil depth layer (mm) */
  10000.      /* soil depth layer (mm) */
],

"fbd_fac" : [1., 1.2, 1.4, 0.], /* fuel bulk density factors */

"soilpar" :
[

/*
ID           Name               Ks     Sf    %sand %silt   %clay hsg  tdiff_0  tdiff_15  tdiff_100  cond_pwp  cond_100  cond_100_ice
----------  ------------------- ----   ----  ------ ------ ----- ---  -------  --------  ---------  --------  --------  -------------*/
{ "id" : Cl,     "name" : "clay",            "Ks" : 3.5,   "Sf" : 243., "sand" : 0.22,  "silt" : 0.20, "clay": 0.58, "hsg" : 1, "tdiff_0" : 0.572, "tdiff_15" : 0.571, "tdiff_100" : 0.555, "cond_pwp" : 1.388, "cond_100" : 1.719, "cond_100_ice" : 3.233, "denit_rate": 6, "a_denit": -73.5, "b_denit": 0.4804, "anion_excl": 0.4, "nitvol_temp_factor": 0.41, "vol_cation_exchange": 0.15, "denit_water_threshold": 2, "a_nit": 0.6, "b_nit": 1.27, "c_nit": 0.0012, "d_nit": 2.84,"cn_ratio": 15},
{ "id" : SiCl,   "name" : "silty clay",      "Ks" : 4.8,   "Sf" : 169., "sand" : 0.06,  "silt" : 0.47, "clay": 0.47, "hsg" : 2, "tdiff_0" : 0.502, "tdiff_15" : 0.503, "tdiff_100" : 0.491, "cond_pwp" : 1.177, "cond_100" : 1.516, "cond_100_ice" : 2.853, "denit_rate": 6, "a_denit": -73.5, "b_denit": 0.4804, "anion_excl": 0.4, "nitvol_temp_factor": 0.41, "vol_cation_exchange": 0.15, "denit_water_threshold": 2, "a_nit": 0.6, "b_nit": 1.27,"c_nit": 0.0012, "d_nit": 2.84,"cn_ratio": 15},
{ "id" : SaCl,   "name" : "sandy clay",      "Ks" : 26.0,  "Sf" :  51., "sand" : 0.52,  "silt" : 0.06, "clay": 0.42, "hsg" : 4, "tdiff_0" : 0.785, "tdiff_15" : 0.791, "tdiff_100" : 0.799, "cond_pwp" : 1.720, "cond_100" : 2.347, "cond_100_ice" : 4.060, "denit_rate": 6, "a_denit": -73.5, "b_denit": 0.4804, "anion_excl": 0.4, "nitvol_temp_factor": 0.41, "vol_cation_exchange": 0.15, "denit_water_threshold": 2, "a_nit": 0.55, "b_nit": 1.7,"c_nit": -0.007, "d_nit": 3.22,"cn_ratio": 15},
{ "id" : ClLo,   "name" : "clay loam",       "Ks" : 8.8,   "Sf" : 139., "sand" : 0.32,  "silt" : 0.34, "clay": 0.34, "hsg" : 3, "tdiff_0" : 0.650, "tdiff_15" : 0.656, "tdiff_100" : 0.653, "cond_pwp" : 1.369, "cond_100" : 1.967, "cond_100_ice" : 3.685, "denit_rate": 6, "a_denit": -73.5, "b_denit": 0.4804, "anion_excl": 0.4, "nitvol_temp_factor": 0.41, "vol_cation_exchange": 0.15, "denit_water_threshold": 2, "a_nit": 0.6, "b_nit": 1.27, "c_nit": 0.0012, "d_nit": 2.84,"cn_ratio": 15},
{ "id" : SiClLo, "name" : "silty clay loam", "Ks" : 7.3,   "Sf" : 324., "sand" : 0.10,  "silt" : 0.56, "clay": 0.34, "hsg" : 2, "tdiff_0" : 0.556, "tdiff_15" : 0.557, "tdiff_100" : 0.542, "cond_pwp" : 1.270, "cond_100" : 1.675, "cond_100_ice" : 3.134, "denit_rate": 6, "a_denit": -73.5, "b_denit": 0.4804, "anion_excl": 0.4, "nitvol_temp_factor": 0.41, "vol_cation_exchange": 0.15, "denit_water_threshold": 2, "a_nit": 0.6, "b_nit": 1.27, "c_nit": 0.0012, "d_nit": 2.84,"cn_ratio": 15},
{ "id" : SaClLo, "name" : "sandy clay loam", "Ks" : 16.0,  "Sf" :  72., "sand" : 0.58,  "silt" : 0.15, "clay": 0.27, "hsg" : 1, "tdiff_0" : 0.780, "tdiff_15" : 0.808, "tdiff_100" : 0.867, "cond_pwp" : 1.498, "cond_100" : 2.527, "cond_100_ice" : 4.360, "denit_rate": 6, "a_denit": -73.5, "b_denit": 0.4804, "anion_excl": 0.4, "nitvol_temp_factor": 0.41, "vol_cation_exchange": 0.15, "denit_water_threshold": 2, "a_nit": 0.55, "b_nit": 1.7, "c_nit": -0.007, "d_nit": 3.22,"cn_ratio": 15},
{ "id" : Lo,     "name" : "loam",            "Ks" : 12.2,  "Sf" : 192., "sand" : 0.43,  "silt" : 0.39, "clay": 0.18, "hsg" : 4, "tdiff_0" : 0.701, "tdiff_15" : 0.740, "tdiff_100" : 0.797, "cond_pwp" : 1.276, "cond_100" : 2.340, "cond_100_ice" : 4.233, "denit_rate": 6, "a_denit": -73.5, "b_denit": 0.4804, "anion_excl": 0.4, "nitvol_temp_factor": 0.41, "vol_cation_exchange": 0.15, "denit_water_threshold": 2, "a_nit": 0.6, "b_nit": 1.27, "c_nit": 0.0012, "d_nit": 2.84,"cn_ratio": 15},
{ "id" : SiLo,   "name" : "silt loam",       "Ks" : 10.1,  "Sf" : 409., "sand" : 0.17,  "silt" : 0.70, "clay": 0.13, "hsg" : 2, "tdiff_0" : 0.637, "tdiff_15" : 0.657, "tdiff_100" : 0.661, "cond_pwp" : 1.219, "cond_100" : 1.999, "cond_100_ice" : 3.803, "denit_rate": 6, "a_denit": -73.5, "b_denit": 0.4804, "anion_excl": 0.4, "nitvol_temp_factor": 0.41, "vol_cation_exchange": 0.15, "denit_water_threshold": 2, "a_nit": 0.6, "b_nit": 1.27, "c_nit": 0.0012, "d_nit": 2.84,"cn_ratio": 15},
{ "id" : SaLo,   "name" : "sandy loam",      "Ks" : 18.8,  "Sf" :  77., "sand" : 0.58,  "silt" : 0.32, "clay": 0.10, "hsg" : 2, "tdiff_0" : 0.640, "tdiff_15" : 0.713, "tdiff_100" : 0.863, "cond_pwp" : 1.053, "cond_100" : 2.530, "cond_100_ice" : 4.547, "denit_rate": 6, "a_denit": -73.5, "b_denit": 0.4804, "anion_excl": 0.4, "nitvol_temp_factor": 0.41, "vol_cation_exchange": 0.15, "denit_water_threshold": 2, "a_nit": 0.55, "b_nit": 1.7, "c_nit": -0.007, "d_nit": 3.22,"cn_ratio": 15},
{ "id" : Si,     "name" : "silt",            "Ks" : 10.1,  "Sf" : 409., "sand" : 0.10,  "silt" : 0.60, "clay": 0.30, "hsg" : 2, "tdiff_0" : 0.637, "tdiff_15" : 0.657, "tdiff_100" : 0.661, "cond_pwp" : 1.219, "cond_100" : 1.999, "cond_100_ice" : 3.803, "denit_rate": 6, "a_denit": -73.5, "b_denit": 0.4804, "anion_excl": 0.4, "nitvol_temp_factor": 0.41, "vol_cation_exchange": 0.15, "denit_water_threshold": 2, "a_nit": 0.6, "b_nit": 1.27, "c_nit": 0.0012, "d_nit": 2.84, "cn_ratio": 15},
{ "id" : LoSa,   "name" : "loamy sand",      "Ks" : 50.7,  "Sf" :  20., "sand" : 0.82, "silt" : 0.12, "clay": 0.06, "hsg" : 2, "tdiff_0" : 0.403, "tdiff_15" : 0.529, "tdiff_100" : 0.850, "cond_pwp" : 0.601, "cond_100" : 2.706, "cond_100_ice" : 4.778, "denit_rate": 6, "a_denit": -73.5, "b_denit": 0.4804, "anion_excl": 0.4, "nitvol_temp_factor": 0.41, "vol_cation_exchange": 0.15, "denit_water_threshold": 2, "a_nit": 0.55, "b_nit": 1.7, "c_nit": -0.007, "d_nit": 3.22,"cn_ratio": 15},
{ "id" : Sa,     "name" : "sand",            "Ks" : 167.8, "Sf" :  39., "sand" : 0.92,  "silt" : 0.05, "clay": 0.03, "hsg" : 2, "tdiff_0" : 0.201, "tdiff_15" : 0.196, "tdiff_100" : 0.896, "cond_pwp" : 0.303, "cond_100" : 3.431, "cond_100_ice" : 5.423, "denit_rate": 6, "a_denit": -73.5, "b_denit": 0.4804, "anion_excl": 0.4, "nitvol_temp_factor": 0.41, "vol_cation_exchange": 0.15, "denit_water_threshold": 2, "a_nit": 0.55, "b_nit": 1.7, "c_nit": -0.007, "d_nit": 3.22,"cn_ratio": 15},
{ "id" : ROCK,   "name" : "rock and ice",    "Ks" : 0.1,   "Sf" :   1., "sand" : 0.99, "silt" : 0.00, "clay": 0.01, "hsg" : 4, "tdiff_0" : 4.137, "tdiff_15" : 4.127, "tdiff_100" : 4.089, "cond_pwp" : 8.768, "cond_100" : 8.657, "cond_100_ice" : 8.727, "denit_rate": 6, "a_denit": -73.5, "b_denit": 0.4804, "anion_excl": 0.4, "nitvol_temp_factor": 0.41, "vol_cation_exchange": 0.15, "denit_water_threshold": 2, "a_nit": 0.6, "b_nit": 1.27, "c_nit": 0.0012, "d_nit": 2.84,"cn_ratio": 15}
],
/*
----------  ------------------- ----   ----  ------ ------ ----- ---  -------  --------  ---------  --------  --------  -------------*/

/* Ks: saturated hydraulic conductivity (mm/h) following Cosby (1984)
 * Sf: Suction head (mm) in Green-Ampt equation following Rawls, Brakensiek and Miller (1982)
 * w_pwp: water content at wilting point following Cosby (1984)
 * w_fc: water content at field cpacity following Cosby (1984)
 * w_sat: water content at saturation following Cosby (1984)
 * tdiff0: thermal diffusivity (mm^2/s) at wilting point (0% whc) following Lawrence and Slater (2008)
 * tdiff15: thermal diffusivity (mm^2/s) at 15% whc following Lawrence and Slater (2008)
 * tdiff100: thermal diffusivity (mm^2/s) at field capacity (100% whc) following Lawrence and Slater (2008)
 * cond_pwp: thermal conductivity (W/m^2/K) at wilting point following Lawrence and Slater (2008)
 * cond_100: thermal conductivity (W/m^2/K) at saturation (all water) following Lawrence and Slater (2008)
 * cond_100_ice: thermal conductivity (W/m^2/K) at saturation (all ice) Lawrence and Slater (2008)
 */
