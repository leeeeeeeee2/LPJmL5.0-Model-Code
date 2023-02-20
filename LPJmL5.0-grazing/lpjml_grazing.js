/**************************************************************************************/
/**                                                                                \n**/
/**                   l  p  j  m  l  .  j  s                                       \n**/
/**                                                                                \n**/
/** Default configuration file for LPJmL C Version 5.1.001                         \n**/
/**                                                                                \n**/
/** Configuration file is divided into five sections:                              \n**/
/**                                                                                \n**/
/**  I.   Simulation description and type section                                  \n**/
/**  II.  Input parameter section                                                  \n**/
/**  III. Input data section                                                       \n**/
/**  IV.  Output data section                                                      \n**/
/**  V.   Run settings section                                                     \n**/
/**                                                                                \n**/
/** (C) Potsdam Institute for Climate Impact Research (PIK), see COPYRIGHT file    \n**/
/** authors, and contributors see AUTHORS file                                     \n**/
/** This file is part of LPJmL and licensed under GNU AGPL Version 3               \n**/
/** or later. See LICENSE file or go to http://www.gnu.org/licenses/               \n**/
/** Contact: https://github.com/PIK-LPJmL/LPJmL                                    \n**/
/**                                                                                \n**/
/**************************************************************************************/

#include "include/conf.h" /* include constant definitions */

//#define DAILY_OUTPUT  /* enables daily output */

/*#define LSU 1    livestock density applied for daily or rotational grazing on mangement grasslands */
/*#define FERT 0  default: 20; fixed fertilizer application rate in gN/m2/yr */
#define MG_SCEN 2 /* grassland management scenario; 0: GS_DEFAULT, 1: GS_MOWING, 2: GS_GRAZING_EXT, 3: GS_GRAZING_INT, 4: GS_NONE} */


//#define OUTFOLDER CONCATENATE(output_lsuha,LSU)


//#define FROM_RESTART

{   /* LPJmL configuration in JSON format */

/*===================================================================*/
/*  I. Simulation description and type section                       */
/*===================================================================*/

  "sim_name" : "LPJmL Run", /* Simulation description */
  "sim_id"   : LPJML,       /* LPJML Simulation type with managed land use */
  "version"  : "5.2",       /* LPJmL version expected */
  "random_prec" : false,     /* Random weather generator for precipitation enabled */
  "random_seed" : 2,        /* seed for random number generator */
  "radiation" : RADIATION,  /* other options: CLOUDINESS, RADIATION, RADIATION_SWONLY, RADIATION_LWDOWN */
  "fire" : FIRE,            /* fire disturbance enabled, other options: NO_FIRE, FIRE, SPITFIRE, SPITFIRE_TMAX */
  "firewood" : false,
  "new_phenology": true,    /* GSI phenology enabled */
  "river_routing" : false,
  "permafrost" : true,
  "const_climate" : false,
  "const_deposition" : false,
#ifdef FROM_RESTART
  "population" : false,
  "landuse" : ONLY_GRASS, //LANDUSE, /* other options: NO_LANDUSE, LANDUSE, CONST_LANDUSE, ALL_CROPS */
  "landuse_year_const" : 2000, /* set landuse year for CONST_LANDUSE case */
  "reservoir" : false,
  "wateruse" : NO_WATERUSE,  /* other options: NO_WATERUSE, WATERUSE, ALL_WATERUSE */
  //"with_nitrogen" : UNLIM_NITROGEN, //AUTO_FERTILIZER, /* other options: NO_NITROGEN, LIM_NITROGEN, UNLIM_NITROGEN */
  "with_nitrogen" : LIM_NITROGEN, /* other options: NO_NITROGEN, LIM_NITROGEN, UNLIM_NITROGEN */
  "equilsoil" : false,
#else
  "population" : false,
  "landuse" : NO_LANDUSE,
  "reservoir" : false,
  "wateruse" : NO_WATERUSE,
  "with_nitrogen" : LIM_NITROGEN, /* other options: NO_NITROGEN, LIM_NITROGEN, UNLIM_NITROGEN */
  "equilsoil" : true,
#endif
  "prescribe_burntarea" : false,
  "prescribe_landcover" : NO_LANDCOVER, /* NO_LANDCOVER, LANDCOVERFPC, LANDCOVEREST */
  "sowing_date_option" : PRESCRIBED_SDATE,   /* NO_FIXED_SDATE, FIXED_SDATE, PRESCRIBED_SDATE */
  "crop_phu_option" : PRESCRIBED_CROP_PHU,    /* PRESCRIBED_CROP_PHU (PHU dataset used, requires PRESCRIBED_SDATE), SEMISTATIC_CROP_PHU (LPJmL4 semi-static PHU approach) */
  "sdate_fixyear" : 1970,               /* year in which sowing dates shall be fixed */
  "intercrop" : false,                   /* intercrops on setaside */
  "black_fallow" : false,
  "no_ndeposition" : false,
  "remove_residuals" : false,           /* remove residuals */
  "residues_fire" : false,              /* fire in residuals */
  "irrigation" : NO_IRRIGATION,        /* NO_IRRIGATION, LIM_IRRIGATION, POT_IRRIGATION, ALL_IRRIGATION */
  "laimax_interpolate" : LAIMAX_PAR, /* laimax values from manage parameter file, */
                                        /* other options: LAIMAX_CFT, CONST_LAI_MAX, LAIMAX_INTERPOLATE, LAIMAX_PAR  */
  "rw_manage" : false,                  /* rain water management */
  "laimax" : 5,                         /* maximum LAI for CONST_LAI_MAX */
  "fertilizer_input" : true,            /* enable fertilizer input */
  "residue_treatment" : READ_RESIDUE_DATA, /* residue options: READ_RESIDUE_DATA, NO_RESIDUE_REMOVE, FIXED_RESIDUE_REMOVE (uses param residues_in_soil) */ 
  "tillage_type" : NO_TILLAGE,          /* Options: TILLAGE (all agr. cells tilled), NO_TILLAGE (no cells tilled) and READ_TILLAGE (tillage dataset used) */
  "manure_input" : false,               /* enable manure input */
  "fix_fertilization" : true,          /* fix fertilizer input */
  "others_to_crop" : false,             /* move PFT type others into PFT crop, maize for tropical, wheat for temperate */
  "grassonly" : false,                  /* set all cropland including others to zero but keep managed grasslands */
  "istimber" : true,
  "grassland_fixed_pft" : false,
  "grass_harvest_options" : false,

/*===================================================================*/
/*  II. Input parameter section                                      */
/*===================================================================*/

#include "par/lpjparam_grazing.js"      /* LPJ parameter file */
#include "par/soil.js"          /* Soil parameter file */
#include "par/pft_phase3_12cft.js"           /* PFT parameter file*/
#include "par/manage_laimax_alphaa_fao_rev4453_20180507.js" /* Management parameter file */
#include "par/manage_reg.js"    /* Management parameter file for regions*/
#include "par/outputvars.js"

/*===================================================================*/
/*  III. Input data section                                          */
/*===================================================================*/

#include "input_grazing.js"    /* Input files of CRU dataset */

/*===================================================================*/
/*  IV. Output data section                                          */
/*===================================================================*/

#ifdef WITH_GRIDBASED
  "pft_output_scaled" : GRIDBASED,
#define SUFFIX grid.bin
#else
  "pft_output_scaled" : PFTBASED,
#define SUFFIX pft.bin
#endif

#define mkstr(s) xstr(s) /* putting string in quotation marks */
#define xstr(s) #s

  "crop_index" : TEMPERATE_CEREALS,  /* CFT for daily output */
  "crop_irrigation" : DAILY_RAINFED, /* irrigation flag for daily output */

#ifdef FROM_RESTART

#ifdef CONSTCO2
  "outpath" : "/p/projects/landuse/users/heinke/Grazing/output_28861_constco2",
#else
  "outpath" : "/p/projects/landuse/users/heinke/Grazing/output_28861",
#endif
  "output" : 
  [

/*
ID                         Fmt                    filename
-------------------------- ---------------------- ----------------------------- */
    { "id" : GRID,               "file" : { "fmt" : RAW, "name" : mkstr(OUTFOLDER/grid.clm) }},
    { "id" : ALITFALLC,          "file" : { "fmt" : RAW, "name" : mkstr(OUTFOLDER/alitterfallc.bin)}},
    { "id" : ALITBURNC,          "file" : { "fmt" : RAW, "name" : mkstr(OUTFOLDER/alitterburnc.bin)}},
    { "id" : VEGC,               "file" : { "fmt" : RAW, "name" : mkstr(OUTFOLDER/vegc.bin)}},
    { "id" : VEGN,               "file" : { "fmt" : RAW, "name" : mkstr(OUTFOLDER/vegn.bin)}},
    { "id" : SOILC,              "file" : { "fmt" : RAW, "name" : mkstr(OUTFOLDER/soilc.bin)}},
    { "id" : LITC,               "file" : { "fmt" : RAW, "name" : mkstr(OUTFOLDER/litc.bin)}},
    { "id" : MEANVEGCMANGRASS,   "file" : { "fmt" : RAW, "name" : mkstr(OUTFOLDER/mean_vegc_mangrass.bin)}},
    { "id" : MGRASS_SOILC,       "file" : { "fmt" : RAW, "name" : mkstr(OUTFOLDER/soilc_mangrass.bin)}},
    { "id" : MGRASS_LITC,        "file" : { "fmt" : RAW, "name" : mkstr(OUTFOLDER/litc_mangrass.bin)}},
    { "id" : MGRASS_SOILN,       "file" : { "fmt" : RAW, "name" : mkstr(OUTFOLDER/soiln_mangrass.bin)}},
    { "id" : MGRASS_LITN,        "file" : { "fmt" : RAW, "name" : mkstr(OUTFOLDER/litn_mangrass.bin)}},
    { "id" : MGRASS_N2O_NIT,     "file" : { "fmt" : RAW, "name" : mkstr(OUTFOLDER/n2o_nit_mangrass.bin)}},
    { "id" : MGRASS_N2O_DENIT,   "file" : { "fmt" : RAW, "name" : mkstr(OUTFOLDER/n2o_denit_mangrass.bin)}},
    { "id" : MGRASS_N2_EMIS,     "file" : { "fmt" : RAW, "name" : mkstr(OUTFOLDER/n2_emis_mangrass.bin)}},
    { "id" : MGRASS_BNF,         "file" : { "fmt" : RAW, "name" : mkstr(OUTFOLDER/bnf_mangrass.bin)}},
    { "id" : MGRASS_LEACHING,    "file" : { "fmt" : RAW, "name" : mkstr(OUTFOLDER/leaching_mangrass.bin)}},
    { "id" : MGRASS_VOLATILIZATION, "file" : { "fmt" : RAW, "name" : mkstr(OUTFOLDER/volatilization_mangrass.bin)}},
    { "id" : CFT_NFERT,          "file" : { "fmt" : RAW, "name" : mkstr(OUTFOLDER/cft_nfert.bin)}},
    { "id" : PFT_HARVESTC,       "file" : { "fmt" : RAW, "name" : mkstr(OUTFOLDER/pft_harvest.SUFFIX)}},
    { "id" : PFT_HARVESTN,       "file" : { "fmt" : RAW, "name" : mkstr(OUTFOLDER/pft_harvestn.SUFFIX)}},
    { "id" : AHARVESTN_MGRASS,   "file" : { "fmt" : RAW, "name" : mkstr(OUTFOLDER/aharvestn_mgrass.bin)}},
    { "id" : ALITFALLN_MGRASS,   "file" : { "fmt" : RAW, "name" : mkstr(OUTFOLDER/alitfalln_mgrass.bin)}},
    { "id" : HARVESTC,           "file" : { "fmt" : RAW, "name" : mkstr(OUTFOLDER/harvestc.bin)}},
    { "id" : HARVESTN,           "file" : { "fmt" : RAW, "name" : mkstr(OUTFOLDER/harvestn.bin)}},
    { "id" : UPTAKEC_MGRASS,     "file" : { "fmt" : RAW, "name" : mkstr(OUTFOLDER/uptakec_mgrass.bin)}},
    { "id" : YIELDC_MGRASS,      "file" : { "fmt" : RAW, "name" : mkstr(OUTFOLDER/yieldc_mgrass.bin)}},
    { "id" : MANUREC_MGRASS,     "file" : { "fmt" : RAW, "name" : mkstr(OUTFOLDER/manurec_mgrass.bin)}},
    { "id" : URINEC_MGRASS,      "file" : { "fmt" : RAW, "name" : mkstr(OUTFOLDER/urinec_mgrass.bin)}},
    { "id" : METHANEC_MGRASS,    "file" : { "fmt" : RAW, "name" : mkstr(OUTFOLDER/methanec_mgrass.bin)}},
    { "id" : RESPC_MGRASS,       "file" : { "fmt" : RAW, "name" : mkstr(OUTFOLDER/respc_mgrass.bin)}},
    { "id" : UPTAKEN_MGRASS,     "file" : { "fmt" : RAW, "name" : mkstr(OUTFOLDER/uptaken_mgrass.bin)}},
    { "id" : YIELDN_MGRASS,      "file" : { "fmt" : RAW, "name" : mkstr(OUTFOLDER/yieldn_mgrass.bin)}},
    { "id" : MANUREN_MGRASS,     "file" : { "fmt" : RAW, "name" : mkstr(OUTFOLDER/manuren_mgrass.bin)}},
    { "id" : URINEN_MGRASS,      "file" : { "fmt" : RAW, "name" : mkstr(OUTFOLDER/urinen_mgrass.bin)}},
    { "id" : ADELTA_NORG_SOIL_MGRASS, "file" : { "fmt" : RAW, "name" : mkstr(OUTFOLDER/adelta_norg_soil_mgrass.bin)}},
    { "id" : ADELTA_NMIN_SOIL_MGRASS, "file" : { "fmt" : RAW, "name" : mkstr(OUTFOLDER/adelta_nmin_soil_mgrass.bin)}},
    { "id" : ADELTA_NVEG_SOIL_MGRASS, "file" : { "fmt" : RAW, "name" : mkstr(OUTFOLDER/adelta_nveg_soil_mgrass.bin)}},
    { "id" : ANUPTAKE_MGRASS,    "file" : { "fmt" : RAW, "name" : mkstr(OUTFOLDER/anuptake_mgrass.bin)}},
    { "id" : ANMINERALIZATION_MGRASS, "file" : { "fmt" : RAW, "name" : mkstr(OUTFOLDER/anmineralization_mgrass.bin)}},
    { "id" : ANIMMOBILIZATION_MGRASS, "file" : { "fmt" : RAW, "name" : mkstr(OUTFOLDER/animmobilization_mgrass.bin)}},
    { "id" : ASEEDN_MGRASS,      "file" : { "fmt" : RAW, "name" : mkstr(OUTFOLDER/aseedn_mgrass.bin)}},
    { "id" : ANFERT_MGRASS,      "file" : { "fmt" : RAW, "name" : mkstr(OUTFOLDER/afert_mgrass.bin)}},
    { "id" : ANDEPO_MGRASS,      "file" : { "fmt" : RAW, "name" : mkstr(OUTFOLDER/adepo_mgrass.bin)}},
    { "id" : MN2O_DENIT,         "file" : { "fmt" : RAW, "name" : mkstr(OUTFOLDER/mn2o_denit.bin)}},
    { "id" : MN2O_NIT,           "file" : { "fmt" : RAW, "name" : mkstr(OUTFOLDER/mn2o_nit.bin)}},
    { "id" : CFT_AIRRIG,         "file" : { "fmt" : RAW, "name" : mkstr(OUTFOLDER/cft_airrig.SUFFIX)}},
    { "id" : MGRASS_GPP,         "file" : { "fmt" : RAW, "name" : mkstr(OUTFOLDER/gpp_mgrass.bin)}},
    { "id" : MGRASS_RA,          "file" : { "fmt" : RAW, "name" : mkstr(OUTFOLDER/ra_mgrass.bin)}},
    { "id" : MGRASS_RH,          "file" : { "fmt" : RAW, "name" : mkstr(OUTFOLDER/rh_mgrass.bin)}},
    { "id" : MGRASS_DELTAC,      "file" : { "fmt" : RAW, "name" : mkstr(OUTFOLDER/deltac_mgrass.bin)}},
    { "id" : MGRASS_SEEDC,       "file" : { "fmt" : RAW, "name" : mkstr(OUTFOLDER/seedc_mgrass.bin)}},
    { "id" : MGRASS_AGBMEAN,     "file" : { "fmt" : RAW, "name" : mkstr(OUTFOLDER/agbmean_mgrass.bin)}},
#ifdef DAILY_OUTPUT
    { "id" : D_PHU,              "file" : { "fmt" : RAW, "name" : mkstr(OUTFOLDER/d_phu.bin)}},
    { "id" : D_FPHU,             "file" : { "fmt" : RAW, "name" : mkstr(OUTFOLDER/d_fphu.bin)}}
#endif
/*------------------------ ---------------------- ------------------------------- */
  ],

#else

  "output" : [],  /* no output written */

#endif

/*===================================================================*/
/*  V. Run settings section                                          */
/*===================================================================*/

  "startgrid" : 28861, /* 27410, 67208 60400 all grid cells */
  "endgrid" : 28861,

#ifdef CONSTCO2
  "restartpath" : "/p/projects/landuse/users/heinke/Grazing/output_28861_constco2",
#else
  "restartpath" : "/p/projects/landuse/users/heinke/Grazing/output_28861",
#endif
    
#ifdef CHECKPOINT
  "checkpoint_filename" : "restart/restart_checkpoint.lpj", /* filename of checkpoint file */
#endif

#ifndef FROM_RESTART

  "nspinup" : 7000,  /* spinup years */
  "nspinyear" : 30,  /* cycle length during spinup (yr) */
  "firstyear": 1901, /* first year of simulation */
  "lastyear" : 1901, /* last year of simulation */
  "restart" : false, /* do not start from restart file */
  "write_restart" : true, /* create restart file: the last year of simulation=restart-year */
  "write_restart_filename" : "restart_1900_nv_stdfire.lpj.lpj", /* filename of restart file */
  "restart_year": 1900 /* write restart at year */

#else

  "nspinup" : 3000, // 390,   /* spinup years */
  "nspinyear" : 30,  /* cycle length during spinup (yr)*/
  "firstyear": 1901, /* first year of simulation */
  "lastyear" : 2016, /* last year of simulation */
  "outputyear": 1901, /* first year output is written  */
  "restart" :  true, /* start from restart file */
  "restart_filename" : "restart_1900_nv_stdfire.lpj.lpj", /* filename of restart file */
  "write_restart" : false, /* create restart file */
  "write_restart_filename" : "", /* filename of restart file */
  "restart_year": 1900 /* write restart at year */

#endif
}
