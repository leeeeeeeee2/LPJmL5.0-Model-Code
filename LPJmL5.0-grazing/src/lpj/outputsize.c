/**************************************************************************************/
/**                                                                                \n**/
/**                o  u  t  p  u  t  s  i  z  e  .  c                              \n**/
/**                                                                                \n**/
/**     C implementation of LPJmL                                                  \n**/
/**                                                                                \n**/
/**     Function calculates number of items per cell in output file                \n**/
/**                                                                                \n**/
/** (C) Potsdam Institute for Climate Impact Research (PIK), see COPYRIGHT file    \n**/
/** authors, and contributors see AUTHORS file                                     \n**/
/** This file is part of LPJmL and licensed under GNU AGPL Version 3               \n**/
/** or later. See LICENSE file or go to http://www.gnu.org/licenses/               \n**/
/** Contact: https://github.com/PIK-LPJmL/LPJmL                                    \n**/
/**                                                                                \n**/
/**************************************************************************************/

#include "lpj.h"

int outputsize(int index,    /**< output index */
               int npft,     /**< number of natural PFTs */
               int nbiomass, /**< number of biomass PFTs */
               int ncft      /**< number of crop PFTs */
              )              /** \return number of items per cell */
{
  switch(index)
  {
    case SDATE: case HDATE: case HDATE2: case SDATE2: case CFT_MSWC:
    case SYEAR: case SYEAR2: case HUSUM: case HUSUM2: case CFT_LEACHING:
    case CFT_N2O_DENIT : case CFT_N2O_NIT: case CFT_N2_EMIS: case CFT_C_EMIS:
    case CFT_N2O_DENIT2: case CFT_N2O_NIT2: case CFT_N2_EMIS2: case CFT_C_EMIS2:
    case CFT_LEACHING2: case CFT_RUNOFF: case CFT_RUNOFF2:
      return ncft*2;
    case PFT_NPP: case PFT_GCGP: case PFT_LAIMAX: case PFT_NLIMIT:
    case PFT_NUPTAKE: case PFT_NDEMAND: case PFT_VEGC: case PFT_VEGN:
    case PFT_CLEAF: case PFT_NLEAF:
    case PFT_CROOT: case PFT_NROOT: case PFT_CSAPW: case PFT_NSAPW:
    case PFT_CHAWO: case PFT_NHAWO:
      return npft-nbiomass+(ncft+NGRASS+NBIOMASSTYPE)*2;
    case PFT_HARVESTC: case PFT_RHARVESTC: 
    case PFT_HARVESTN: case PFT_RHARVESTN: 
    case CFT_CONSUMP_WATER_G: 
    case CFT_CONV_LOSS_EVAP: case CFT_CONV_LOSS_DRAIN:
    case CFT_CONSUMP_WATER_B: case CFTFRAC: case CFT_AIRRIG: case CFT_FPAR: 
    case CFT_RETURN_FLOW_B:
    case LUC_IMAGE: case CFT_INTERC: case CFT_INTERC_B: case CFT_NIR: 
    case CFT_TRANSP: case CFT_TRANSP_B:
    case CFT_EVAP: case CFT_EVAP_B: case CFT_IRRIG_EVENTS:
    case PFT_HARVESTC2: case PFT_RHARVESTC2:
    case PFT_HARVESTN2: case PFT_RHARVESTN2:
    case CFT_INTERC2: 
    case CFTFRAC2: case CFT_AIRRIG2: 
    case CFT_TRANSP2: case CFT_NIR2:
    case CFT_EVAP2:
      return (ncft+NGRASS+NBIOMASSTYPE)*2;
    case FPC:
      return npft-nbiomass+1;
    case MSOILTEMP: case MSWC:
      return NSOILLAYER;
    case SOILC_LAYER: case SOILN_LAYER: case SOILNO3_LAYER: case SOILNH4_LAYER: case SOILC_AGR_LAYER:
    case RESPONSE_LAYER_AGR: case RESPONSE_LAYER_NV: case CSHIFT_FAST_NV: case CSHIFT_SLOW_NV:
      return LASTLAYER;
    case GROWING_PERIOD: case CFT_TEMP:case CFT_PREC:
    case CFT_SRAD: case CFT_ABOVEGBMC: case CFT_ABOVEGBMN:
    case GROWING_PERIOD2:case CFT_TEMP2:case CFT_PREC2:
    case CFT_SRAD2: case CFT_ABOVEGBMC2: case CFT_ABOVEGBMN2:
    case CFT_PET: case CFT_PET2:
      return (ncft+NGRASS)*2;
    default:
      return 1;
  } /* of 'switch' */
} /* of 'outputsize' */
