/**************************************************************************************/
/**                                                                                \n**/
/**        f  r  e  e  o  u  t  p  u  t  .  c                                      \n**/
/**                                                                                \n**/
/**     C implementation of LPJmL                                                  \n**/
/**                                                                                \n**/
/** (C) Potsdam Institute for Climate Impact Research (PIK), see COPYRIGHT file    \n**/
/** authors, and contributors see AUTHORS file                                     \n**/
/** This file is part of LPJmL and licensed under GNU AGPL Version 3               \n**/
/** or later. See LICENSE file or go to http://www.gnu.org/licenses/               \n**/
/** Contact: https://github.com/PIK-LPJmL/LPJmL                                    \n**/
/**                                                                                \n**/
/**************************************************************************************/

#include "lpj.h"

void freeoutput(Output *output /**< Output data */
               )
{
  free(output->sdate);
  free(output->hdate);
  free(output->husum);
  free(output->pft_harvest);
  free(output->cft_consump_water_g);
  free(output->cft_consump_water_b);
  free(output->growing_period);  
  free(output->pft_npp);
  free(output->pft_nuptake);
  free(output->pft_ndemand);
  free(output->pft_gcgp);
  free(output->gcgp_count);
  free(output->fpc);
  free(output->cftfrac);
  free(output->cft_pet);
  free(output->cft_transp);
  free(output->cft_transp_b);
  free(output->cft_evap);
  free(output->cft_evap_b);
  free(output->cft_interc);
  free(output->cft_interc_b);
  free(output->cft_return_flow_b);
  free(output->cft_nir);
  free(output->cft_fpar);
  free(output->cft_temp);
  free(output->cft_prec);
  free(output->cft_srad);
  free(output->cft_aboveground_biomass);
  free(output->cft_airrig);
  free(output->cft_conv_loss_evap);
  free(output->cft_conv_loss_drain);
  free(output->cft_luc_image);
  free(output->cft_irrig_events);
  free(output->cft_leaf);
  free(output->cft_root);
  free(output->cft_veg);
  free(output->cft_nlimit);
  free(output->cft_laimax);
  free(output->cft_mswc);
  free(output->nday_month);
  free(output->cft_runoff);
  free(output->cft_n2o_denit);
  free(output->cft_n2o_nit);
  free(output->cft_n2_emis);
  free(output->cft_leaching);
  free(output->cft_c_emis);
  free(output->cft_nfert);
#ifdef DOUBLE_HARVEST
  free(output->cftfrac2);
  free(output->sdate2);
  free(output->hdate2);
  free(output->husum2);
  free(output->pft_harvest2);
  free(output->cft_pet2);
  free(output->cft_transp2);
  free(output->cft_evap2);
  free(output->cft_interc2);
  free(output->cft_nir2);
  free(output->cft_temp2);
  free(output->cft_prec2);
  free(output->cft_srad2);
  free(output->cft_aboveground_biomass2);
  free(output->growing_period2);
  free(output->cft_airrig2);
  free(output->syear);
  free(output->syear2);
  free(output->pft_nuptake2);
  free(output->cft_runoff2);
  free(output->cft_n2o_denit2);
  free(output->cft_n2o_nit2);
  free(output->cft_n2_emis2);
  free(output->cft_leaching2);
  free(output->cft_c_emis2);
  free(output->cft_nfert2);

  output->sdate2=output->hdate2=output->syear=output->syear2=NULL;
  output->husum2=NULL;
  output->cft_transp2=output->cft_evap2=output->cft_interc2=output->cft_nir2=
    output->cft_pet2=output->cftfrac2=output->cft_airrig2=NULL;
  output->pft_harvest2=NULL;
  output->cft_temp2=output->cft_prec2=output->cft_srad2=NULL;
  output->cft_aboveground_biomass2=NULL;
  output->cft_runoff2=output->cft_n2o_denit2=output->cft_n2o_nit2=output->cft_n2_emis2=
    output->cft_leaching2=output->cft_c_emis2=output->cft_nfert2=NULL;

#endif
  output->sdate=output->hdate=NULL; 
  output->husum=NULL;
  output->pft_harvest=NULL;
  output->cft_consump_water_g=output->cft_consump_water_b=output->cft_pet=NULL;
  output->cft_transp=output->cft_transp_b=output->cft_evap=output->cft_evap_b=output->cft_interc=output->cft_interc_b=
    output->cft_return_flow_b=output->cft_nir=output->cft_temp=output->cft_prec=output->cft_srad=output->cft_fpar=NULL;
  output->cft_aboveground_biomass=NULL;
  output->cft_conv_loss_evap=output->cft_conv_loss_drain=NULL;
  output->growing_period=NULL;
  output->cft_irrig_events=NULL;
  output->pft_npp=output->fpc=output->cftfrac=output->cft_airrig=output->cft_luc_image=NULL;
  output->cft_mswc=NULL;
  output->cft_runoff=output->cft_n2o_denit=output->cft_n2o_nit=output->cft_n2_emis=
    output->cft_leaching=output->cft_c_emis=output->cft_nfert=NULL;
} /* of 'freeoutput' */
