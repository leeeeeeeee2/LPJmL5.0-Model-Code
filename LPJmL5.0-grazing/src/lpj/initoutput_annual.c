/**************************************************************************************/
/**                                                                                \n**/
/**         i  n  i  t  o  u  t  p  u  t  _  a  n  n  u  a  l  .  c                \n**/
/**                                                                                \n**/
/**     C implementation of LPJmL                                                  \n**/
/**                                                                                \n**/
/**     Function initializes annual output data to zero                            \n**/
/**                                                                                \n**/
/** (C) Potsdam Institute for Climate Impact Research (PIK), see COPYRIGHT file    \n**/
/** authors, and contributors see AUTHORS file                                     \n**/
/** This file is part of LPJmL and licensed under GNU AGPL Version 3               \n**/
/** or later. See LICENSE file or go to http://www.gnu.org/licenses/               \n**/
/** Contact: https://github.com/PIK-LPJmL/LPJmL                                    \n**/
/**                                                                                \n**/
/**************************************************************************************/

#include "lpj.h"

void initoutput_annual(Output *output, /**< Output data */
                       int npft,       /**< number of natural PFTs */
                       int nbiomass,   /**< number of biomass PFTs */
                       int ncft        /**< number of crop PFTs */
                      )
{
  int i;
  output->neg_fluxes.carbon= output->neg_fluxes.nitrogen=
  output->fire.carbon=output->fire.nitrogen=output->firef=
  output->flux_harvest.carbon=output->flux_harvest.nitrogen=
  output->flux_estab.carbon=output->flux_estab.nitrogen=0;
  output->input_lake=output->flux_firewood.carbon=output->flux_firewood.nitrogen=
  output->flux_rharvest_burnt.carbon=output->flux_rharvest_burnt.nitrogen=
  output->flux_rharvest_burnt_in_field.carbon= output->flux_rharvest_burnt_in_field.nitrogen=0;
  output->atransp=output->aevap=output->ainterc=output->airrig=output->aconv_loss_evap=output->aconv_loss_drain=output->awateruse_hil=0;
  output->awd_unsustainable=output->aevap_lake=output->aevap_res=0;
  output->soil_storage=output->aburntarea=0;
  output->soil_storage=output->alittfall.carbon=output->alittfall.nitrogen=0;
  output->runoff_surf=output->runoff_lat=output->anpp = output->anpp_agr = output->arh = output->arh_agr = 0;
  output->abnf_agr=output->anfert_agr=output->anmanure_agr=output->andepo_agr=output->anmineralization_agr=output->animmobilization_agr=
  output->anuptake_agr=output->anleaching_agr=output->an2o_denit_agr=output->an2o_nit_agr=output->anh3_agr=output->an2_agr=
  output->andepo_mgrass=output->anfert_mgrass=output->aharvestn_mgrass=output->alitfalln_mgrass=
  output->animmobilization_mgrass=output->anmineralization_mgrass=output->anuptake_mgrass=output->aseedn_mgrass=
  output->alitfalln_agr=output->aharvestn_agr=output->aseedn_agr=output->adelta_norg_soil_agr=output->adelta_nmin_soil_agr=
  output->adelta_norg_soil_mgrass=output->adelta_nmin_soil_mgrass=output->adelta_nveg_soil_mgrass=
  output->adelta_nveg_soil_agr=output->cellfrac_agr=output->mgrass_n2o_denit=output->mgrass_n2o_nit=output->mgrass_n2_emis=0;
  output->mgrass_leaching=output->mgrass_bnf=output->mgrass_volatilization=0;
  output->flux_nfert=0;
  output->deforest_emissions.carbon=output->deforest_emissions.nitrogen=output->fburn=output->ftimber=output->timber_harvest.carbon=output->timber_harvest.nitrogen=0;
  output->trad_biofuel=output->mean_vegc_mangrass=0;
  output->alittfall_wood.carbon=output->alittfall_wood.nitrogen=
  output->alitburnc=output->alitburnc_wood=0;
  output->yieldc_mgrass=output->manurec_mgrass=output->urinec_mgrass=
  output->yieldn_mgrass=output->manuren_mgrass=output->urinen_mgrass=
  output->uptakec_mgrass=output->uptaken_mgrass=
  output->deficit_lsu_mp=output->deficit_lsu_ne=
  output->methanec_mgrass=output->respc_mgrass=0.0;
  output->mgrass_gpp=output->mgrass_ra=output->mgrass_rh=output->mgrass_deltac=output->mgrass_seedc=0;
  output->mgrass_agbmean=0;

  /* these need to be initialized to one! */
  output->decay_wood_agr=output->decay_wood_nv=output->decay_leaf_agr=output->decay_leaf_nv=1;

#ifdef IMAGE
  output->prod_turnover=0;
  output->product_pool_fast=output->product_pool_slow=0;
#else
  output->prod_turnover.carbon=output->prod_turnover.nitrogen=0;
  output->product_pool_fast.carbon=output->product_pool_slow.carbon=output->product_pool_fast.nitrogen=output->product_pool_slow.nitrogen=0;
#endif

  /* memory allocation now in newgrid.c */

  for(i=0;i<NSOILLAYER;i++)
    output->response_agr[i]=output->response_nv[i]=
    output->cshift_fast_nv[i]=output->cshift_slow_nv[i]=0;

  for(i=0;i<(ncft+NGRASS+NBIOMASSTYPE)*2;i++)
    output->pft_harvest[i].harvest.carbon=output->pft_harvest[i].residual.carbon=
    output->pft_harvest[i].harvest.nitrogen=output->pft_harvest[i].residual.nitrogen=
    output->cftfrac[i]=output->cft_laimax[i]=
    output->cft_consump_water_g[i]=output->cft_consump_water_b[i]=
    output->cft_transp[i]=output->cft_transp_b[i]=output->cft_evap[i]=output->cft_evap_b[i]=
    output->cft_interc[i]=output->cft_interc_b[i]=output->cft_return_flow_b[i]=output->cft_nir[i]=
    output->cft_nfert[i]=
#ifdef DOUBLE_HARVEST
    output->pft_harvest2[i].harvest.carbon=output->pft_harvest2[i].residual.carbon=
    output->pft_harvest2[i].harvest.nitrogen=output->pft_harvest2[i].residual.nitrogen=
    output->cftfrac2[i]=
    output->cft_transp2[i]=output->cft_evap2[i]=output->cft_interc2[i]=
    output->cft_nir2[i]=output->cft_airrig2[i]=output->cft_nfert2[i]=
#endif
    output->cft_nlimit[i]=output->cft_laimax[i]=
    output->cft_airrig[i]=output->cft_fpar[i]=output->cft_luc_image[i]=output->cft_conv_loss_evap[i]=output->cft_conv_loss_drain[i]=
    output->cft_leaf[i].nitrogen=output->cft_leaf[i].carbon=output->cft_root[i].carbon=output->cft_root[i].nitrogen=
    output->cft_veg[i].carbon=output->cft_veg[i].nitrogen=output->cft_irrig_events[i]=0;
  for(i=0;i<(ncft+NGRASS)*2;i++)
    output->growing_period[i]=output->cft_pet[i]=
    output->cft_temp[i]=output->cft_prec[i]=output->cft_srad[i]=
#ifdef DOUBLE_HARVEST
    output->growing_period2[i]=output->cft_pet2[i]=
    output->cft_temp2[i]=output->cft_prec2[i]=output->cft_srad2[i]=
    output->cft_aboveground_biomass2[i].carbon=output->cft_aboveground_biomass2[i].nitrogen=
#endif
    output->cft_aboveground_biomass[i].carbon=output->cft_aboveground_biomass[i].nitrogen=0;
  for(i=0;i<(ncft*2);i++){
#ifdef DOUBLE_HARVEST
    output->sdate2[i]=output->hdate2[i]=
    output->syear[i]=output->syear2[i]=0;
    output->husum2[i]=0.;
    output->cft_runoff2[i]=output->cft_n2o_denit2[i]=output->cft_n2o_nit2[i]=output->cft_n2_emis2[i]=
      output->cft_c_emis2[i]=output->cft_leaching2[i]=0.0;
#endif
    output->sdate[i]=output->hdate[i]=0;
    output->husum[i]=0;
    output->cft_runoff[i]=output->cft_n2o_denit[i]=output->cft_n2o_nit[i]=output->cft_n2_emis[i]=
      output->cft_c_emis[i]=output->cft_leaching[i]=0.0;
  }
  for(i=0;i<(ncft+NGRASS+NBIOMASSTYPE)*2+npft-nbiomass;i++)
  {
    output->pft_npp[i]=0;
    output->pft_gcgp[i]=0;
    output->gcgp_count[i]=0;
    output->pft_nuptake[i]=0;
    output->pft_ndemand[i]=0;
#ifdef DOUBLE_HARVEST
    output->pft_nuptake2[i]=0;
#endif
  }
  for (i=0; i<npft-nbiomass+1;++i)
    output->fpc[i] = 0;
} /* of 'initoutput_annual' */
