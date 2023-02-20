/**************************************************************************************/
/**                                                                                \n**/
/**         h  a  r  v  e  s  t  _  c  r  o  p  .  c                               \n**/
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
#include "agriculture.h"
#include "crop.h"

/*
 *  Called in function daily_agriculture when crop
 *  is harvested
 *
 */
void harvest_crop(Output *output,        /**< Output data */
                  Stand *stand,          /**< pointer to crop stand */
                  Pft *pft,              /**< PFT variables */
                  int npft,              /**< number of natural PFTs */
                  int ncft,              /**< number of crop PFTs */
                  int year,              /* year of harvest */
                  int residue_treatment, /** treatment of residues */
                  Real agrfrac,          /** agricultural cell fraction */
                  Bool residues_fire,    /**< fire in residuals after harvest (TRUE/FALSE) */
                  Bool pft_output_scaled /**< pft-specific output scaled with
                                              stand->frac (TRUE/FALSE) */
                 )
{
  Pftcroppar *croppar;
  Harvest harvest;
  Pftcrop *crop;
  Irrigation *data;
  Real fuelratio,bifratio,factor;
  data=stand->data;
  crop=pft->data;
  Real res_onfield, res_remove;

  if(residue_treatment<READ_RESIDUE_DATA)
    res_onfield = residue_treatment==FIXED_RESIDUE_REMOVE ? param.residues_in_soil : 1 ;
  else
    res_onfield=stand->cell->ml.residue_on_field[data->irrigation].crop[pft->par->id-npft];
  if (year<param.till_startyear)
    res_onfield=0.1;
  res_remove = (1-res_onfield);
  stand->soil.litter.item[pft->litter].ag.leaf.carbon += (crop->ind.leaf.carbon + crop->ind.pool.carbon)*res_onfield;
  stand->cell->output.alittfall_agr.carbon += (crop->ind.leaf.carbon + crop->ind.pool.carbon)*res_onfield*stand->frac/agrfrac;
  stand->soil.litter.item[pft->litter].ag.leaf.nitrogen += (crop->ind.leaf.nitrogen + crop->ind.pool.nitrogen)*res_onfield;
  if (!residues_fire)
  {
    harvest.residuals_burnt.carbon = harvest.residuals_burntinfield.carbon =
      harvest.residuals_burnt.nitrogen = harvest.residuals_burntinfield.nitrogen = 0;
  }
  else
  {
    fuelratio = stand->cell->ml.manage.regpar->fuelratio; /* burn outside of field */
    bifratio = stand->cell->ml.manage.regpar->bifratio; /* burn in field */
    if (bifratio + fuelratio > res_remove)
    {
      bifratio *= res_remove;
      fuelratio *= res_remove;
    }
    res_remove -= fuelratio - bifratio;
    harvest.residuals_burnt.carbon = (crop->ind.leaf.carbon + crop->ind.pool.carbon)*fuelratio;
    harvest.residuals_burntinfield.carbon = (crop->ind.leaf.carbon + crop->ind.pool.carbon)*bifratio;
    harvest.residuals_burnt.nitrogen = (crop->ind.leaf.nitrogen + crop->ind.pool.nitrogen)*fuelratio;
    harvest.residuals_burntinfield.nitrogen = (crop->ind.leaf.nitrogen + crop->ind.pool.nitrogen)*bifratio;
  }
  harvest.residual.carbon = (crop->ind.leaf.carbon + crop->ind.pool.carbon)*res_remove;
  harvest.residual.nitrogen = (crop->ind.leaf.nitrogen + crop->ind.pool.nitrogen)*res_remove;
  harvest.harvest=crop->ind.so;
  stand->soil.litter.item[pft->litter].bg.carbon+=crop->ind.root.carbon;
  stand->cell->output.alittfall.carbon+=crop->ind.root.carbon*stand->frac;
  stand->cell->output.alittfall_agr.carbon+= crop->ind.root.carbon*stand->frac/agrfrac;
  stand->soil.litter.item[pft->litter].bg.nitrogen+=crop->ind.root.nitrogen;
  stand->cell->output.alittfall.nitrogen+=crop->ind.root.nitrogen*stand->frac;
#ifdef DOUBLE_HARVEST
  if(pft_output_scaled)
  {
    double_harvest(output->syear2[pft->par->id-npft+data->irrigation*ncft],
      &(output->pft_harvest[pft->par->id-npft+data->irrigation*(ncft+NGRASS+NBIOMASSTYPE)].harvest.carbon),
      &(output->pft_harvest2[pft->par->id-npft+data->irrigation*(ncft+NGRASS+NBIOMASSTYPE)].harvest.carbon),harvest.harvest.carbon*stand->frac);
    double_harvest(output->syear2[pft->par->id-npft+data->irrigation*ncft],
      &(output->pft_harvest[pft->par->id-npft+data->irrigation*(ncft+NGRASS+NBIOMASSTYPE)].harvest.nitrogen),
      &(output->pft_harvest2[pft->par->id-npft+data->irrigation*(ncft+NGRASS+NBIOMASSTYPE)].harvest.nitrogen),harvest.harvest.nitrogen*stand->frac);
    double_harvest(output->syear2[pft->par->id-npft+data->irrigation*ncft],
      &(output->pft_harvest[pft->par->id-npft+data->irrigation*(ncft+NGRASS+NBIOMASSTYPE)].residual.carbon),
      &(output->pft_harvest2[pft->par->id-npft+data->irrigation*(ncft+NGRASS+NBIOMASSTYPE)].residual.carbon),(harvest.residual.carbon+harvest.residuals_burnt.carbon+harvest.residuals_burntinfield.carbon)*stand->frac);
    double_harvest(output->syear2[pft->par->id-npft+data->irrigation*ncft],
      &(output->pft_harvest[pft->par->id-npft+data->irrigation*(ncft+NGRASS+NBIOMASSTYPE)].residual.nitrogen),
      &(output->pft_harvest2[pft->par->id-npft+data->irrigation*(ncft+NGRASS+NBIOMASSTYPE)].residual.nitrogen),(harvest.residual.nitrogen+harvest.residuals_burnt.nitrogen+harvest.residuals_burntinfield.nitrogen)*stand->frac);
  }
  else
  {
    double_harvest(output->syear2[pft->par->id-npft+data->irrigation*ncft],
      &(output->pft_harvest[pft->par->id-npft+data->irrigation*(ncft+NGRASS+NBIOMASSTYPE)].harvest.carbon),
      &(output->pft_harvest2[pft->par->id-npft+data->irrigation*(ncft+NGRASS+NBIOMASSTYPE)].harvest.carbon),harvest.harvest.carbon);
    double_harvest(output->syear2[pft->par->id-npft+data->irrigation*ncft],
      &(output->pft_harvest[pft->par->id-npft+data->irrigation*(ncft+NGRASS+NBIOMASSTYPE)].harvest.nitrogen),
      &(output->pft_harvest2[pft->par->id-npft+data->irrigation*(ncft+NGRASS+NBIOMASSTYPE)].harvest.nitrogen),harvest.harvest.nitrogen);
    double_harvest(output->syear2[pft->par->id-npft+data->irrigation*ncft],
      &(output->pft_harvest[pft->par->id-npft+data->irrigation*(ncft+NGRASS+NBIOMASSTYPE)].residual.carbon),
      &(output->pft_harvest2[pft->par->id-npft+data->irrigation*(ncft+NGRASS+NBIOMASSTYPE)].residual.carbon),harvest.residual.carbon+harvest.residuals_burnt.carbon+harvest.residuals_burntinfield.carbon);
    double_harvest(output->syear2[pft->par->id-npft+data->irrigation*ncft],
      &(output->pft_harvest[pft->par->id-npft+data->irrigation*(ncft+NGRASS+NBIOMASSTYPE)].residual.nitrogen),
      &(output->pft_harvest2[pft->par->id-npft+data->irrigation*(ncft+NGRASS+NBIOMASSTYPE)].residual.nitrogen),harvest.residual.nitrogen+harvest.residuals_burnt.nitrogen+harvest.residuals_burntinfield.nitrogen);
  }
  /* harvested area */
  double_harvest(output->syear2[pft->par->id-npft+data->irrigation*ncft],
    output->cftfrac+pft->par->id-npft+data->irrigation*(ncft+NGRASS+NBIOMASSTYPE),
    output->cftfrac2+pft->par->id-npft+data->irrigation*(ncft+NGRASS+NBIOMASSTYPE),stand->frac);
  if(output->syear2[pft->par->id-npft+data->irrigation*ncft]>0)
    output->sdate2[pft->par->id-npft+data->irrigation*ncft]=crop->sdate;
  else
    output->sdate[pft->par->id-npft+data->irrigation*ncft]=crop->sdate;
#else
  if(pft_output_scaled)
  {
    output->pft_harvest[pft->par->id-npft+data->irrigation*(ncft+NGRASS+NBIOMASSTYPE)].harvest.carbon+=harvest.harvest.carbon*stand->frac;
    output->pft_harvest[pft->par->id-npft+data->irrigation*(ncft+NGRASS+NBIOMASSTYPE)].harvest.nitrogen+=harvest.harvest.nitrogen*stand->frac;
    output->pft_harvest[pft->par->id-npft+data->irrigation*(ncft+NGRASS+NBIOMASSTYPE)].residual.carbon+=(harvest.residual.carbon+harvest.residuals_burnt.carbon+harvest.residuals_burntinfield.carbon)*stand->frac;
    output->pft_harvest[pft->par->id-npft+data->irrigation*(ncft+NGRASS+NBIOMASSTYPE)].residual.nitrogen+=(harvest.residual.nitrogen+harvest.residuals_burnt.nitrogen+harvest.residuals_burntinfield.nitrogen)*stand->frac;
  }
  else
  {
    output->pft_harvest[pft->par->id-npft+data->irrigation*(ncft+NGRASS+NBIOMASSTYPE)].harvest.carbon+=harvest.harvest.carbon;
    output->pft_harvest[pft->par->id-npft+data->irrigation*(ncft+NGRASS+NBIOMASSTYPE)].harvest.nitrogen+=harvest.harvest.nitrogen;
    output->pft_harvest[pft->par->id-npft+data->irrigation*(ncft+NGRASS+NBIOMASSTYPE)].residual.carbon+=(harvest.residual.carbon+harvest.residuals_burnt.carbon+harvest.residuals_burntinfield.carbon);
    output->pft_harvest[pft->par->id-npft+data->irrigation*(ncft+NGRASS+NBIOMASSTYPE)].residual.nitrogen+=(harvest.residual.nitrogen+harvest.residuals_burnt.nitrogen+harvest.residuals_burntinfield.nitrogen);
  }
  /* harvested area */
  output->cftfrac[pft->par->id-npft+data->irrigation*(ncft+NGRASS+NBIOMASSTYPE)]=stand->frac;
#endif

  output->cft_leaf[pft->par->id-npft+data->irrigation*(ncft+NGRASS+NBIOMASSTYPE)].nitrogen+=crop->ind.leaf.nitrogen;
  output->cft_nlimit[pft->par->id-npft+data->irrigation*(ncft+NGRASS+NBIOMASSTYPE)]+=pft->nlimit;
  output->cft_leaf[pft->par->id-npft+data->irrigation*(ncft+NGRASS+NBIOMASSTYPE)].carbon+=crop->ind.leaf.carbon;
  output->cft_root[pft->par->id-npft+data->irrigation*(ncft+NGRASS+NBIOMASSTYPE)].nitrogen+=crop->ind.root.nitrogen;
  output->cft_root[pft->par->id-npft+data->irrigation*(ncft+NGRASS+NBIOMASSTYPE)].carbon+=crop->ind.root.carbon;
  output->cft_veg[pft->par->id-npft+data->irrigation*(ncft+NGRASS+NBIOMASSTYPE)].nitrogen+=vegn_sum(pft);
  output->cft_veg[pft->par->id-npft+data->irrigation*(ncft+NGRASS+NBIOMASSTYPE)].carbon+=vegc_sum(pft);

  output->flux_harvest.carbon+=(harvest.harvest.carbon+harvest.residual.carbon+harvest.residuals_burnt.carbon+harvest.residuals_burntinfield.carbon)*stand->frac;
  output->flux_harvest.nitrogen+=(harvest.harvest.nitrogen+harvest.residual.nitrogen+harvest.residuals_burnt.nitrogen+harvest.residuals_burntinfield.nitrogen)*stand->frac;
  output->dcflux+=(harvest.harvest.carbon+harvest.residual.carbon+harvest.residuals_burnt.carbon+harvest.residuals_burntinfield.carbon)*stand->frac;
  output->flux_rharvest_burnt.carbon+=harvest.residuals_burnt.carbon*stand->frac;
  output->flux_rharvest_burnt.nitrogen+=harvest.residuals_burnt.nitrogen*stand->frac;
  output->flux_rharvest_burnt_in_field.carbon+=harvest.residuals_burntinfield.carbon*stand->frac;
  output->flux_rharvest_burnt_in_field.nitrogen+=harvest.residuals_burntinfield.nitrogen*stand->frac;
  croppar=pft->par->data;
  if(data->irrigation)
    stand->cell->ml.cropdates[pft->par->id-npft].fallow_irrig=croppar->fallow_days;
  else
    stand->cell->ml.cropdates[pft->par->id-npft].fallow=croppar->fallow_days;
} /* of 'harvest_crop' */
