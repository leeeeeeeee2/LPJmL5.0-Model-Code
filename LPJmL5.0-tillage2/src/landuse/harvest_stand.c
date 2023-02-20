/**************************************************************************************/
/**                                                                                \n**/
/**         h  a  r  v  e  s  t  _  s  t  a  n  d  .  c                            \n**/
/**                                                                                \n**/
/**     C implementation of LPJmL                                                  \n**/
/**                                                                                \n**/
/**     Function harvests grassland stand                                          \n**/
/**                                                                                \n**/
/** (C) Potsdam Institute for Climate Impact Research (PIK), see COPYRIGHT file    \n**/
/** authors, and contributors see AUTHORS file                                     \n**/
/** This file is part of LPJmL and licensed under GNU AGPL Version 3               \n**/
/** or later. See LICENSE file or go to http://www.gnu.org/licenses/               \n**/
/** Contact: https://github.com/PIK-LPJmL/LPJmL                                    \n**/
/**                                                                                \n**/
/**************************************************************************************/

#include "lpj.h"
#include "grass.h"
#include "agriculture.h"
#include "grassland.h"

Harvest harvest_grass(Stand *stand, /**< pointer to stand */
                      Real hfrac    /**< harvest fractions */
                     )              /** \return harvested grass (gC/m2) */
{
  Harvest harvest;
  Harvest sum={{0,0},{0,0},{0,0},{0,0}};
  Pftgrass *grass;
  Pft *pft;
  int p;

  foreachpft(pft,p,&stand->pftlist)
  {
    grass=pft->data;
    harvest.harvest.carbon=grass->ind.leaf.carbon*hfrac;
    harvest.harvest.nitrogen=grass->ind.leaf.nitrogen*hfrac*0.25; /*0.25*/
    pft->stand->soil.NH4[0]+=grass->ind.leaf.nitrogen*hfrac*0.75;
    grass->ind.leaf.carbon*=(1-hfrac);
    grass->ind.leaf.nitrogen*=(1-hfrac);
    sum.harvest.carbon+=harvest.harvest.carbon;
    sum.harvest.nitrogen+=harvest.harvest.nitrogen;
    grass->max_leaf = grass->ind.leaf.carbon;
    pft->phen=1; /*0.3;*/
    pft->gdd=30;
  }
  return sum;
} /* of 'harvest_grass' */

static Harvest harvest_grass_mowing(Stand *stand)
{
  Harvest harvest;
  Harvest sum={{0,0},{0,0},{0,0},{0,0}};
  Pftgrass *grass;
  Pft *pft;
  int p;

  foreachpft(pft,p,&stand->pftlist)
  {
    grass=pft->data;
    harvest.harvest.carbon = grass->ind.leaf.carbon - STUBBLE_HEIGHT_MOWING;
    harvest.harvest.nitrogen = harvest.harvest.carbon/grass->ind.leaf.carbon*grass->ind.leaf.nitrogen;

    grass->ind.leaf.carbon = STUBBLE_HEIGHT_MOWING;
    grass->ind.leaf.nitrogen -= harvest.harvest.nitrogen;
    sum.harvest.carbon+=harvest.harvest.carbon;
    sum.harvest.nitrogen+=harvest.harvest.nitrogen;
    pft->gdd=pft->phen=0.0; // change -> relative from ind.leaf
  }
  return sum;
} /* of 'harvest_grass_mowing' */

static Harvest harvest_grass_grazing_ext(Stand *stand)
{
  Harvest sum={{0,0},{0,0},{0,0},{0,0}};
  Pftgrass *grass;
  Pft *pft;
  int p;
  Real bm_grazed;
  Real fact;
  Stocks bm_tot =  {0.0,0.0};
  Stocks bm_grazed_pft;

  foreachpft(pft,p,&stand->pftlist)
  {
    grass=pft->data;
    bm_tot.carbon += grass->ind.leaf.carbon;
    bm_tot.nitrogen+= grass->ind.leaf.nitrogen;
  }
//  bm_grazed = stand->cell->ml.nr_of_lsus_ext * DEMAND_COW_EXT; 
  bm_grazed = 1e-4* stand->cell->ml.nr_of_lsus_ext * DEMAND_COW_EXT;

  foreachpft(pft,p,&stand->pftlist)
  {
    grass=pft->data;
    if (bm_tot.carbon < 1e-5) // to avoid division by zero!
      fact = 1;
    else
      fact = grass->ind.leaf.carbon/bm_tot.carbon;

    bm_grazed_pft.carbon   = bm_grazed * fact;
    if (grass->ind.leaf.carbon - bm_grazed_pft.carbon < GRAZING_STUBBLE)
      bm_grazed_pft.carbon = grass->ind.leaf.carbon - GRAZING_STUBBLE;
    if (bm_grazed_pft.carbon < 0)
      bm_grazed_pft.carbon = 0;

    pft->gdd = (1-(bm_grazed_pft.carbon/grass->ind.leaf.carbon)) * pft->gdd;


    /* Nitrogen */
    //bm_grazed_pft.nitrogen = bm_grazed * fact;
    bm_grazed_pft.nitrogen = bm_grazed_pft.carbon/grass->ind.leaf.carbon*grass->ind.leaf.nitrogen;

    grass->ind.leaf.carbon -= bm_grazed_pft.carbon;
    sum.harvest.carbon     += (1-MANURE)*bm_grazed_pft.carbon;                       // 60% atmosphere, 15% cows
    stand->soil.pool->fast.carbon += MANURE * bm_grazed_pft.carbon;             // 25% back to soil

    grass->ind.leaf.nitrogen -=  bm_grazed_pft.nitrogen;
    sum.harvest.nitrogen     += (1-MANURE)*bm_grazed_pft.nitrogen;                       // 60% atmosphere, 15% cows
    stand->soil.pool->fast.nitrogen += MANURE * bm_grazed_pft.nitrogen;             // 25% back to soil 
    // pft->phen recalculated in phenology_grass
  }
  return sum;
} /* of 'harvest_grass_grazing_ext' */

static Harvest harvest_grass_grazing_int(Stand *stand)
{
  Harvest sum={{0,0},{0,0},{0,0},{0,0}};
  Pftgrass *grass;
  Pft *pft;
  int p;
  Real rotation_len;
  Real fact;
  Real bm_grazed;
  Stocks bm_tot = {0,0};
  Stocks bm_grazed_pft;
  Rotation *rotation;

  rotation = &(stand->cell->ml.rotation);
  foreachpft(pft,p,&stand->pftlist)
  {
    grass=pft->data;
    bm_tot.carbon += grass->ind.leaf.carbon;
    bm_tot.nitrogen += grass->ind.leaf.nitrogen;
  }

  if (rotation->rotation_mode == RM_UNDEFINED) //initial calculate grazing days and recovery days
  {
    rotation_len = (bm_tot.carbon - GRAZING_STUBBLE) / (1e4*stand->cell->ml.nr_of_lsus_int * DEMAND_COW_INT) ;
    if (rotation_len > MAX_ROTATION_LENGTH)
      rotation_len = MAX_ROTATION_LENGTH;

    if (rotation_len > MIN_ROTATION_LENGTH) // otherwise wait for more growth 
    {
      rotation->grazing_days = (int)ceil(rotation_len/MAX_PADDOCKS);
      rotation->paddocks = (int)floor((rotation_len/rotation->grazing_days) + 0.5);
      rotation->recovery_days = (rotation->paddocks-1) * rotation->grazing_days;
      rotation->rotation_mode = RM_GRAZING;
    }
  }

  if (rotation->rotation_mode == RM_GRAZING)
  {
    bm_grazed = stand->cell->ml.nr_of_lsus_int * DEMAND_COW_INT * rotation->paddocks;
    foreachpft(pft,p,&stand->pftlist)
    {
      grass=pft->data;
      fact = grass->ind.leaf.carbon / bm_tot.carbon;
      bm_grazed_pft.carbon = bm_grazed * fact;

      if (grass->ind.leaf.carbon - bm_grazed_pft.carbon < GRAZING_STUBBLE)
        bm_grazed_pft.carbon = grass->ind.leaf.carbon - GRAZING_STUBBLE;

      if (bm_grazed_pft.carbon < 0)
        bm_grazed_pft.carbon =0;

      pft->gdd = (1-(bm_grazed_pft.carbon/grass->ind.leaf.carbon)) * pft->gdd;

      /* Nitrogen */
      bm_grazed_pft.nitrogen = bm_grazed_pft.carbon/grass->ind.leaf.carbon*grass->ind.leaf.nitrogen;

      grass->ind.leaf.carbon -= bm_grazed_pft.carbon;
      sum.harvest.carbon     += (1-MANURE)*bm_grazed_pft.carbon;              // 60% atmosphere, 15% cows
      stand->soil.pool->fast.carbon += MANURE * bm_grazed_pft.carbon;    // 25% back to soil

      grass->ind.leaf.nitrogen -= bm_grazed_pft.nitrogen;
      sum.harvest.nitrogen     += (1-MANURE)*bm_grazed_pft.nitrogen;              // 60% atmosphere, 15% cows
      stand->soil.pool->fast.nitrogen += MANURE * bm_grazed_pft.nitrogen;    // 25% back to soil
    }

    rotation->grazing_days -= 1;
    if (rotation->grazing_days == 0)
      rotation->rotation_mode = RM_RECOVERY;
  }
  else if (rotation->rotation_mode == RM_RECOVERY)
  {
    rotation->recovery_days -= 1;
    if (rotation->recovery_days == 0)
      rotation->rotation_mode = RM_UNDEFINED;
  }
  return sum;
} /* of 'harvest_grass_grazing_int' */


Harvest harvest_stand(Output *output, /**< Output data */
                      Stand *stand,   /**< pointer to grassland stand */
                      Real hfrac      /**< harvest fraction */
                     )                /** \return harvested carbon (gC/m2) */
{
  Harvest harvest;
  if (stand->type->landusetype == GRASSLAND)
  {
    switch (stand->cell->ml.grass_scenario)
    {
      case GS_DEFAULT: // default
        harvest=harvest_grass(stand,hfrac);
        break;
      case GS_MOWING: // mowing
        harvest=harvest_grass_mowing(stand);
        break;
      case GS_GRAZING_EXT: // ext. grazing
        harvest=harvest_grass_grazing_ext(stand);
        break;
      case GS_GRAZING_INT: // int. grazing
        harvest=harvest_grass_grazing_int(stand);
        break;
    }
  }
  else /* option for biomass_grass */
  {
    harvest=harvest_grass(stand,hfrac);
  }
  output->flux_harvest.carbon+=(harvest.harvest.carbon+harvest.residual.carbon)*stand->frac;
  output->flux_harvest.nitrogen+=(harvest.harvest.nitrogen+harvest.residual.nitrogen)*stand->frac;
  output->dcflux+=(harvest.harvest.carbon+harvest.residual.carbon)*stand->frac;
  stand->growing_days=0;
  return harvest;
}
 /* of 'harvest_stand' */
