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
  Output *output;
  output=&stand->cell->output;

  foreachpft(pft,p,&stand->pftlist)
  {
    grass=pft->data;
    harvest.harvest.carbon=grass->ind.leaf.carbon*hfrac;
    harvest.harvest.nitrogen=grass->ind.leaf.nitrogen*hfrac*0.5; /*0.25*/
    stand->soil.NH4[0]+=grass->ind.leaf.nitrogen*hfrac*0.5*pft->nind;
    grass->ind.leaf.carbon*=(1-hfrac);
    grass->ind.leaf.nitrogen*=(1-hfrac);
    stand->soil.litter.item[pft->litter].bg.carbon+=grass->ind.root.carbon*hfrac*param.rootreduction*pft->nind;
    output->alittfall.carbon+=grass->ind.root.carbon*hfrac*param.rootreduction*pft->nind*stand->frac;
    stand->soil.litter.item[pft->litter].bg.nitrogen+=grass->ind.root.nitrogen*hfrac*param.rootreduction*pft->nind*pft->par->fn_turnover;
    output->alittfall.nitrogen+=grass->ind.root.nitrogen*hfrac*param.rootreduction*pft->nind*stand->frac*pft->par->fn_turnover;
    output->alitfalln_mgrass+=grass->ind.root.nitrogen*hfrac*param.rootreduction*pft->nind*stand->frac*pft->par->fn_turnover;
    pft->bm_inc.nitrogen+=grass->ind.root.nitrogen*hfrac*param.rootreduction*pft->nind*pft->nind*(1-pft->par->fn_turnover);
    grass->ind.root.carbon*=(1-hfrac*param.rootreduction);
    grass->ind.root.nitrogen*=(1-hfrac*param.rootreduction);
    sum.harvest.carbon+=harvest.harvest.carbon*pft->nind;
    sum.harvest.nitrogen+=harvest.harvest.nitrogen*pft->nind;
    grass->max_leaf = grass->ind.leaf.carbon;
    pft->phen=1; /*0.3;*/
    pft->gdd=30;
  }
  return sum;
} /* of 'harvest_grass' */

Harvest harvest_grass_mowing(Stand *stand)
{
  Harvest harvest;
  Harvest sum={{0,0},{0,0},{0,0},{0,0}};
  Pftgrass *grass;
  Pft *pft;
  int p;
  Real fpc_sum=0.0;
  Real hfrac;
  Output *output;
  output=&stand->cell->output;
  foreachpft(pft,p,&stand->pftlist)
    fpc_sum+=pft->fpc;

  foreachpft(pft,p,&stand->pftlist)
  {
    grass=pft->data;
    harvest.harvest.carbon = grass->ind.leaf.carbon - STUBBLE_HEIGHT_MOWING*pft->fpc/fpc_sum;
    hfrac=harvest.harvest.carbon/grass->ind.leaf.carbon;
    harvest.harvest.nitrogen = hfrac*grass->ind.leaf.nitrogen*0.25;
    stand->soil.NH4[0]+=hfrac*grass->ind.leaf.nitrogen*0.75*pft->nind;
    stand->soil.litter.item[pft->litter].bg.carbon+=grass->ind.root.carbon*hfrac*param.rootreduction*pft->nind;
    output->alittfall.carbon+=grass->ind.root.carbon*hfrac*param.rootreduction*pft->nind*stand->frac;
    stand->soil.litter.item[pft->litter].bg.nitrogen+=grass->ind.root.nitrogen*hfrac*param.rootreduction*pft->nind*pft->par->fn_turnover;
    output->alittfall.nitrogen+=grass->ind.root.nitrogen*hfrac*param.rootreduction*pft->nind*stand->frac*pft->par->fn_turnover;
    output->alitfalln_mgrass+=grass->ind.root.nitrogen*hfrac*param.rootreduction*pft->nind*stand->frac*pft->par->fn_turnover;
    pft->bm_inc.nitrogen+=grass->ind.root.nitrogen*hfrac*param.rootreduction*pft->nind*pft->nind*(1-pft->par->fn_turnover);

    grass->ind.root.carbon*=(1-hfrac*param.rootreduction);
    grass->ind.root.nitrogen*=(1-hfrac*param.rootreduction);
    sum.harvest.carbon+=harvest.harvest.carbon*pft->nind;
    sum.harvest.nitrogen+=harvest.harvest.nitrogen*pft->nind;
    pft->gdd=0.0; // change -> relative from ind.leaf
  }
  return sum;
} /* of 'harvest_grass_mowing' */

/* Version of feed intake and transformation into milk and waste streams from 5.12.2021*/
Harvest harvest_grass_grazing_ext_livestock(Stand *stand)
{
  Harvest sum={{0,0},{0,0},{0,0},{0,0}};
  Pftgrass *grass;
  Pft *pft;
  int p;
  Stocks bm_tot =  {0.0,0.0};
  Stocks bm_gra =  {0.0,0.0};
  Real sum_lai = 0;
  Real hfrac,n2cn,lsu2area,K;
  Real w_C_in_CP,w_C_in_CHO,w_C_in_L,w_C_in_FA,w_N_in_CP,w_FA_in_milk,w_CP_in_milk,w_CHO_in_milk,nreq_milk,nfrac_milk,cfrac_milk;
  Real cfrac_methane,n2c_urine,fcfrac,fnfrac,cp,dmi_max,fdmi,dmi,ne,nem,dign,ne4milk,n4milk,milk,digc;
  Real c_intake,c_methane,c_feces,c_urine,c_milk,c_co2;
  Real n_intake,n_milk,n_feces,n_urine;
  Real Wone=500.0;
  Output *output;
  output=&stand->cell->output;
  foreachpft(pft,p,&stand->pftlist)
  {
    grass=pft->data;
    bm_tot.carbon += grass->ind.leaf.carbon*pft->nind; 
    bm_tot.nitrogen+= grass->ind.leaf.nitrogen*pft->nind;
    sum_lai += actual_lai_grass(pft);
  }

  n2cn = bm_tot.nitrogen/(bm_tot.nitrogen+bm_tot.carbon);
  lsu2area = stand->cell->ml.nr_of_lsus_ext * 1e3 * 1e-4;

  w_C_in_CP  = 0.53; // mass fraction of carbon in crude protein (estimated from chemical composition)
  w_C_in_CHO = 0.44; // mass fraction of carbon in carbohydrates (estimated from chemical composition)
  w_C_in_L   = 0.66; // mass fraction of carbon in lignin (estimated from chemical composition)
  w_C_in_FA  = 0.78; // mass fraction of carbon in fat (estimated from chemical composition)
  w_N_in_CP  = 0.16; // mass fraction of nitrogen in crude protein (Santos and Huber, 2002)
  w_FA_in_milk = 0.04;
  w_CP_in_milk = 0.032;
  w_CHO_in_milk = 0.0485;
  nreq_milk  = w_CP_in_milk / 0.67 * w_N_in_CP; // N requirement for milk kg kg-1
  nfrac_milk = w_CP_in_milk * w_N_in_CP; // mass fraction of nitrogen in milk kg kg-1
  cfrac_milk = w_FA_in_milk * w_C_in_FA + w_CP_in_milk * w_C_in_CP + w_CHO_in_milk * w_C_in_CHO; // mass fraction of carbon in milk kg kg-1
  cfrac_methane = 12.0 / 16.0; // mass fraction of carbon in methane kg kg-1
  n2c_urine = 1; // ratio of N to C in urine; urea=2/1, purine=5/4, creatine=3/4, hippuric_acid=9/1
  K = 0.229 * pow(Wone,0.36); // LAI forh half maximum dmi  

  /* carbon fraction */
  fcfrac = 0.424;
  /* net energy */
  ne  = 0.921 + 8.433 * n2cn;
  /* Fraction of digestible nitrogen */
  dign = 0.914 - 0.494 * exp(-59.559 * n2cn);
  /* Fraction of digestible carbon */
  digc = 0.561 + 2.190 * n2cn;
  
   /* carbon fraction */
  fnfrac = fcfrac * n2cn / (1 - n2cn);
  /* crude protein */
  cp = fnfrac / w_N_in_CP;
  
  /* dry matter intake */
  dmi_max = max(Wone * 1.33 * (2.35 - 3.85 * exp(-0.32 * cp * 100)) / 100, 0);
  /* fraction that can be fulfilled */
  fdmi = pow(sum_lai/K,3)/(1+pow(sum_lai/K,3));
  dmi = dmi_max * fdmi;
  c_intake = dmi * fcfrac;
  n_intake = dmi * fnfrac;

  /* fecal N and C with account for metabolic fecal protein */
  n_feces = n_intake * (1 - dign) + min(0.03 * dmi * w_N_in_CP, n_intake * dign);
  c_feces = c_intake * (1 - digc) + min(0.03 * dmi * w_C_in_CP, c_intake * digc);

  /* add dialy net energy required for maintenance in Mcal day-1 to deficit */
  stand->cell->ml.deficit_lsu_ne += 0.08 * pow(Wone,0.75);
  /* constrain deficit to 1 year of maintenance requirements */
  stand->cell->ml.deficit_lsu_ne = min(stand->cell->ml.deficit_lsu_ne, 365 * 0.08 * pow(Wone,0.75));
  /* fulfill defictit */
  if(dmi * ne > stand->cell->ml.deficit_lsu_ne)
  {
    ne4milk = dmi * ne - stand->cell->ml.deficit_lsu_ne;
    stand->cell->ml.deficit_lsu_ne = 0;
  }
  else
  {
    ne4milk = 0;
    stand->cell->ml.deficit_lsu_ne -= dmi * ne;
  }

  /* add metabolized protein required for maintenance in kgN day-1 to deficit */
  stand->cell->ml.deficit_lsu_mp += 0.0041 * pow(Wone,0.5) * w_N_in_CP;
  /* constrain deficit to 1 year of maintenance requirements */
  stand->cell->ml.deficit_lsu_mp = min(stand->cell->ml.deficit_lsu_mp, 365 * 0.0041 * pow(Wone,0.5) * w_N_in_CP);
  /* fulfill defictit */
  if(n_intake - n_feces > stand->cell->ml.deficit_lsu_mp)
  {
    n4milk = n_intake - n_feces - stand->cell->ml.deficit_lsu_mp;
    stand->cell->ml.deficit_lsu_mp = 0;
  }
  else
  {
    n4milk = 0;
    stand->cell->ml.deficit_lsu_mp -= n_intake - n_feces;
  }

  /* milk production is minimum of milk from N and milk from NE */
  milk = min(n4milk / nreq_milk, ne4milk / (0.360 + 9.69 * w_FA_in_milk));
  n_milk = milk * nfrac_milk;
  c_milk = milk * cfrac_milk;

  /* all N for maintenance and all unused N is excreted */
  n_urine = n_intake - n_feces - n_milk;
  c_urine = n_urine/n2c_urine;
  
  c_methane = dmi * 18.4 * 0.065 / 55.6 * cfrac_methane;
  /* the remainder is respired */
  c_co2 = c_intake - c_feces - c_milk - c_urine - c_methane;
  /* end of calculation of uptake and transformation for 1 LSU */
  
  /* conversion from kg/LSU to g/m2 */
  c_intake *= lsu2area; 
  c_methane *= lsu2area;
  c_feces *= lsu2area;
  c_milk *= lsu2area;
  c_co2 *= lsu2area;
  c_urine *= lsu2area;
  n_intake *= lsu2area;
  n_milk  *= lsu2area;
  n_feces *= lsu2area;
  n_urine *= lsu2area;

/******************************************/

  if(c_intake>epsilon)
  {
    hfrac = c_intake / bm_tot.carbon; 
    foreachpft(pft,p,&stand->pftlist)
    {
      grass=pft->data;
      grass->ind.leaf.carbon *= (1 - hfrac);
      grass->ind.leaf.nitrogen *= (1 - hfrac);

      stand->soil.litter.item[pft->litter].bg.carbon+=grass->ind.root.carbon*hfrac*param.rootreduction*pft->nind;
      output->alittfall.carbon+=grass->ind.root.carbon*hfrac*param.rootreduction*pft->nind*stand->frac;
      stand->soil.litter.item[pft->litter].bg.nitrogen+=grass->ind.root.nitrogen*hfrac*param.rootreduction*pft->nind*pft->par->fn_turnover;
      output->alittfall.nitrogen+=grass->ind.root.nitrogen*hfrac*param.rootreduction*pft->nind*stand->frac*pft->par->fn_turnover;
      output->alitfalln_mgrass+=grass->ind.root.nitrogen*hfrac*param.rootreduction*pft->nind*stand->frac*pft->par->fn_turnover;
      pft->bm_inc.nitrogen+=grass->ind.root.nitrogen*hfrac*param.rootreduction*pft->nind*pft->nind*(1-pft->par->fn_turnover);
      
      grass->ind.root.carbon*=(1-hfrac*param.rootreduction);
      grass->ind.root.nitrogen*=(1-hfrac*param.rootreduction);
      // grass leaves after grazing
      bm_gra.carbon += grass->ind.leaf.carbon*pft->nind; 
      bm_gra.nitrogen+= grass->ind.leaf.nitrogen*pft->nind;
    }

    /* checking sums */
    if (fabs(bm_tot.carbon-bm_gra.carbon-c_intake)>1e-6)
    {
      fprintf(stdout,"grassland harvest C check, c_leaf before: %.4f, after: %.4f, c_intake: %.4f, balance: %.4f \n", bm_tot.carbon,bm_gra.carbon,c_intake,bm_tot.carbon-bm_gra.carbon - c_intake);
      fprintf(stdout,"grassland harvest C components, c_in: %.4f, c_methane: %.4f, c_co2: %.4f, c_feces: %.4f, c_urine: %.4f c_milk: %.4f \n", c_intake,c_methane,c_co2,c_feces,c_urine,c_milk);
      fprintf(stdout,"grassland harvest N components, n_in: %.4f, n_feces: %.4f, n_urine: %.4f n_milk: %.4f \n", n_intake,n_feces,n_urine,n_milk);
      fprintf(stdout,"grassland harvest sums, c_in: %.f, c_sum: %.4f, n_in: %.4f, n_sum: %.4f \n", c_intake,  c_intake+c_methane+c_co2+c_feces+c_urine+c_milk,n_intake,n_feces+n_urine+n_milk);
    }

    sum.harvest.carbon   += c_methane+c_milk+c_co2;       // sum of losses
    stand->soil.pool->fast.carbon += c_feces + c_urine;   // manure back to soil
    sum.harvest.nitrogen += n_milk;                       // sum of losses
    stand->soil.pool->fast.nitrogen += n_feces;             // back to soil
    stand->soil.NH4[0] += n_urine;             // urine to ammonia

    /* writing new mgrass outputs */
    output->uptakec_mgrass += c_intake;
    output->uptaken_mgrass += n_intake;
    output->deficit_lsu_mp = max(output->deficit_lsu_mp, stand->cell->ml.deficit_lsu_mp);
    output->deficit_lsu_ne = max(output->deficit_lsu_ne, stand->cell->ml.deficit_lsu_ne);
    output->yieldc_mgrass += c_milk;
    output->yieldn_mgrass += n_milk;
    output->manurec_mgrass += c_feces;
    output->manuren_mgrass += n_feces;
    output->urinec_mgrass += c_urine;
    output->urinen_mgrass += n_urine;
    output->respc_mgrass += c_co2;
    output->methanec_mgrass += c_methane;
  }
  
  return sum;
} /* of 'harvest_grass_grazing_ext_livestock' */

Harvest harvest_grass_grazing_ext(Stand *stand)
{
  Harvest sum={{0,0},{0,0},{0,0},{0,0}};
  Pftgrass *grass;
  Pft *pft;
  int p;
  Real bm_grazed;
  Real fact;
  Stocks bm_tot =  {0.0,0.0};
  Stocks bm_grazed_pft;
  Real fpc_sum=0.0;
  Real hfrac;
  Output *output;
  output=&stand->cell->output;
  foreachpft(pft,p,&stand->pftlist)
  {
    grass=pft->data;
    bm_tot.carbon += grass->ind.leaf.carbon;
    bm_tot.nitrogen+= grass->ind.leaf.nitrogen;
    fpc_sum+=pft->fpc;
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
    if (grass->ind.leaf.carbon - bm_grazed_pft.carbon < GRAZING_STUBBLE*pft->fpc/fpc_sum)
      bm_grazed_pft.carbon = grass->ind.leaf.carbon - GRAZING_STUBBLE*pft->fpc/fpc_sum;
    if (bm_grazed_pft.carbon < 0)
      bm_grazed_pft.carbon = 0;
    hfrac=bm_grazed_pft.carbon/grass->ind.leaf.carbon;
    pft->gdd = (1-hfrac) * pft->gdd;


    /* Nitrogen */
    //bm_grazed_pft.nitrogen = bm_grazed * fact;
    bm_grazed_pft.nitrogen = hfrac*grass->ind.leaf.nitrogen;

    grass->ind.leaf.carbon -= bm_grazed_pft.carbon;
    sum.harvest.carbon     += (1-MANURE)*bm_grazed_pft.carbon*pft->nind;                       // 60% atmosphere, 15% cows
    stand->soil.pool->fast.carbon += MANURE * bm_grazed_pft.carbon*pft->nind;             // 25% back to soil

    grass->ind.leaf.nitrogen -=  bm_grazed_pft.nitrogen;
    sum.harvest.nitrogen     += (1-MANURE)*bm_grazed_pft.nitrogen*pft->nind;                       // 60% atmosphere, 15% cows
    stand->soil.pool->fast.nitrogen += MANURE * bm_grazed_pft.nitrogen*pft->nind;             // 25% back to soil

    stand->soil.litter.item[pft->litter].bg.carbon+=grass->ind.root.carbon*hfrac*param.rootreduction*pft->nind;
    output->alittfall.carbon+=grass->ind.root.carbon*hfrac*param.rootreduction*pft->nind*stand->frac;
    stand->soil.litter.item[pft->litter].bg.nitrogen+=grass->ind.root.nitrogen*hfrac*param.rootreduction*pft->nind*pft->par->fn_turnover;
    output->alittfall.nitrogen+=grass->ind.root.nitrogen*hfrac*param.rootreduction*pft->nind*stand->frac*pft->par->fn_turnover;
    output->alitfalln_mgrass+=grass->ind.root.nitrogen*hfrac*param.rootreduction*pft->nind*stand->frac*pft->par->fn_turnover;
    pft->bm_inc.nitrogen+=grass->ind.root.nitrogen*hfrac*param.rootreduction*pft->nind*pft->nind*(1-pft->par->fn_turnover);

    grass->ind.root.carbon*=(1-hfrac*param.rootreduction);
    grass->ind.root.nitrogen*=(1-hfrac*param.rootreduction);
    // pft->phen recalculated in phenology_grass
  }
  return sum;
} /* of 'harvest_grass_grazing_ext' */

Harvest harvest_grass_grazing_int(Stand *stand)
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
  Real fpc_sum=0.0;
  Real hfrac;
  Output *output;
  output=&stand->cell->output;
  rotation = &(stand->cell->ml.rotation);
  foreachpft(pft,p,&stand->pftlist)
  {
    grass=pft->data;
    bm_tot.carbon += grass->ind.leaf.carbon;
    bm_tot.nitrogen += grass->ind.leaf.nitrogen;
    fpc_sum+=pft->fpc;
  }

  if (rotation->rotation_mode == RM_UNDEFINED) //initial calculate grazing days and recovery days
  {
    rotation_len = (bm_tot.carbon - GRAZING_STUBBLE) / (1e-4*stand->cell->ml.nr_of_lsus_int * DEMAND_COW_INT) ;
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
    bm_grazed = 1e-4*stand->cell->ml.nr_of_lsus_int * DEMAND_COW_INT * rotation->paddocks;
    foreachpft(pft,p,&stand->pftlist)
    {
      grass=pft->data;
      fact = grass->ind.leaf.carbon / bm_tot.carbon;
      bm_grazed_pft.carbon = bm_grazed * fact;

      if (grass->ind.leaf.carbon - bm_grazed_pft.carbon < GRAZING_STUBBLE*pft->fpc/fpc_sum)
        bm_grazed_pft.carbon = grass->ind.leaf.carbon - GRAZING_STUBBLE*pft->fpc/fpc_sum;

      if (bm_grazed_pft.carbon < 0)
        bm_grazed_pft.carbon =0;

      hfrac=bm_grazed_pft.carbon/grass->ind.leaf.carbon;
      pft->gdd = (1-hfrac) * pft->gdd;

      /* Nitrogen */
      bm_grazed_pft.nitrogen = hfrac*grass->ind.leaf.nitrogen;

      grass->ind.leaf.carbon -= bm_grazed_pft.carbon;
      sum.harvest.carbon     += (1-MANURE)*bm_grazed_pft.carbon*pft->nind;              // 60% atmosphere, 15% cows
      stand->soil.pool->fast.carbon += MANURE * bm_grazed_pft.carbon*pft->nind;    // 25% back to soil

      grass->ind.leaf.nitrogen -= bm_grazed_pft.nitrogen;
      sum.harvest.nitrogen     += (1-MANURE)*bm_grazed_pft.nitrogen*pft->nind;              // 60% atmosphere, 15% cows
      stand->soil.pool->fast.nitrogen += MANURE * bm_grazed_pft.nitrogen*pft->nind;    // 25% back to soil

      stand->soil.litter.item[pft->litter].bg.carbon+=grass->ind.root.carbon*hfrac*param.rootreduction*pft->nind;
      output->alittfall.carbon+=grass->ind.root.carbon*hfrac*param.rootreduction*pft->nind*stand->frac;
      stand->soil.litter.item[pft->litter].bg.nitrogen+=grass->ind.root.nitrogen*hfrac*param.rootreduction*pft->nind*pft->par->fn_turnover;
      output->alittfall.nitrogen+=grass->ind.root.nitrogen*hfrac*param.rootreduction*pft->nind*stand->frac*pft->par->fn_turnover;
      output->alitfalln_mgrass+=grass->ind.root.nitrogen*hfrac*param.rootreduction*pft->nind*stand->frac*pft->par->fn_turnover;
      pft->bm_inc.nitrogen+=grass->ind.root.nitrogen*hfrac*param.rootreduction*pft->nind*pft->nind*(1-pft->par->fn_turnover);

      grass->ind.root.carbon*=(1-hfrac*param.rootreduction);
      grass->ind.root.nitrogen*=(1-hfrac*param.rootreduction);
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
        if (param.newlivestock<1){
          harvest=harvest_grass_grazing_ext(stand);
          //fprintf(stdout,"old grazing harvest routine param: %i\n",param.newlivestock);
        } else {
          harvest=harvest_grass_grazing_ext_livestock(stand);
          //fprintf(stdout,"new grazing harvest routine param: %i\n",param.newlivestock);
        }
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
  output->aharvestn_mgrass+=(harvest.harvest.nitrogen+harvest.residual.nitrogen+
    +harvest.residuals_burnt.nitrogen+harvest.residuals_burntinfield.nitrogen)*stand->frac;
  output->dcflux+=(harvest.harvest.carbon+harvest.residual.carbon)*stand->frac;
  stand->growing_days=0;
  return harvest;
}
 /* of 'harvest_stand' */
