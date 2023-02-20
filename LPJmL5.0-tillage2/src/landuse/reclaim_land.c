/**************************************************************************************/
/**                                                                                \n**/
/**          r  e  c  l  a  i  m  _  l a  n  d  .  c                               \n**/
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
#include "grass.h"
#include "tree.h"

static void remove_vegetation_copy(Soil *soil, /* soil pointer */
                                   const Stand *stand, /* stand pointer */
                                   Cell *cell, /* cell pointer */
                                   Real standfrac, /* stand fraction (0..1) */
                                   Bool istimber /* IMAGE coupling (TRUE/FALSE) */
                                  )
{
  int p;
  Pft *pft;
  Real nind;
  Real ftimber; /* fraction harvested for timber */
  Stocks harvest;
  Stocks stocks;
#ifdef IMAGE
  Bool tharvest=FALSE;
  ftimber=min(1,cell->ml.image_data->timber_frac/standfrac);
#else
  Poolpar frac1,frac2;
  ftimber=1;
  frac1.fast=frac2.fast=1.0;
  frac1.slow=frac2.slow=0.0;
#endif

  foreachpft(pft,p,&stand->pftlist)
  {
    nind = pft->nind;
    /* if plot is deforested, wood is returned to litter, harvested or burnt
    * allows for mixed use, first harvesting a fraction of the stand,
    * then burning a fraction, then returning the rest to the litter pools
    */
    if(istimber)
    {
      if(pft->par->type==TREE)
      {
#ifdef DEBUG_IMAGE
        if/*(ftimber>0 ||
          (cell->coord.lon-.1<-43.25 && cell->coord.lon+.1>-43.25 && cell->coord.lat-.1<-11.75 && cell->coord.lat+.1>-11.75)||
          (cell->coord.lon-.1<94.25 && cell->coord.lon+.1>94.25 && cell->coord.lat-.1<22.25 && cell->coord.lat+.1>22.25))*/
          (cell->coord.lon>102.2 && cell->coord.lon < 102.3 && cell->coord.lat >28.7 && cell->coord.lat< 28.8)
        {
          printf("A %g/%g timber_harvest: %s ftimber %g standfrac %g littersum %g vegcsum %g nind %g %g\n",cell->coord.lon,cell->coord.lat,pft->par->name,ftimber,
                 standfrac,littersum(&soil->litter),vegc_sum(pft),pft->nind,nind);
          fflush(stdout);
        }
#endif
        /* harvesting timber */
        cell->output.ftimber=ftimber;
#ifdef IMAGE
        tharvest=TRUE;
        harvest=timber_harvest(pft,soil,&cell->ml.image_data->timber,
                               cell->ml.image_data->timber_f,ftimber,standfrac,&nind,&cell->output.trad_biofuel);
#ifdef DEBUG_IMAGE
        if(ftimber>0 ||
          (cell->coord.lon-.1<-43.25 && cell->coord.lon+.1>-43.25 && cell->coord.lat-.1<-11.75 && cell->coord.lat+.1>-11.75)||
          (cell->coord.lon-.1<94.25 && cell->coord.lon+.1>94.25 && cell->coord.lat-.1<22.25 && cell->coord.lat+.1>22.25))*/
        if(cell->coord.lon>102.2 && cell->coord.lon < 102.3 && cell->coord.lat >28.7 && cell->coord.lat< 28.8)
        {
          printf("B %g/%g timber_burn: %s fburn %g standfrac %g littersum %g vegcsum %g nind %g %g\n",cell->coord.lon,cell->coord.lat,pft->par->name,cell->ml.image_data->fburnt,
                 standfrac,littersum(&soil->litter),vegc_sum(pft),pft->nind,nind);
          fflush(stdout);
        }
#endif
#else
        harvest=timber_harvest(pft,soil,&frac1,frac2,param.ftimber,standfrac,&nind,&cell->output.trad_biofuel);
#endif
        cell->output.timber_harvest.carbon+=harvest.carbon;
        cell->output.timber_harvest.nitrogen+=harvest.nitrogen;
#ifdef IMAGE
        /* burning wood */
        cell->output.fburn=cell->ml.image_data->fburnt;
#ifdef DEBUG_IMAGE
        printf("fburnt %g %g\n",cell->output.fburn,cell->ml.image_data->fburnt);
        fflush(stdout);
#endif
        stocks=timber_burn(pft,cell->ml.image_data->fburnt,&soil->litter,nind);
#else
        stocks=timber_burn(pft,param.fburnt,&soil->litter,nind);
#endif
        cell->output.deforest_emissions.carbon+=stocks.carbon*standfrac;
        cell->output.deforest_emissions.nitrogen+=stocks.nitrogen*(1-param.q_ash)*standfrac;
        soil->NO3[0]+=stocks.nitrogen*param.q_ash*standfrac;
      } /* if tree */
    } /* is timber */
#ifdef DEBUG_IMAGE
    if(cell->coord.lon>102.2 && cell->coord.lon < 102.3 && cell->coord.lat >28.7 && cell->coord.lat< 28.8)
    {
      printf("C %g/%g timber_burn: %s fburn %g standfrac %g littersum %g vegcsum %g nind %g %g\n",cell->coord.lon,cell->coord.lat,pft->par->name,cell->ml.image_data->fburnt,
             standfrac,littersum(&soil->litter),vegc_sum(pft),pft->nind,nind);
      fflush(stdout);
    }
#endif
    /* rest goes to litter */
    litter_update(&soil->litter,pft,nind);
#ifdef DEBUG_IMAGE
    if(ftimber>0 ||
      (cell->coord.lon-.1<-43.25 && cell->coord.lon+.1>-43.25 && cell->coord.lat-.1<-11.75 && cell->coord.lat+.1>-11.75)||
      (cell->coord.lon-.1<94.25 && cell->coord.lon+.1>94.25 && cell->coord.lat-.1<22.25 && cell->coord.lat+.1>22.25))
    {
      printf("D %g/%g timber_burn: %s fburn %g standfrac %g littersum %g vegcsum %g nind %g %g\n",cell->coord.lon,cell->coord.lat,pft->par->name,cell->ml.image_data->fburnt,
             standfrac,littersum(&soil->litter),vegc_sum(pft),pft->nind,nind);
      fflush(stdout);
    }
#endif

  } /* of foreachpft */
#ifdef IMAGE
  if(tharvest)
    cell->ml.image_data->timber_frac-=ftimber;
#endif
  /* There should be no vegetation on stand2, only the soil carbon was copied to stand2
   * therefore delpft is not necessary here */

}/* of 'remove_vegetation_copy' */

void reclaim_land(const Stand *stand1,Stand *stand2,Cell *cell,Bool istimber, int ntotpft)
{
  int l,p;
  Soil *soil;
  soil=&stand2->soil;
  stand2->fire_sum=stand1->fire_sum;

  for (l=0;l<LASTLAYER;l++)
  {
    soil->c_shift_fast[l]=newvec(Real,ntotpft);
    check(soil->c_shift_fast[l]);
    soil->c_shift_slow[l]=newvec(Real,ntotpft);
    check(soil->c_shift_slow[l]);
  }
  for (p=0;p<ntotpft;p++)
  {
    soil->c_shift_fast[0][p]=1;
    soil->c_shift_slow[0][p]=1;
  }
  for (l=1;l<LASTLAYER;l++)
    for (p=0;p<ntotpft;p++)
    {
      soil->c_shift_fast[l][p]=0;
      soil->c_shift_slow[l][p]=0;
    }

  copysoil(&stand2->soil,&stand1->soil,ntotpft);
  /*allow one tillage event on new stand upon cultivation, after deforest or landuseexpansion*/
//  tillage(&stand2->soil, param.residue_frac); //now in landusechange.c
  for(l=0;l<NSOILLAYER;l++)
    stand2->frac_g[l]=stand1->frac_g[l];
  remove_vegetation_copy(&stand2->soil,stand1,cell,stand2->frac,
                         istimber);
}/* of 'reclaim_land' */

/*
- called in cultivate(), deforest(), regrowth(), landexpansion(),
  grasslandreduction()
- copies the values of the variables of one stand to another stand except of the
  pftlist
- calls local function copysoil()
  -> copysoil() copies the values of all variables in struct Soil (see soil.h)
- calls local function remove_vegetation_crop()
  -> remove_vegetation_crop() removes the pfts from the new stand by updating the
     litter pools of the new stand
*/
