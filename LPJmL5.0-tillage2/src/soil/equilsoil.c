/**************************************************************************************/
/**                                                                                \n**/
/**                     e  q  u  i  l  s  o  i  l  . c                             \n**/
/**                                                                                \n**/
/**     C implementation of LPJmL                                                  \n**/
/**                                                                                \n**/
/**     SOIL DECOMPOSITION EQUILIBRIUM CALCULATION                                 \n**/
/**                                                                                \n**/
/**     Analytical solution of differential flux equations for fast and            \n**/
/**     slow soil carbon pools.  Implemented after (soil_equil_year)               \n**/
/**     simulation years, when annual litter inputs should be close to             \n**/
/**     equilibrium.  Assumes average climate (temperature and soil                \n**/
/**     moisture) from all years up to soil_equil_year.                            \n**/
/**                                                                                \n**/
/**       Rate of change of soil pool size = litter input - decomposition          \n**/
/**         (5) dc/dt = litter_decom - kc                                          \n**/
/**       At equilibrium,                                                          \n**/
/**         (6) dc/dt = 0                                                          \n**/
/**       From (5) & (6),                                                          \n**/
/**         (7) c = litter_decom / k                                               \n**/
/**                                                                                \n**/
/** (C) Potsdam Institute for Climate Impact Research (PIK), see COPYRIGHT file    \n**/
/** authors, and contributors see AUTHORS file                                     \n**/
/** This file is part of LPJmL and licensed under GNU AGPL Version 3               \n**/
/** or later. See LICENSE file or go to http://www.gnu.org/licenses/               \n**/
/** Contact: https://github.com/PIK-LPJmL/LPJmL                                    \n**/
/**                                                                                \n**/
/**************************************************************************************/

#include "lpj.h"

void equilsoil(Soil *soil,            /**< pointer to soil data */
               int ntotpft,           /**< total number of PFTs */
               const Pftpar pftpar[] /**< PFT parameter array */
              )                       /** \return void         */
{
  int l,p,f;
  Real sumlitter,pftlitter,wood=0,socfraction;
  Poolpar *k_mean,*c0,*sum;
  Poolpar cn_ratio[LASTLAYER];
  Pool c_before[LASTLAYER];
  Real total1=0, total2=0;
  k_mean=newvec(Poolpar,ntotpft);
  check(k_mean);
  c0=newvec(Poolpar,ntotpft);
  check(c0);
  sum=newvec(Poolpar,ntotpft);
  check(sum);
  sumlitter=littercarbon(&soil->litter);

  soil->decomp_litter_mean.carbon/=(soil_equil_year-veg_equil_year+1);
  soil->decomp_litter_mean.nitrogen/=(soil_equil_year-veg_equil_year+1);
  forrootsoillayer(l)
  {
    soil->k_mean[l].fast/=(soil_equil_year-veg_equil_year+1);
    soil->k_mean[l].slow/=(soil_equil_year-veg_equil_year+1);
    c_before[l].fast.carbon=soil->pool[l].fast.carbon;
    c_before[l].slow.carbon=soil->pool[l].slow.carbon;
  }
  for(p=0;p<ntotpft;p++)
    k_mean[p].fast=k_mean[p].slow=sum[p].fast=sum[p].slow=0.0;
  forrootsoillayer(l)
    for(p=0;p<ntotpft;p++)
    {
      socfraction=pow(10,pftpar[p].soc_k*logmidlayer[l])
                 - (l>0 ? pow(10,pftpar[p].soc_k*logmidlayer[l-1]): 0);
      k_mean[p].fast+=soil->k_mean[l].fast*socfraction;
      k_mean[p].slow+=soil->k_mean[l].slow*socfraction;
    }
  for(p=0;p<ntotpft;p++)
  {
    c0[p].fast=k_mean[p].fast>epsilon ? (1-param.atmfrac)*param.fastfrac*(soil->decomp_litter_mean.carbon)/k_mean[p].fast : 0;
    c0[p].slow=k_mean[p].slow>epsilon ? (1-param.atmfrac)*(1.0-param.fastfrac)*(soil->decomp_litter_mean.carbon)/k_mean[p].slow : 0;
  }

  if(soil->decomp_litter_mean.carbon>epsilon)
  {
    forrootsoillayer(l)
    {
      cn_ratio[l].slow=soil->pool[l].slow.carbon/soil->pool[l].slow.nitrogen;
      cn_ratio[l].fast=soil->pool[l].fast.carbon/soil->pool[l].fast.nitrogen;
      soil->pool[l].slow.carbon=soil->pool[l].fast.carbon=0;
      soil->pool[l].slow.nitrogen=soil->pool[l].fast.nitrogen=0;
    }
    forrootsoillayer(l)
    {
      for(p=0;p<soil->litter.n;p++)
      {
        wood=0;
        for(f=0;f<NFUELCLASS;f++)
          wood+=soil->litter.item[p].ag.wood[f].carbon+soil->litter.item[p].agsub.wood[f].carbon;
        pftlitter=soil->litter.item[p].bg.carbon+soil->litter.item[p].ag.leaf.carbon+soil->litter.item[p].agsub.leaf.carbon+wood;
        socfraction=pow(10,soil->litter.item[p].pft->soc_k*logmidlayer[l])
                    - (l>0 ? pow(10,soil->litter.item[p].pft->soc_k*logmidlayer[l-1]): 0);
            if(sumlitter>0)
              soil->pool[l].slow.carbon+=c0[soil->litter.item[p].pft->id].slow*socfraction*pftlitter/sumlitter;
            if(sumlitter>0)
              soil->pool[l].fast.carbon+=c0[soil->litter.item[p].pft->id].fast*socfraction*pftlitter/sumlitter;
      }
      soil->pool[l].slow.nitrogen=soil->pool[l].slow.carbon/soil->par->cn_ratio;
      soil->pool[l].fast.nitrogen=soil->pool[l].fast.carbon/soil->par->cn_ratio;
      //soil->pool[l].slow.nitrogen=soil->pool[l].slow.carbon/cn_ratio[l].slow;
      //soil->pool[l].fast.nitrogen=soil->pool[l].fast.carbon/cn_ratio[l].fast;
      //soil->NH4[l]=0;//soil->pool[l].slow.nitrogen/10000;
      //soil->NO3[l]=0;//soil->pool[l].slow.nitrogen/10000;
      //if (soil->decomp_litter_mean.nitrogen<0.2)
      if (soil->decomp_litter_mean.nitrogen<epsilon)
      {
        soil->NH4[l]=0;
        soil->NO3[l]=0;
        soil->pool[l].slow.nitrogen=0;
        soil->pool[l].fast.nitrogen=0;
      }
    }
    forrootsoillayer(l)
    {
      for(p=0;p<ntotpft;p++)
      {
        socfraction=pow(10,pftpar[p].soc_k*logmidlayer[l])
                  - (l>0 ? pow(10,pftpar[p].soc_k*logmidlayer[l-1]): 0);
        soil->c_shift_fast[l][p]=socfraction*soil->k_mean[l].fast/k_mean[p].fast;
        soil->c_shift_slow[l][p]=socfraction*soil->k_mean[l].slow/k_mean[p].slow;
        sum[p].fast+=soil->c_shift_fast[l][p];
        sum[p].slow+=soil->c_shift_slow[l][p];
      }

      soil->k_mean[l].slow=soil->k_mean[l].fast=0.0;
    }
    for(p=0;p<ntotpft;p++)
    {
      if(sum[p].fast<1.0 && sum[p].fast>0.0)
      //if(sum[p].fast>=epsilon)
        for (l=0;l<LASTLAYER;l++) soil->c_shift_fast[l][p]=soil->c_shift_fast[l][p]/sum[p].fast;
      else if (sum[p].fast<epsilon)
      {
        soil->c_shift_fast[0][p]=1.0;
        for (l=1;l<LASTLAYER;l++) soil->c_shift_fast[l][p]=0;
      }
      if(sum[p].slow<1.0 && sum[p].slow>0.0)
      //if(sum[p].slow>=epsilon)
        for (l=0;l<LASTLAYER;l++) soil->c_shift_slow[l][p]=soil->c_shift_slow[l][p]/sum[p].slow;
      else if (sum[p].slow<epsilon)
      {
        soil->c_shift_slow[0][p]=1.0;
        for (l=1;l<LASTLAYER;l++) soil->c_shift_slow[l][p]=0;
      }
    }
    soil->count=0;
    soil->decomp_litter_mean.carbon=soil->decomp_litter_mean.nitrogen=0.0;
  }
  else
  { /* if there is no carbon, remove mineral N from initialization */
    forrootsoillayer(l)
    {
      soil->NH4[l]=soil->NO3[l]=0.0;
      soil->pool[l].slow.nitrogen=0;
      soil->pool[l].fast.nitrogen=0;
    }
  }
  /* avoiding soilC increases in regions with low thawdepth (=permafrost) */
  if (soil->mean_maxthaw<1500)
  {
    forrootsoillayer(l)
    {
      total1+=soil->pool[l].fast.carbon+soil->pool[l].slow.carbon;
      total2+=c_before[l].fast.carbon+c_before[l].slow.carbon;
    }
    /* if soilcarbon is to increase by 10% in permafrost regions, reset to original */
    if (total1>total2*1.1)
    {
      forrootsoillayer(l)
      {
        soil->pool[l].fast.carbon=c_before[l].fast.carbon;
        soil->pool[l].slow.carbon=c_before[l].slow.carbon;
        soil->pool[l].slow.nitrogen=soil->pool[l].slow.carbon/cn_ratio[l].slow;
        soil->pool[l].fast.nitrogen=soil->pool[l].fast.carbon/cn_ratio[l].fast;

      }
    }
  }
  /* freeing memory */
  free(k_mean);
  free(c0);
  free(sum);
} /* of 'equilsoil' */
