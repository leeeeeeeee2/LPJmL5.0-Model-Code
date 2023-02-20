/**************************************************************************************/
/**                                                                                \n**/
/**                s  t  a  n  d  s  t  o  c  k  s  .  c                           \n**/
/**                                                                                \n**/
/**     C implementation of LPJmL                                                  \n**/
/**                                                                                \n**/
/**     Function computes total carbon and nitrogen in stand                       \n**/
/**                                                                                \n**/
/** (C) Potsdam Institute for Climate Impact Research (PIK), see COPYRIGHT file    \n**/
/** authors, and contributors see AUTHORS file                                     \n**/
/** This file is part of LPJmL and licensed under GNU AGPL Version 3               \n**/
/** or later. See LICENSE file or go to http://www.gnu.org/licenses/               \n**/
/** Contact: https://github.com/PIK-LPJmL/LPJmL                                    \n**/
/**                                                                                \n**/
/**************************************************************************************/

#include "lpj.h"

Stocks standstocks(const Stand *stand /**< pointer to stand */
                  )                   /** \return stocks sum (gC/m2,gN/m2) */
{
  int p;
  const Pft *pft;
  Stocks tot;
  tot=soilstocks(&stand->soil); /* get stocks in soil */
  foreachpft(pft,p,&stand->pftlist)
  {
    tot.carbon+=vegc_sum(pft); /* sum up carbon in PFTs */
    tot.nitrogen+=vegn_sum(pft); /* sum up nitrogen in PFTs */
  }
  return tot;
} /* of 'standstocks' */
