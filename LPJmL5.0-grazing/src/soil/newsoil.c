/**************************************************************************************/
/**                                                                                \n**/
/**                n  e  w  _  s  o  i  l  .  c                                    \n**/
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

void newsoil(Soil *soil /**< pointer to soil data */)
{
  int l;
  soil->litter.n=0;
  soil->litter.item=NULL;
  soil->litter.agtop_wcap=soil->litter.agtop_moist=soil->litter.agtop_cover=soil->litter.agtop_temp=0;
  forrootsoillayer(l)
    soil->pool[l].fast.carbon=soil->pool[l].slow.carbon=soil->pool[l].fast.nitrogen=soil->pool[l].slow.nitrogen=soil->YEDOMA=0.0;
  soil->snowheight=soil->snowfraction=0;
} /* of 'newsoil' */
