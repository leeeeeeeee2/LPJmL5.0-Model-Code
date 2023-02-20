/**************************************************************************************/
/**                                                                                \n**/
/**                l  i  t  t  e  r  _  a  g  _  s  u  m  .  c                     \n**/
/**                                                                                \n**/
/**     C implementation of LPJmL                                                  \n**/
/**                                                                                \n**/
/**     Function computes sum of all above-ground litter carbon pools              \n**/
/**                                                                                \n**/
/** (C) Potsdam Institute for Climate Impact Research (PIK), see COPYRIGHT file    \n**/
/** authors, and contributors see AUTHORS file                                     \n**/
/** This file is part of LPJmL and licensed under GNU AGPL Version 3               \n**/
/** or later. See LICENSE file or go to http://www.gnu.org/licenses/               \n**/
/** Contact: https://github.com/PIK-LPJmL/LPJmL                                    \n**/
/**                                                                                \n**/
/**************************************************************************************/

#include "lpj.h"

Real litter_ag_sum(const Litter *litter /**< pointer to litter data */
                  )                     /** \return aboveground litter (gC/m2) */
{
  int i,l;
  Real sum;
  sum=0;
  for(l=0;l<litter->n;l++)
  {
    sum+=litter->item[l].ag.leaf.carbon;
    for(i=0;i<NFUELCLASS;i++)
      sum+=litter->item[l].ag.wood[i].carbon;
  }
  return sum;
} /* of litter_ag_sum */

Real litter_ag_sum_n(const Litter *litter)
{
  int i,l;
  Real sum;
  sum=0;
  for(l=0;l<litter->n;l++)
  {
    sum+=litter->item[l].ag.leaf.nitrogen;
    for(i=0;i<NFUELCLASS;i++)
      sum+=litter->item[l].ag.wood[i].nitrogen;
  }
  return sum;
} /* of litter_ag_sum_n */

Real litter_agsub_sum(const Litter *litter /**< pointer to litter data */
                  )                     /** \return aboveground litter (gC/m2) */
{
  int i,l;
  Real sum;
  sum=0;
  for(l=0;l<litter->n;l++)
  {
    sum+=litter->item[l].agsub.leaf.carbon;
    for(i=0;i<NFUELCLASS;i++)
      sum+=litter->item[l].agsub.wood[i].carbon;
  }
  return sum;
} /* of litter_agsub_sum */

Real litter_agsub_sum_n(const Litter *litter)
{
  int i,l;
  Real sum;
  sum=0;
  for(l=0;l<litter->n;l++)
  {
    sum+=litter->item[l].agsub.leaf.nitrogen;
    for(i=0;i<NFUELCLASS;i++)
      sum+=litter->item[l].agsub.wood[i].nitrogen;
  }
  return sum;
} /* of litter_agsub_sum_n */

