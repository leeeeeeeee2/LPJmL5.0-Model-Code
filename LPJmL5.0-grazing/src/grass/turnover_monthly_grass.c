/**************************************************************************************/
/**                                                                                \n**/
/**        t u r n o v e r _ m o n t h l y _ g r a s s . c                         \n**/
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

void turnover_monthly_grass(Litter *litter,Pft *pft)
{
  Pftgrass *grass;
  Pftgrasspar *grasspar;
  Output *output;
  if(pft->stand->type->landusetype==NATURAL)
  {
    grasspar=pft->par->data;
    grass=pft->data;
    output=&pft->stand->cell->output;
    grass->turn.root.carbon+=grass->ind.root.carbon*grasspar->turnover.root/NMONTH;
    grass->turn.root.nitrogen+=grass->ind.root.nitrogen*grasspar->turnover.root/NMONTH;
    litter->item[pft->litter].bg.carbon+=grass->ind.root.carbon*grasspar->turnover.root/NMONTH*pft->nind;
    output->alittfall.carbon+=grass->ind.root.carbon*grasspar->turnover.root/NMONTH*pft->nind*pft->stand->frac;
    litter->item[pft->litter].bg.nitrogen+=grass->ind.root.nitrogen*grasspar->turnover.root/NMONTH*pft->nind*pft->par->fn_turnover;
    output->alittfall.nitrogen+=grass->ind.root.nitrogen*grasspar->turnover.root/NMONTH*pft->nind*pft->par->fn_turnover*pft->stand->frac;
  }
} /* of 'turnover_monthly_grass' */
