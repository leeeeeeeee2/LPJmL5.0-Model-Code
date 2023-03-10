/**************************************************************************************/
/**                                                                                \n**/
/**    p  r  o  d  u  c  t  _  t  u  r  n  o  v  e  r  .  c                        \n**/
/**                                                                                \n**/
/**     extension of LPJ to compute turnover of product pools                      \n**/
/**                                                                                \n**/
/** (C) Potsdam Institute for Climate Impact Research (PIK), see COPYRIGHT file    \n**/
/** authors, and contributors see AUTHORS file                                     \n**/
/** This file is part of LPJmL and licensed under GNU AGPL Version 3               \n**/
/** or later. See LICENSE file or go to http://www.gnu.org/licenses/               \n**/
/** Contact: https://github.com/PIK-LPJmL/LPJmL                                    \n**/
/**                                                                                \n**/
/**************************************************************************************/


#include "lpj.h"

#ifdef IMAGE

Real product_turnover(Image_data *data)
{
  data->prod_turn_fast=(data->timber.fast)*PRODUCT_FAST;
  data->prod_turn_slow=(data->timber.slow)*PRODUCT_SLOW;
  
  data->timber.fast*=(1.0-PRODUCT_FAST);
  data->timber.slow*=(1.0-PRODUCT_SLOW);
  return data->prod_turn_fast+data->prod_turn_slow;
} /* of 'product_turnover' */

#else
void product_turnover(Cell *cell)
{
  cell->output.prod_turnover.carbon+=cell->ml.product.fast.carbon*PRODUCT_FAST;
  cell->output.prod_turnover.carbon+=cell->ml.product.slow.carbon*PRODUCT_SLOW;
  cell->ml.product.fast.carbon*=(1.0-PRODUCT_FAST);
  cell->ml.product.slow.carbon*=(1.0-PRODUCT_SLOW);
  cell->output.prod_turnover.nitrogen+=cell->ml.product.fast.nitrogen*PRODUCT_FAST;
  cell->output.prod_turnover.nitrogen+=cell->ml.product.slow.nitrogen*PRODUCT_SLOW;
  cell->ml.product.fast.nitrogen*=(1.0-PRODUCT_FAST);
  cell->ml.product.slow.nitrogen*=(1.0-PRODUCT_SLOW);
}
#endif
