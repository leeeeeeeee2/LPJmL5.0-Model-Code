/**************************************************************************************/
/**                                                                                \n**/
/**                 c  l  i  m  b  u  f  .  h                                      \n**/
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

#ifndef CLIMBUF_H /* Already included? */
#define CLIMBUF_H

#include "date.h"

#define CLIMBUFSIZE 20 /* size of climate buffer */
#define NDAYS 31       /* maximum number of days in a month */
#define NDAYS_PREC 10

/* Definitions of datatypes */

typedef struct
{
  Real gdd5; /* number of days with temp > 5 deg C */
  Real temp[NDAYS];
  Real dval_prec[NDAYS+1]; /* daily precipitation values (mm) */
  Real temp_min; /* minimum annual temperature (deg C) */
  Real temp_max; /* maximum annual temperature (deg C) */
  Real atemp_mean;
  Real atemp; /* annual mean temperature, used to update atemp_mean20 */
  Real atemp_mean20; /* 20-year mean of annual mean temperature */
  Real atemp_mean20_fix; /* 20-year mean of annual mean temperature, fixed for phen_variety computation*/
  Real mtemp;
  Real mprec;
  Real mpet;
  Buffer min,max;
  Real mprec20[NMONTH]; /* 20-year average monthly precip */
  Real mpet20[NMONTH]; /* 20-year average monthly precip */
  Real mtemp20[NMONTH]; /* 20-year average monthly precip */
  Real mtemp_min20; /*20-year average of coldest month temperature*/
  Real aetp_mean;       /*20-year average of evapotranspiration of the cell*/
  Real aprec;     /**< annual sum of 20-year average monthly precip */
  Real *V_req;            /* 20-year average of vernalization requirements in days */
  Real *V_req_a;          /* annual vernalization requirements in days, used to calculate V_req */

} Climbuf;

/* Definitions of macros */

#define initgdd5(climbuf) climbuf.gdd5=0
#define updategdd5(climbuf,temp) if(temp>5) (*climbuf).gdd5++

/* Declaration of functions */

extern Bool new_climbuf(Climbuf *,int);
extern void init_climbuf(Climbuf *,int);
extern void daily_climbuf(Climbuf *,Real);
extern void getmintemp20_n(const Climbuf *,Real [],int);
extern void monthly_climbuf(Climbuf *,Real,Real,Real,int);
extern void annual_climbuf(Climbuf *,Real,int,int,int,int);
extern Bool fwriteclimbuf(FILE *,const Climbuf *,int);
extern Bool freadclimbuf(FILE *,Climbuf *,Bool,int);
extern void freeclimbuf(Climbuf *);
extern void annual_climbuf(Climbuf *,Real,int,int,int,int); 
#endif
