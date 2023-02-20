/**************************************************************************************/
/**                                                                                \n**/
/**                         s  o  i  l  .  h                                       \n**/
/**                                                                                \n**/
/**     C implementation of LPJmL                                                  \n**/
/**                                                                                \n**/
/**     Header file for soil functions                                             \n**/
/**                                                                                \n**/
/** (C) Potsdam Institute for Climate Impact Research (PIK), see COPYRIGHT file    \n**/
/** authors, and contributors see AUTHORS file                                     \n**/
/** This file is part of LPJmL and licensed under GNU AGPL Version 3               \n**/
/** or later. See LICENSE file or go to http://www.gnu.org/licenses/               \n**/
/** Contact: https://github.com/PIK-LPJmL/LPJmL                                    \n**/
/**                                                                                \n**/
/**************************************************************************************/

#ifndef SOIL_H /* Already included? */
#define SOIL_H

/* Definition of constants */

#define NSOILLAYER 6 /* Number of soil layers */
#define LASTLAYER (NSOILLAYER-1)
#define FRACGLAYER NSOILLAYER
#define BOTTOMLAYER (NSOILLAYER-1)
#define NFUELCLASS 4 /* Number of fuel classes */
#define TOPLAYER 0
#define NTILLLAYER 1 /* number of layers to be tilled */

#define SNOWLAYER NSOILLAYER

#define veg_equil_year 2400
#define soil_equil_year veg_equil_year+1320
#define cshift_year 120
#define snow_skin_depth 40.0 /* snow skin layer depth (mm water equivalent)*/
#define c_water 4.2e6 /* J/m3/K */
#define c_ice   2.1e6 /* J/m3/K */
#define c_snow  2.1e6 /* snow heat capacity (J/ton/K) */
#define m_heat  4.0e4  /*J/gC microb. heating*/
#define lambda_snow 0.2
#define th_diff_snow (lambda_snow/snowheatcap) /* thermal diffusivity of snow [m2/s]*/
#define snow_density     0.3   /*ton/m3*/
#define snowheatcap     6.3e5  /*c_snow*snow_density [J/m3/K]*/
#define c_soilheatcapdry 1.2e6/*1.2e6   J/m3/K */
#define c_albsnow        0.65 /* Albedo of snow (0-1) */
#define c_albsoil        0.30 /* Albedo of bare soil (0-1). Should be soil and soil moisture dependent */
#define c_albwater       0.10 /* Albedo of water (0-1). */
#define c_snowcanopyfrac  0.40 /* Maximum snow cover fraction in the canopy (with leaves) */
#define c_watertosnow    6.70 /* Conversion factor from water to snowdepth, i.e. 1 cm water equals 6.7 cm of snow. Taken from RAKMO and personal communication with KNMI. */
#define c_mineral 1.9259e6 /* [J/m3/K] after deVries 1963*/
#define c_water2ice 0.3e9 /* the energy that is needed/released during water/ice conversion (J/m3)*/
#define tsnow 0.0
/* assuming lambda_air=0.025 and c_air = 1300; admit_air=sqrt(lambda_air*c_air)*/
#define admit_air 5.7009
#define T_zero 0   /* [deg C] */
#define maxheatsteps 100
#define NSTEP_DAILY 1
#define LAG_CONV (NDAYYEAR*0.5*M_1_PI)  /* conversion factor for oscillation
                                lag from angular units to days (=365/(2*PI))*/
#define MINERALDENS 2700 /* mineral density in kg/m3*/
#define ORGDENS 1400     /* density of organic soil substances [kg/m3]*/
#define PRIESTLEY_TAYLOR 1.32 /* Priestley-Taylor coefficient */
#define SOILDEPTH_IRRIG 500 /*size of layer considered for calculation of irrigation ammount*/
#define CDN 1.2         /* shape factor for denitrification from SWAT; beta_denit eq 3:1.4.1, SWAT Manual 2009 , take smaller value as it seems to be to high 1.4 originally*/

/* Declaration of variables */

extern Real soildepth[NSOILLAYER];
extern Real layerbound[NSOILLAYER];
extern Real midlayer[NSOILLAYER];
extern Real logmidlayer[NSOILLAYER];
extern Real fbd_fac[NFUELCLASS];

#include "soilpar.h"

/* Definition of datatypes */

typedef enum {NOSTATE,BELOW_T_ZERO,AT_T_ZERO,ABOVE_T_ZERO,FREEZING,
              THAWING} Soilstate;

typedef struct
{
  Real fast;       /**< fast-decomposing component */
  Real slow;       /**< slow-decomposing component */
} Poolpar;

typedef struct
{
  Stocks fast;       /**< fast-decomposing component */
  Stocks slow;       /**< slow-decomposing component */
} Pool;

typedef struct
{
  Real leaf;       /**< leaf litter (gC/m2)  */
  Real wood;       /**< woody litter (gC/m2) */
  Real q10_wood;
} Traitpar;

typedef struct
{
  Stocks leaf;             /**< leaf litter (gC/m2, gN/m2)  */
  Stocks wood[NFUELCLASS]; /**< woody litter (gC/m2, gN/m2) */
} Trait;

typedef struct
{
  const struct Pftpar *pft; /**< PFT id for litter */
  Trait ag;                 /**< above-ground leaf and wood litter */
  Trait agsub;              /**< above-ground leaf and wood litter incorporated to first layer through tillage */
  Stocks bg;                /**< below-ground litter (gC/m2, gN/m2) */
} Litteritem;

typedef struct
{
  Real avg_fbd[NFUELCLASS+1]; /**< average fuel bulk densities */
  Litteritem *item;           /**< litter list for PFTs */
  int n;                      /**< number of litter pools */
  Real agtop_wcap;            /**< capacity of ag litter to store water in mm */
  Real agtop_moist;           /**< amount of water stored in ag litter in mm */
  Real agtop_cover;           /**< fraction of soil coverd by ag litter */
  Real agtop_temp;            /**< temperature of ag litter */
#ifdef MICRO_HEATING
  Real decomC;                /**< litter decomposed (gC/m2) */
#endif
} Litter;

typedef struct
{
  Stocks harvest,residual,residuals_burnt,residuals_burntinfield;
} Harvest;

typedef struct
{
  int type;   /**< unique soil type id */
  char *name; /**< soil name */
  Real Ks;    /**< hydraulic conductivity (mm/h)*/
  Real Sf;    /**< Suction head (mm)*/
  Real sand;  /**< fraction of sand content in soil texture*/
  Real silt;  /**< fraction of silt content in soil texture*/
  Real clay;  /**< fraction of clay content in soil texture*/
  int hsg;        /**< hydrological soil group for CN */
  Real tdiff_0;   /**< thermal diffusivity (mm^2/s) at wilting point (0% whc) */
  Real tdiff_15;  /**< thermal diffusivity (mm^2/s) at 15% whc */
  Real tdiff_100; /**< thermal diffusivity (mm^2/s) at field capacity (100% whc) */
  Real tcond_pwp, tcond_100, tcond_100_ice; /**< thermal conductivity [W/m/K]*/
  Real z_nit;
  Real a_nit;
  Real b_nit;
  Real c_nit;
  Real d_nit;
  Real m_nit;
  Real n_nit;
  Real a_denit;
  Real b_denit;
  Real denit_rate;
  Real anion_excl; /* fraction of porosity from which anions are excluded (from SWAT) */
  Real cn_ratio; /* C:N ration in soil pools */
} Soilpar;  /* soil parameters */

typedef struct
{
  const Soilpar *par; /**< pointer to soil parameters */
  Pool pool[LASTLAYER];          /**< fast and slow carbon pool for all layers*/
  Poolpar k_mean[LASTLAYER];     /**< fast and slow decay constant */
  Real *c_shift_fast[LASTLAYER]; /**< shifting rate of carbon matter to the different layer*/
  Real *c_shift_slow[LASTLAYER]; /**< shifting rate of carbon matter to the different layer*/
  Real NO3[LASTLAYER];           /**< NO3 per soillayer (gN/m2) */
  Real NH4[LASTLAYER];           /**< NH4 per soillayer (gN/m2) */
  Real w[NSOILLAYER],            /**< soil water as fraction of whc (fractional water holding capacity) */
    w_fw[NSOILLAYER];            /**< free water or gravitational water (mm), absolute water content between field capacity and saturation */
  Real w_evap;                   /**< soil moisture content which is not transpired and can evaporate? correct? */
  Real perc_energy[NSOILLAYER];  /**< energy transfer by percolation */
#ifdef MICRO_HEATING
  Real micro_heating[NSOILLAYER]; /**< energy of microbiological heating*/
  Real decomC[NSOILLAYER];
#endif
  Real alag,amp,meanw1;
  Real snowpack;
  Real snowheight; /**< height of snow */
  Real snowfraction;  /**< fraction of snow-covered ground */
  Real temp[NSOILLAYER+1];      /**< [deg C]; last layer=snow*/
  Real Ks[NSOILLAYER];    /**< saturated hydraulic conductivity (mm/h) per layer*/
  Real wpwp[NSOILLAYER];  /**< relative water content at wilting point */
  Real wfc[NSOILLAYER];   /**< relative water content at field capacity */
  Real wsat[NSOILLAYER];  /**< relative water content at saturation */
  Real whcs_all;
  Real whc[NSOILLAYER];   /**< water holding capacity (fraction), whc = wfc - wpwp */
  Real wsats[NSOILLAYER]; /**< absolute water content at saturation (mm), wsats = wsat * soildepth */
  Real whcs[NSOILLAYER];  /**< absolute water holding capacity (mm), whcs = whc * soildepth */
  Real wpwps[NSOILLAYER]; /**< water at permanent wilting point in mm, depends on soildepth */
  Real ice_depth[NSOILLAYER];   /**< mm */
  Real ice_fw[NSOILLAYER];      /**< mm */
  Real freeze_depth[NSOILLAYER]; /**< mm */
  Real ice_pwp[NSOILLAYER];      /**< fraction of water below pwp frozen */
  Real k_dry[NSOILLAYER];        /**< thermal conductivity of dry soil */
  Real bulkdens[NSOILLAYER];     /**<  bulk density of soil [kg/m3]*/
  Real df_tillage[NTILLLAYER];
  Real beta_soil[NSOILLAYER];
  short state[NSOILLAYER];
  Real maxthaw_depth;
  Real mean_maxthaw;
  Stocks decomp_litter_mean;
  int count;
  Real YEDOMA;       /**< g/m2 */
  Litter litter;     /**< Litter pool */
  Real rw_buffer;    /**< available rain water amount in buffer (mm) */
} Soil;

struct Pftpar; /* forward declaration */
struct Dailyclimate; /* forward declaration */

/* Declaration of functions
 * NOTE: initsoiltemp() in cell.h */

extern int addlitter(Litter *,const struct Pftpar *);
extern void convert_water(Soil*,int,Real*);
extern void copysoil(Soil *,const Soil *, int);
extern int findlitter(const Litter *,const struct Pftpar *);
extern Real fire_prob(const Litter *,Real);
extern unsigned int fscansoilpar(LPJfile *,Soilpar **,Verbosity);
extern Bool fscanpoolpar(LPJfile *,Poolpar *,const char *,Verbosity);
extern Bool freadlitter(FILE *,Litter *,const struct Pftpar *,int,Bool);
extern Bool freadsoil(FILE *,Soil *,const Soilpar *,const struct Pftpar *,int,Bool);
extern Bool freadsoilcode(FILE *,unsigned int *,Bool,Type);
extern void freesoil(Soil *);
extern void freelitter(Litter *);
extern void freesoilpar(Soilpar [],int);
extern void fprintlitter(FILE *,const Litter *,int);
extern void fprintsoilpar(FILE *,const Soilpar [],int,int);
extern void fprintsoil(FILE *,const Soil *,const struct Pftpar *,int,int);
extern FILE *fopensoilcode(const Filename *,Bool *,size_t *,Type *,unsigned int,Bool);
extern Bool fwritesoil(FILE *,const Soil *,int);
extern Bool fwritelitter(FILE *,const Litter *);
extern void getlag(Soil *,int);
extern int getnsoilcode(const Filename *,unsigned int,Bool);
extern Soilstate getstate(Real *); /*temperature above/below/at T_zero?*/
extern Bool initsoil(Stand *soil,const Soilpar *, int,Bool);
extern Real litter_ag_sum(const Litter *);
extern Real litter_ag_sum_n(const Litter *);
extern Real litter_agsub_sum(const Litter *);
extern Real litter_agsub_sum_n(const Litter *);
extern Real litter_ag_grass(const Litter *);
extern Real litter_ag_sum_quick(const Litter *);
extern Stocks littersom(Stand *,Real [NSOILLAYER],int,int);
extern Real littercarbon(const Litter *);
extern Stocks litterstocks(const Litter *);
extern Real moistfactor(const Litter *);
extern void moisture2soilice(Soil *,Real *,int);
extern void newsoil(Soil *);
extern int seeksoilcode(FILE *,int,size_t,Type);
extern Real snow(Soil *,Real *,Real *,Real, Real *,Real *);
extern Real snow_old(Real *,Real *,Real *,Real);
extern void soiltemp(Soil*, Real,Bool);
#ifdef COUPLING_WITH_FMS
extern void laketemp(Cell*, const struct Dailyclimate *);
#endif
extern Real soiltemp_lag(const Soil *,const Climbuf *);
extern Real soilcarbon(const Soil *);
extern Real soilcarbon_slow(const Soil *);
extern Stocks soilstocks(const Soil *);
extern Real soilwater(const Soil *);
extern Real soilconduct(const Soil *,int);
extern Real soilheatcap(const Soil *,int);
extern void soilice2moisture(Soil *, Real *,int);
extern Real temp_response(Real);
extern Real litter_ag_tree(const Litter *,int);
extern Real litter_ag_nitrogen_tree(const Litter *,int);
extern Real biologicalnfixation(const Stand *);
extern void leaching(Soil *,const Real);
extern Real volatilization(Real,Real,Real,Real,Real);
extern Real nuptake_temp_fcn(Real);
extern void denitrification(Stand *,int,int);
extern void getrootdist(Real [],const Real[],Real);
extern Stocks checklitter(Litter *);
extern void updatelitterproperties(Stand *,Real);
extern void pedotransfer(Stand *, Real *, Real *,Real);

/* Definition of macros */

#define getsoilpar(soil,var) (soil)->par->var
#define foreachsoillayer(l) for(l=0;l<NSOILLAYER;l++)
#define forrootmoist(l) for(l=0;l<3;l++)                   /* here defined for the first 1 m*/
#define forrootsoillayer(l) for(l=0;l<LASTLAYER;l++)
#define allice(soil,l) (soil->ice_depth[l]+soil->wpwps[l]*soil->ice_pwp[l]+soil->ice_fw[l])
#define allwater(soil,l) (soil->w[l]*soil->whcs[l]+soil->wpwps[l]*(1-soil->ice_pwp[l])+soil->w_fw[l])
#define timestep2sec(timestep,steps) (24.0*3600.0*((timestep)/(steps))) /* convert timestep --> sec */
#define fprintpool(file,pool) fprintf(file,"%.2f %.2f",pool.slow,pool.fast)
#define f_temp(soiltemp) exp(-(soiltemp-18.79)*(soiltemp-18.79)/(2*5.26*5.26)) /* Parton et al 2001*/
#define f_NH4(nh4) (1-exp(-0.0105*(nh4))) /* Parton et al 1996 */
#endif /* SOIL_H */
