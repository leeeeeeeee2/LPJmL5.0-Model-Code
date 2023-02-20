/**************************************************************************************/
/**                                                                                \n**/
/**                         c  r  o  p  .  h                                       \n**/
/**                                                                                \n**/
/**     C implementation of LPJmL                                                  \n**/
/**                                                                                \n**/
/**     Header file for crop PFTs                                                  \n**/
/**                                                                                \n**/
/** (C) Potsdam Institute for Climate Impact Research (PIK), see COPYRIGHT file    \n**/
/** authors, and contributors see AUTHORS file                                     \n**/
/** This file is part of LPJmL and licensed under GNU AGPL Version 3               \n**/
/** or later. See LICENSE file or go to http://www.gnu.org/licenses/               \n**/
/** Contact: https://github.com/PIK-LPJmL/LPJmL                                    \n**/
/**                                                                                \n**/
/**************************************************************************************/

#ifndef CROP_H /* Already included? */
#define CROP_H

/* Definition of constants */

#define MINLGP 4              /* minimum length of growing period, used in calc_seasonality */
#define DEFAULT_MONTH 1       /* default setting if no sowing month can be found in calc_seasonality*/
#define MIN_PREC 0.1          /* minimum daily precipitation (mm) at a wet day - definition from CRU*/

/* Declaration of datatypes */

typedef struct
{
  Real root,so,pool;
} Cropratio;

typedef struct
{
  int sdatenh,sdatesh;
} Initdate;

typedef struct
{
  Real leaf,root,so,pool;
} Cropphys;

typedef struct
{
  Stocks root,so,pool,leaf;
} Cropphys2;

typedef struct
{
  int calcmethod_sdate;     /**< different methods for determining the crop dates */
  Initdate initdate;        /**< init sowing date for n/shemisphere */
  int hlimit;               /**< maximum length of crop cycle */
  int fallow_days;          /**< length of fallow period between crop cycles */
  Real temp_fall;           /**< threshold for decreasing temperature to determine the crop date */
  Real temp_spring;         /**< threshold for increasing temperature to determine the crop date */
  Real temp_vern;           /**< threshold for increasing temperature to determine the crop date */
  Limit tv_eff;             /**< lower and upper temperature thresholds for effective vernalization (deg C) */
  Limit tv_opt;             /**< lower and upper temperature thresholds for optimal vernalization (deg C) */
  Real pvd_max;             /**< maximum number of vernalization days required */
  Real psens;               /**< sensitivity to the photoperiod effect [0-1] (1 means no sensitivity) */
  Real pb;                  /**< basal photoperiod (h) (pb<ps for longer days plants) */
  Real ps;                  /**< saturating photoperiod (h) (ps<pb for shorter days plants) */
  Limit phus;               /**< minimum/maximum potential heat units required for plant maturity (deg C), winter types*/
  Limit phuw;               /**< minimum/maximum potential heat units required for plant maturity (deg C), summer types*/
  Real phu_par;             /**< parameter for determining the variability of phu */
  Limit basetemp;           /**< minimum/maximum base temperature */
  Real fphuc;               /**< fraction of growing season 1 [0-1] */
  Real flaimaxc;            /**< fraction of plant maximal LAI 1 [0-1] */
  Real fphuk;               /**< fraction of growing season 2 [0-1] */
  Real flaimaxk;            /**< fraction of plant maximal LAI 2 [0-1] */
  Real fphusen;             /**< fraction of growing period at which LAI starts decreasing [0-1] */
  Real flaimaxharvest;      /**< fraction of plant maximal LAI still present at harvest [0-1]*/
  Real laimax;              /**< plant maximal LAI (m2leaf/m2)*/
  Real laimin;              /**< plant minimal LAI (m2leaf/m2)*/
  Real hiopt;               /**< optimum harvest index HI reached at harvest*/
  Real himin;               /**< minimum harvest index HI reached at harvest*/
  Real shapesenescencenorm; /**< parameter for calculating the fraction of maximal LAI */
  Cropphys cn_ratio;        /**< C:N mass ratio for root, storage organ, and pool */
  Cropratio ratio;          /**< relative C:N ratio for root, storage organ and pool */
} Pftcroppar;


typedef struct
{
  Bool wtype;               /**< distinguish between winter and summer crop */
  int growingdays;          /**< counter for the days of the crop cycle */
  Real pvd;                 /**< actually required vernalization days (only if STATIC_PHU) */
  Real phu;                 /**< required phenological heat units from emergence to maturity (STATIC_PHU, PRESCRIBED_PHU, CALCULATED_PHU) */
  Real basetemp;            /**< base temperature for phenological development */
  Bool senescence;          /**< current senescence period */
  Bool senescence0;         /**< senescence period of yesterday */
  Real husum;               /**< sum of heat units */
  Real vdsum;               /**< sum of vernalization days */
  Real fphu;                /**< fraction of phenological heat unit */
  Cropphys2 ind;
  Real flaimax;             /**< fraction of maximum lai */
  Real lai;                 /**< current leaf area index */
  Real lai000;              /**< leaf area index of yesterday */
  Real laimax_adjusted;     /**< adjusted maximum lai */
  Real lai_nppdeficit;      /**< LAI reduction due to insufficient NPP */
  Real demandsum;
  Real ndemandsum;
  Real nuptakesum;
  Real nfertilizer;         /* fertilizer amount */
  Real nmanure;             /* manure ammount */
  Real vscal_sum;
  Bool frostkill;           /* set to TRUE in daily_agriculture if tmin<-5 and 0.2<fphu<0.95 */
  Real supplysum;
#ifdef DOUBLE_HARVEST
  Real petsum;
  Real evapsum;
  Real transpsum;
  Real nfertsum;
  Real intercsum;
  Real precsum;
  Real sradsum;
  Real irrig_apply;
  Real tempsum;
  Real nirsum;
  Real lgp;
  Real runoffsum;
  Real n2o_denitsum;
  Real n2o_nitsum;
  Real n2_emissum;
  Real leachingsum;
  Real c_emissum;
  int sdate;
  int sowing_year;
#endif
} Pftcrop;

/* Declaration of functions */

extern void new_crop(Pft *,int,int,int);
extern void allocation_daily_crop(Pft *,Real, Real,int,Daily_outputs *);
extern Real npp_crop(Pft *,Real,Real,Real,Bool *,Real,int,Daily_outputs *);
extern Real fpc_crop(Pft *);
extern Real fpar_crop(const Pft *);
extern Real alphaa_crop(const Pft *,int,int);
extern void litter_update_crop(Litter *,Pft *,Real);
extern Real lai_crop(const Pft *);
extern Real actual_lai_crop(const Pft *);
extern Bool phenology_crop(Pft *,Real,Real,Real,int,const Config *);
extern void laimax_manage(Manage *,const Pftpar [],int,int,int);
extern Bool fwrite_crop(FILE *,const Pft *);
extern void fprint_crop(FILE *,const Pft *,int);
extern Bool fread_crop(FILE *,Pft *,Bool);
extern Bool fscanpft_crop(LPJfile *,Pftpar *,Verbosity);
extern Stocks establishment_crop(Pft *,Real,Real,int);
extern void init_crop(Pft *);
extern Real vegc_sum_crop(const Pft *);
extern Real vegn_sum_crop(const Pft *);
extern Real agb_crop(const Pft *);
extern void free_crop(Pft *);
extern void phen_variety(Pft *,int,Real,int,Bool,const Config *,int,int);
extern void harvest_crop(Output *,Stand *,Pft *,int,int,int,int,Bool,Bool);
extern void adapt_crop_type(Real [],Real,const Pftpar [],int,int,int);
extern Real wdf_crop(Pft *,Real,Real);
extern void fprintpar_crop(FILE *,const Pftpar *);
extern void output_daily_crop(Daily_outputs *,const Pft *,Real,Real);
extern void calc_seasonality(Cell *,int,int,const Config *);
extern void albedo_crop(Pft *,Real,Real);
extern void double_harvest(int, Real *, Real *, Real);
extern Real nuptake_crop(Pft *,Real *,Real *,int,int,int,int,Bool);
extern Real ndemand_crop(const Pft *,Real *,Real,Real,Real);
extern Real vmaxlimit_crop(const Pft *,Real,Real);


/* Definitions of macros */

#define iscrop(pft) (getpftpar(pft,type)==CROP)
#define phys_sum_crop(ind) (ind.leaf.carbon+ind.root.carbon+ind.so.carbon+ind.pool.carbon)
#define phys_sum_crop_n(ind) (ind.leaf.nitrogen+ind.root.nitrogen+ind.so.nitrogen+ind.pool.nitrogen)
#define fprintcropphys2(file,phys,nind) fprintf(file,"%.2f %.2f %.2f %.2f (gC/m2) %.2f %.2f %.2f %.2f (gN/m2)",phys.leaf.carbon*nind,phys.so.carbon*nind,phys.pool.carbon*nind,phys.root.carbon*nind,phys.leaf.nitrogen*nind,phys.so.nitrogen*nind,phys.pool.nitrogen*nind,phys.root.nitrogen*nind)
#define fprintcropphys2carbon(file,phys,nind) fprintf(file,"%.2f %.2f %.2f %.2f (gC/m2)",phys.leaf.carbon*nind,phys.so.carbon*nind,phys.pool.carbon*nind,phys.root.carbon*nind)

#endif
