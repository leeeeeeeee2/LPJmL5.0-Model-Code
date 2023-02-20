/**************************************************************************************/
/**                                                                                \n**/
/**                 l  a  n  d  u  s  e  .  h                                      \n**/
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

#ifndef LANDUSE_H /* Already included? */
#define LANDUSE_H

/* Definitions of datatypes */

typedef enum {NATURAL,SETASIDE_RF,SETASIDE_IR,AGRICULTURE,MANAGEDFOREST,
              GRASSLAND,BIOMASS_TREE,BIOMASS_GRASS,KILL} Landusetype;

typedef struct landuse *Landuse;

typedef struct
{
  Real *crop;
  Real grass[NGRASS];
  Real biomass_grass;
  Real biomass_tree;
} Landfrac;

typedef struct
{
  int *crop;
  int grass[NGRASS];
  int biomass_grass;
  int biomass_tree;
} Irrig_system;

typedef enum {NO_SEASONALITY, PREC, PRECTEMP, TEMP, TEMPPREC} Seasonality;
typedef enum {GS_DEFAULT, GS_MOWING, GS_GRAZING_EXT, GS_GRAZING_INT, GS_NONE} GrassScenarioType;
typedef enum {RM_UNDEFINED, RM_GRAZING, RM_RECOVERY} RotationModeType;

typedef struct
{
  int grazing_days;
  int recovery_days;
  int paddocks;
  RotationModeType rotation_mode;
} Rotation;

typedef struct
{
  Landfrac *landfrac;     /**< land use fractions */
  Landfrac *fertilizer_nr;   /**< reactive nitrogen fertilizer */
  Landfrac *manure_nr;     /* manure nitrogen fertilizer */
  Landfrac *residue_on_field;     /* fraction of residues left on field */
  Irrig_system *irrig_system; /**< irrigation system type (SURF=1,SPRINK=2,DRIP=3) */
  Manage manage;
  Cropdates *cropdates;
  Real cropfrac_rf;       /**< rain-fed crop fraction (0..1) */
  Real cropfrac_ir;       /**< irrigated crop fraction (0..1) */
  int *sowing_month;      /**< sowing month (index of month, 1..12), rainfed, irrigated*/
  int *gs;                /**< length of growing season (number of consecutive months, 0..11)*/
  Seasonality seasonality_type;  /**< seasonality type (0..4) 0:no seasonality, 1 and
                                    2:precipitation seasonality, 3 and 4:temperature
                                    seasonality*/
  int *sdate_fixed;       /**< array to store fixed or prescribed sowing dates */
  Real *crop_phu_fixed;   /**< array to store prescribed crop phus */
  int *tilltypes;         /* array to store tillage types */
  Resdata *resdata;       /**< Reservoir data */
  Real *fraction;
  Real reservoirfrac;     /**< reservoir fraction (0..1) */
  Real mdemand;           /**< monthly irrigation demand */
  Bool dam;               /**< dam inside cell (TRUE/FALSE) */
  Bool with_tillage;      /* simulation with tillage implementation */
  int fixed_grass_pft;              /**< fix C3 or C4 for GRASS pft */
  GrassScenarioType grass_scenario; /* 0=default, 1=mowing, 2=ext.grazing, 3=int.grazing */
  Real nr_of_lsus_ext;              /* nr of livestock units for extensive grazing */
  Real nr_of_lsus_int;              /* nr of livestock units for intensive grazing */
  Real deficit_lsu_ne;              /* NE deficit per livestock unit */
  Real deficit_lsu_mp;              /* MP deficit per livestock unit */
  Rotation rotation;                /* rotation mode and parameters for intensive grazing */

#ifdef IMAGE
  Image_data *image_data; /**< pointer to IMAGE data structure */
#else
  Pool product;
#endif
} Managed_land;

/* Definitions of macros */

#define mixpool(pool1,pool2,frac1,frac2) pool1=(pool1*frac1+pool2*frac2)/(frac1+frac2)
#define rothers(ncft) ncft
#define rmgrass(ncft) (ncft+1)
#define rbgrass(ncft) (ncft+2)
#define rbtree(ncft) (ncft+3)

/* Declaration of functions */

extern Landuse initlanduse(int,const Config *);
extern void freelanduse(Landuse,Bool);
extern Bool getintercrop(const Landuse);
extern void newlandfrac(Landfrac [2],int);
extern void initlandfrac(Landfrac [2],int);
extern void scalelandfrac(Landfrac [2],int,Real);
extern void freelandfrac(Landfrac [2]);
extern Bool fwritelandfrac(FILE *,const Landfrac [2],int);
extern Bool freadlandfrac(FILE *,Landfrac [2],int,Bool);
extern Real landfrac_sum(const Landfrac [2],int,Bool);
extern Real crop_sum_frac(Landfrac *,int,Real,Bool);
extern Stocks cultivate(Cell *,const Pftpar *,int,Real,Bool,int,Bool,Stand *,
                      Bool,Bool,const Config *,int,int,int,int);
extern void reclaim_land(const Stand *, Stand *,Cell *,Bool,int);
extern Bool getlanduse(Landuse,Cell *,int,int,int,const Config *);
extern void landusechange(Cell *,const Pftpar[],int,int,int,Bool,Bool,Bool,int,const Config *);
extern Bool setaside(Cell *,Stand *,const Pftpar[],Bool,Bool,int,Bool,int,int);
extern Stocks sowing_season(Cell *,int,int,int,Real,int,const Config *);
extern Stocks sowing_prescribe(Cell *,int,int,int,int,const Config *);
extern Stocks sowing(Cell *,Real,int,int,int,int,const Config *);
extern void deforest(Cell *,Real,const Pftpar [],Bool,int,Bool,Bool,Bool,int,int,Real,int);
extern Stocks woodconsum(Stand*,Real);
extern void calc_nir(Stand *,Real,Real [],Real);
extern Real rw_irrigation(Stand *,Real,const Real [],Real);
extern void irrig_amount_river(Cell *,const Config *);
extern void irrig_amount(Stand *,Bool,int,int);
extern void mixsetaside(Stand *,Stand *,Bool);
extern void set_irrigsystem(Stand *,int,int,Bool);
extern void tillage(Soil *, Real);

/* Declaration of variables */

extern Real tinyfrac;

#endif
