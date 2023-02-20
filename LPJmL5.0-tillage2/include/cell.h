/**************************************************************************************/
/**                                                                                \n**/
/**                      c  e  l  l  .  h                                          \n**/
/**                                                                                \n**/
/**     C implementation of LPJmL                                                  \n**/
/**                                                                                \n**/
/**     Declaration of cell datatype                                               \n**/
/**                                                                                \n**/
/**     Hierarchy:                                                                 \n**/
/**                                                                                \n**/
/**     --cell                                                                     \n**/
/**       +-Coord                                                                  \n**/
/**       +-Managed_land                                                           \n**/
/**       | +-Landfrac                                                             \n**/
/**       | +-Resdata                                                              \n**/
/**       | | +-Reservoir                                                          \n**/
/**       | +-Cropdates                                                            \n**/
/**       | +-Manage                                                               \n**/
/**       | +-Image_data (only for version with -DIMAGE set)                       \n**/
/**       |   +-Pool                                                               \n**/
/**       +-Ignition                                                               \n**/
/**       +-Discharge                                                              \n**/
/**       | +-Queue                                                                \n**/
/**       +-Output                                                                 \n**/
/**       +-Climbuf                                                                \n**/
/**       | +-Buffer                                                               \n**/
/**       +-Balance                                                                \n**/
/**       +-Standlist                                                              \n**/
/**         +-Stand                                                                \n**/
/**           +-Pftlist                                                            \n**/
/**           | +-Pft                                                              \n**/
/**           +-Soil                                                               \n**/
/**             +-Pool                                                             \n**/
/**             +-Litter                                                           \n**/
/**               +-Litteritem                                                     \n**/
/**                 +-Trait                                                        \n**/
/**                                                                                \n**/
/** (C) Potsdam Institute for Climate Impact Research (PIK), see COPYRIGHT file    \n**/
/** authors, and contributors see AUTHORS file                                     \n**/
/** This file is part of LPJmL and licensed under GNU AGPL Version 3               \n**/
/** or later. See LICENSE file or go to http://www.gnu.org/licenses/               \n**/
/** Contact: https://github.com/PIK-LPJmL/LPJmL                                    \n**/
/**                                                                                \n**/
/**************************************************************************************/

#ifndef CELL_H /* Already included? */
#define CELL_H

/* Definition of datatypes */


typedef struct
{
  Real nep;                 /**< annual NEP (gC/m2) */
  Real awater_flux;         /**< annual water flux (mm) */
  Real aprec;               /**< annual precipitation (mm) */
  Stocks tot;                /**< total carbon and nitrogen (g/m2) */
  Stocks biomass_yield;      /**< harvested biomass (gC/m2, gN/m2)*standfrac*/
  Stocks estab_storage_tree[2];  /**< carbon and nitrogen taken out from annual NPP to satisfy tree establishment rate */
  Stocks estab_storage_grass[2]; /**< carbon and nitrogen taken out from annual NPP to satisfy grass establishment rate */
  Real totw;                /**< total water (mm) */
  Real surface_storage;     /**< total water in surface storages (dm3) */
  Real excess_water;        /**< excess water (mm) */
  Real soil_storage;        /**< total water in soil storages (dm3) */
  Real total_reservoir_out; /**< total water extracted from reservoirs (dm3) */
  Real total_irrig_from_reservoir; /**< total water added to fields from reservoirs (dm3)*/
  Real n_influx;            /**< all N inputs: deposition, fertilizer, BNF */
  Real n_outflux;           /**< all N outputs: n2onit, n2odenit, n2denit, leaching */
  Real n_demand;            /**< N demand by plants (gN)*/
  Real n_uptake;            /**< N uptake by plants (gN) */

} Balance;

typedef struct celldata *Celldata;

struct cell
{
  Coord coord;              /**< Cell coordinate and area */
  Standlist standlist;      /**< Stand list */
  Climbuf climbuf;
  Ignition ignition;
  Real afire_frac;          /**< fraction of grid cell burnt this year */
  Real *gdd;                /**< Growing degree days array */
  Real lakefrac;            /**< lake fraction (0..1) */
#ifdef COUPLING_WITH_FMS
  Real laketemp;            /**< temperatures in this cell, for now only one layer */
  Real day_after_snowfall;  /**< days after snowfall, to compute the albedo of lakes following:

Albedo models for snow and ice on a freshwater lake
Heather E. Henneman, Heinz G. Stefan
St. Anthony Falls Laboratory, Department of Civil Engineering, University of Minnesota, Minneapolis, MN 55414, USA
Received 19 November 1997; accepted 15 January 1999*/

  Real albedo_lake;         /**< for a changing albedo of lakes*/
  Real snowpool_above_lake; /**< temporarily storing the snow mass on a lake*/

#endif
  Real soilph;              /**< soil pH */
  Bool skip;                /**< Invalid soil code in cell (TRUE/FALSE) */
  Managed_land ml;          /**< Managed land */
  Output output;            /**< Output data */
  Discharge discharge;
  int elevation;            /**< cell elevation (m) */
  Balance balance;          /**< balance checks */
};

/* Declaration of functions */

extern void freegrid(Cell [],int,const Config *);
extern void freecell(Cell *,int,Bool);
extern void update_daily(Cell *,Real,Real,Dailyclimate,int,
                         int,int,int,int, Bool,Bool,const Config *);
extern void update_annual(Cell *,int,int,
                          Real,int,const Real *,Bool,Bool,const Config *);
extern void update_monthly(Cell *,Real,Real,int);
extern void init_annual(Cell *,int,int,int);
extern int fwritecell(FILE *,long long [],const Cell [],int,int,int,int,Bool);
extern void fprintcell(FILE *,const Cell [],int,int,int,const Config *);
extern Bool freadcell(FILE *,Cell *,int,int,const Soilpar *,
                      const Standtype [],int,Bool,const Config *);
extern int writecoords(Outputfile *,int,const Cell [],const Config *);
extern int writecountrycode(Outputfile *,int,const Cell [],const Config *);
extern int writeregioncode(Outputfile *,int,const Cell [],const Config *);
extern int iterate(Outputfile *,Cell [],Input,
                   int,int,const Config *);
extern void iterateyear(Outputfile *,Cell [],Input,
                        Real,int,int,int,const Config *);
extern void fwriteoutput_annual(Outputfile *,const Cell [],int,const Config *);
extern void fwriteoutput_monthly(Outputfile *,const Cell [],int,int,const Config *);
extern void fwriteoutput_daily(Outputfile *,const Cell [],int,int,const Config *);
extern void fwriteoutput_pft(Outputfile *,Cell [],int,int,int,const Config *);
extern void equilsom(Cell *,int, const Pftpar []);
extern void equilveg(Cell *);
extern void check_fluxes(Cell *,int,int,const Config *);
extern void check_balance(Flux,int,const Config *);
extern Bool initdrain(Cell [],Config *);
extern void drain(Cell [],int,const Config *);
extern Real nep_sum(const Cell [],const Config *);
extern Real cflux_sum(const Cell [],const Config *);
extern Real flux_sum(Flux *,Cell [],const Config *);
extern Bool getwateruse(Wateruse, Cell [],int,const Config *);
extern Wateruse initwateruse(const Config *);
extern void freewateruse(Wateruse,Bool);
extern void killstand(Cell *,const Pftpar [],int,Bool,Bool,int);
extern Bool initsoiltemp(Climate *, Cell*,const Config *);
extern Celldata opencelldata(Config *);
extern Bool seekcelldata(Celldata,int);
extern Bool readcelldata(Celldata,Coord *,unsigned int *,Real *,Intcoord *,int,Config *);
extern void closecelldata(Celldata);
extern Real albedo(Cell *, Real , Real );

/* Definition of macros */

#define printcell(cell,n,npft,ncft,config) fprintcell(stdout,cell,n,npft,ncft,config)

#endif
