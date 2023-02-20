/**************************************************************************************/
/**                                                                                \n**/
/**                   f  s  c  a  n  c  o  n  f  i  g  .  c                        \n**/
/**                                                                                \n**/
/**     C implementation of LPJmL                                                  \n**/
/**                                                                                \n**/
/**     Function reads configuration from an ascii file                            \n**/
/**                                                                                \n**/
/** (C) Potsdam Institute for Climate Impact Research (PIK), see COPYRIGHT file    \n**/
/** authors, and contributors see AUTHORS file                                     \n**/
/** This file is part of LPJmL and licensed under GNU AGPL Version 3               \n**/
/** or later. See LICENSE file or go to http://www.gnu.org/licenses/               \n**/
/** Contact: https://github.com/PIK-LPJmL/LPJmL                                    \n**/
/**                                                                                \n**/
/**************************************************************************************/

#include "lpj.h"

#define fscanint2(file,var,name) if(fscanint(file,var,name,FALSE,verbose)) return TRUE;
#define fscanreal2(file,var,name) if(fscanreal(file,var,name,FALSE,verbose)) return TRUE;
#define fscanbool2(file,var,name) if(fscanbool(file,var,name,FALSE,verbose)) return TRUE;
#define fscanname(file,var,name) {              \
    if(fscanstring(file,var,name,FALSE,verbose)) {                 \
      if(verbose) readstringerr(name);  \
      return TRUE;                              \
    }                                              \
  }

#define scanfilename(file,var,path,what) {                              \
    if(readfilename2(file,var,what,path,verbose)) {                  \
      if(verbose) fprintf(stderr,"ERROR209: Cannot read input filename for '%s'.\n",what); \
      return TRUE;                                                      \
    }                                                                   \
    if(verbose>=VERB)                      \
      printf("%s %s\n", what, (var)->name);                             \
  }

#define scanclimatefilename(file,var,path,isfms,what) {                 \
    if(readclimatefilename(file,var,what,path,isfms,verbose)) {      \
      if(verbose) fprintf(stderr,"ERROR209: Cannot read input filename for '%s'.\n",what); \
      return TRUE;                                                      \
    }                                                                   \
    if(verbose>=VERB)                      \
      printf("%s %s\n", what, (var)->name);                             \
  }

static Bool readfilename2(LPJfile *file,Filename *name,const char *key,const char *path,Verbosity verbose)
{
  if(readfilename(file,name,key,path,FALSE,verbose))
    return TRUE;
  if(name->fmt==CDF)
  {
    if(verbose)
      fprintf(stderr,"ERROR197: NetCDF is not supported for input '%s' in this version of LPJmL.\n",name->name);
    return TRUE;
  }
  else if(name->fmt==TXT)
  {
    if(verbose)
      fprintf(stderr,"ERROR197: text file is not supported for input '%s' in this version of LPJmL.\n",name->name);
    return TRUE;
  }
  return FALSE;
} /* of 'readfilename2' */

static Bool readclimatefilename(LPJfile *file,Filename *name,const char *key,const char *path,Bool isfms,Verbosity verbose)
{
  if(readfilename(file,name,key,path,TRUE,verbose))
    return TRUE;
  if(!isfms && name->fmt==FMS)
  {
    if(verbose)
      fprintf(stderr,"ERROR197: FMS coupler not allowed for input '%s'.\n",name->name);
    return TRUE;
  }
  if(name->fmt==TXT)
  {
    if(verbose)
      fprintf(stderr,"ERROR197: text file is not supported for input '%s' in this version of LPJmL.\n",name->name);
    return TRUE;
  }
  return FALSE;
} /* of 'readclimatefilename' */


static void divide(int *start, /**< index of first grid cell */
                   int *end,   /**< index of last grid cell */
                   int rank,   /**< my rank id */
                   int ntask   /**< total number of tasks */
                  )
{
/*
 * Function is used in the parallel MPI version to distribute the cell grid
 * equally on ntask tasks
 * On return start and end are set to the local boundaries of each task
 */
  int i,lo,hi,n;
  n=*end-*start+1;
  lo=*start;
  hi=*start+n/ntask-1;
  if(n % ntask)
    hi++;
  for(i=1;i<=rank;i++)
  {
    lo=hi+1;
    hi=lo+n/ntask-1;
    if(n % ntask>i)
      hi++;
  }
  *start=lo;
  *end=hi;
} /* of 'divide' */

Bool fscanconfig(Config *config,    /**< LPJ configuration */
                 LPJfile *file,        /**< File pointer to LPJ file */
                 Fscanpftparfcn scanfcn[], /**< array of PFT-specific scan
                                              functions */
                 int ntypes,        /**< Number of PFT classes */
                 int nout_max       /**< maximum number of output files */
                )                   /** \return TRUE on error */
 {
  String name;
  LPJfile input;
  int restart,endgrid,israndom,grassfix,grassharvest;
  Verbosity verbose;

  verbose=(isroot(*config)) ? config->scan_verbose : NO_ERR;

  /*=================================================================*/
  /* I. Reading type section                                         */
  /*=================================================================*/

  if (verbose>=VERB) puts("// I. type section");
  fscanbool2(file,&israndom,"random_prec");
  if(israndom)
  {
    config->seed=RANDOM_SEED;
    if(fscanint(file,&config->seed,"random_seed",TRUE,verbose))
      return TRUE;
    if(config->seed==RANDOM_SEED)
      config->seed=time(NULL);
  }
  else
    config->seed=0;
  config->with_nitrogen=LIM_NITROGEN;
  if(fscanint(file,&config->with_nitrogen,"with_nitrogen",TRUE,verbose))
    return TRUE;
  fscanint2(file,&config->with_radiation,"radiation");
  if(config->with_radiation<CLOUDINESS || config->with_radiation>RADIATION_LWDOWN)
  {
    if(verbose)
      fprintf(stderr,"ERROR219: Invalid radiation model %d.\n",config->with_radiation);
    return TRUE;
  }
#ifdef IMAGE
  if(config->sim_id==LPJML_IMAGE && config->with_radiation)
  {
    if(verbose)
      fputs("ERROR208: Radiation data not supported by IMAGE coupler.\n",
            stderr);
    return TRUE;
  }
#endif
  fscanint2(file,&config->fire,"fire");
  if(config->fire<NO_FIRE || config->fire>SPITFIRE_TMAX)
  {
    if(verbose)
      fprintf(stderr,"ERROR166: Invalid value for fire=%d.\n",config->fire);
    return TRUE;
  }
  if(config->sim_id==LPJ)
    config->firewood=NO_FIREWOOD;
  else
  {
    fscanbool2(file,&config->firewood,"firewood");
  }
  fscanbool2(file,&config->ispopulation,"population");
  config->prescribe_burntarea=FALSE;
  if(fscanbool(file,&config->prescribe_burntarea,"prescribe_burntarea",TRUE,verbose))
    return TRUE;
  if(config->prescribe_burntarea && config->fire!=SPITFIRE && config->fire!=SPITFIRE_TMAX)
  {
    if(verbose)
      fputs("WARNING029: Prescribed burnt area can only by set for SPITFIRE, will be disabled.\n",stderr);
    config->prescribe_burntarea=FALSE;
  }
  config->prescribe_landcover=NO_LANDCOVER;
  if(fscanint(file,&config->prescribe_landcover,"prescribe_landcover",TRUE,verbose))
    return TRUE;
  if(config->prescribe_landcover<NO_LANDCOVER || config->prescribe_landcover>LANDCOVERFPC)
  {
    if(verbose)
      fprintf(stderr,"ERROR166: Invalid value for prescribe landcover=%d.\n",
              config->prescribe_landcover);
    return TRUE;
  }
  fscanbool2(file,&config->new_phenology,"new_phenology");
  fscanbool2(file,&config->river_routing,"river_routing");
  config->reservoir=FALSE;
  fscanbool2(file,&config->permafrost,"permafrost");
  config->sdate_option=NO_FIXED_SDATE;
  config->rw_manage=FALSE;
  config->const_climate=FALSE;
  config->tillage_type=NO_TILLAGE;
  config->residue_treatment=NO_RESIDUE_REMOVE;
  config->no_ndeposition=FALSE;
  if(fscanbool(file,&config->const_climate,"const_climate",TRUE,verbose))
    return TRUE;
  config->const_deposition=FALSE;
  if(config->with_nitrogen==LIM_NITROGEN)
  {
    if(fscanbool(file,&config->const_deposition,"const_deposition",TRUE,verbose))
      return TRUE;
    if(fscanbool(file,&config->no_ndeposition,"no_ndeposition",TRUE,verbose))
      return TRUE;
  }
  if(config->sim_id!=LPJ)
  {
    fscanint2(file,&config->withlanduse,"landuse");
    if(config->withlanduse<NO_LANDUSE || config->withlanduse>ONLY_CROPS)
    {
      if(verbose)
        fprintf(stderr,"ERROR166: Invalid value for landuse=%d.\n",
                config->withlanduse);
      return TRUE;
    }
    if(config->withlanduse!=NO_LANDUSE)
    {
      if(config->withlanduse==CONST_LANDUSE || config->withlanduse==ONLY_CROPS)
        fscanint2(file,&config->landuse_year_const,"landuse_year_const");
      fscanint2(file,&config->sdate_option,"sowing_date_option");
      if(config->sdate_option<0 || config->sdate_option>PRESCRIBED_SDATE)
      {
        if(verbose)
          fprintf(stderr,"ERROR166: Invalid value for sowing date option=%d.\n",
                  config->sdate_option);
        return TRUE;
      }
      if(config->sdate_option==FIXED_SDATE || config->sdate_option==PRESCRIBED_SDATE)
        fscanint2(file,&config->sdate_fixyear,"sdate_fixyear");
      fscanint2(file,&config->irrig_scenario,"irrigation");
      if(config->irrig_scenario<0 || config->irrig_scenario>ALL_IRRIGATION)
      {
        if(verbose)
          fprintf(stderr,"ERROR166: Invalid value for irrigation scenario=%d.\n",
                  config->irrig_scenario);
        return TRUE;
      }
      fscanbool2(file,&config->intercrop,"intercrop");
      if(config->with_nitrogen)
      {
        config->fertilizer_input=TRUE;
        if(fscanbool(file,&config->fertilizer_input,"fertilizer_input",TRUE,verbose))
          return TRUE;
        config->manure_input=FALSE;
        if (fscanbool(file,&config->manure_input,"manure_input",TRUE,verbose))
          return TRUE;
        config->fix_fertilization=FALSE;
        if(fscanbool(file,&config->fix_fertilization,"fix_fertilization",TRUE,verbose))
          return TRUE;
      }
      config->others_to_crop = FALSE;
      if (fscanbool(file, &config->others_to_crop, "others_to_crop", TRUE, verbose))
        return TRUE;
      config->grassonly = FALSE;
      if (fscanbool(file, &config->grassonly, "grassonly", TRUE, verbose))
        return TRUE;
      config->istimber=FALSE;
      if(fscanbool(file,&config->istimber,"istimber",TRUE,verbose))
        return TRUE;
      config->residues_fire=FALSE;
      if(fscanbool(file,&config->residues_fire,"residues_fire",TRUE,verbose))
        return TRUE;
      if(fscanbool(file,&config->rw_manage,"rw_manage",TRUE,verbose))
        return TRUE;
      fscanint2(file,&config->laimax_interpolate,"laimax_interpolate");
      if(config->laimax_interpolate<0 || config->laimax_interpolate>LAIMAX_PAR)
      {
        if(verbose)
          fprintf(stderr,"ERROR166: Invalid value for laimax_interpolate=%d.\n",
                  config->laimax_interpolate);
        return TRUE;
      }
      if(config->laimax_interpolate==CONST_LAI_MAX)
        fscanreal2(file,&config->laimax,"laimax");
      if(config->river_routing)
        fscanbool2(file,&config->reservoir,"reservoir");
      grassfix=FALSE;
      if(fscanbool(file,&grassfix,"grassland_fixed_pft",TRUE,verbose))
        return TRUE;
      grassharvest=FALSE;
      if(fscanbool(file,&grassharvest,"grass_harvest_options",TRUE,verbose))
        return TRUE;
      fscanint2(file,&config->tillage_type,"tillage_type");
      fscanint2(file,&config->residue_treatment,"residue_treatment");
    }
    config->black_fallow=FALSE;
    if(fscanbool(file,&config->black_fallow,"black_fallow",TRUE,verbose))
      return TRUE;
    if(config->black_fallow)
    {
      fscanbool2(file,&config->till_fallow,"till_fallow");
      fscanbool2(file,&config->prescribe_residues,"prescribe_residues");
    }
    if(isboolean(file,"wateruse"))
    {
      if(verbose)
        fputs("WARNING028: Type of 'wateruse' is boolean, converted to int.\n",stderr);
      fscanbool2(file,&config->wateruse,"wateruse");
    }
    else
    {
      fscanint2(file,&config->wateruse,"wateruse");
      if(config->wateruse<NO_WATERUSE || config->wateruse>ALL_WATERUSE)
      {
        if(verbose)
          fprintf(stderr,"ERROR166: Invalid value for wateruse=%d.\n",
                  config->wateruse);
        return TRUE;
      }
    }
    if(config->wateruse && config->withlanduse==NO_LANDUSE)
    {
      if(verbose)
        fputs("ERROR224: Wateruse without landuse set.\n",stderr);
      return TRUE;
    }
  }
  else
  {
    config->withlanduse=NO_LANDUSE;
    config->wateruse=NO_WATERUSE;
  }
  /*=================================================================*/
  /* II. Reading input parameter section                             */
  /*=================================================================*/

  if (verbose>=VERB) puts("// II. input parameter section");

  /* ntypes is set to natural vegetation must be considered
   * in light and establishment
   * crops must have last id-number */
  /* Read PFT parameter array */
  if(fscanparam(file,config))
  {
    if(verbose)
      fputs("ERROR230: Cannot read LPJ parameter.\n",stderr);
    return TRUE;
  }
  if((config->nsoil=fscansoilpar(file,&config->soilpar,verbose))==0)
  {
    if(verbose)
      fputs("ERROR230: Cannot read soil parameter.\n",stderr);
    return TRUE;
  }
  if((config->npft=fscanpftpar(file,&config->pftpar,scanfcn,ntypes,verbose))==NULL)
  {
    if(verbose)
      fputs("ERROR230: Cannot read PFT parameter.\n",stderr);
    return TRUE;
  }
  config->ntypes=ntypes;
  config->nbiomass=getnbiomass(config->pftpar,config->npft[GRASS]+config->npft[TREE]);
  /* Read soil paramater array */
  if(config->withlanduse!=NO_LANDUSE)
  {
    /* landuse enabled */
    if((config->ncountries=fscancountrypar(file,&config->countrypar,config->rw_manage,(config->laimax_interpolate==LAIMAX_CFT) ? config->npft[CROP] : 0,verbose))==0)
    {
      if(verbose)
        fputs("ERROR230: Cannot read country parameter.\n",stderr);
      return TRUE;
    }
    if((config->nregions=fscanregionpar(file,&config->regionpar,verbose))==0)
    {
      if(verbose)
        fputs("ERROR230: Cannot read region parameter.\n",stderr);
      return TRUE;
    }
  }
  else
  {
    config->countrypar=NULL;
    config->regionpar=NULL;
  }
  config->compress=0;
  if(fscanint(file,&config->compress,"compress",TRUE,verbose))
    return TRUE;
#ifdef USE_NETCDF
  if(config->compress && verbose)
    fputs("WARNING403: Compression of NetCDF files is not supported in this version of NetCDF.\n",stderr);
#endif
  config->missing_value=MISSING_VALUE_FLOAT;
  if(fscanfloat(file,&config->missing_value,"missing_value",TRUE,verbose))
    return TRUE;
  config->outnames=fscanoutputvar(file,NOUT,verbose);
  if(config->outnames==NULL)
  {
    if(verbose)
      fputs("ERROR230: Cannot read output description.\n",stderr);
    return TRUE;
  }
  /*=================================================================*/
  /* III. Reading input data section                                 */
  /*=================================================================*/

  if (verbose>=VERB) puts("// III. input data section");
  config->check_climate=FALSE;
  if(iskeydefined(file,"check_climate"))
  {
    fscanbool2(file,&config->check_climate,"check_climate");
  }
  if(iskeydefined(file,"inpath"))
  {
    if(fscanstring(file,name,"inpath",FALSE,verbose))
      return TRUE;
    free(config->inputdir);
    config->inputdir=strdup(name);
  }
  if(fscanstruct(file,&input,"input",verbose))
    return TRUE;
  scanclimatefilename(&input,&config->soil_filename,config->inputdir,FALSE,"soil");
  if(config->soil_filename.fmt!=CDF)
  {
    scanfilename(&input,&config->coord_filename,config->inputdir,"coord");
  }
  if(config->withlanduse!=NO_LANDUSE)
  {
    scanclimatefilename(&input,&config->countrycode_filename,config->inputdir,FALSE,"countrycode");
    if(config->countrycode_filename.fmt==CDF)
    {
      scanclimatefilename(&input,&config->regioncode_filename,config->inputdir,FALSE,"regioncode");
    }
    scanclimatefilename(&input,&config->landuse_filename,config->inputdir,FALSE,"landuse");
    if(config->sdate_option==PRESCRIBED_SDATE)
    {
      scanclimatefilename(&input,&config->sdate_filename,config->inputdir,FALSE,"sdate");
    }
    if(config->with_nitrogen && config->fertilizer_input)
      scanclimatefilename(&input,&config->fertilizer_nr_filename,config->inputdir,FALSE,"fertilizer_nr");
    if (config->with_nitrogen && config->manure_input)
      scanclimatefilename(&input,&config->manure_nr_filename,config->inputdir,FALSE,"manure_nr");
    if (config->tillage_type==READ_TILLAGE)
      scanclimatefilename(&input,&config->with_tillage_filename,config->inputdir,FALSE,"with_tillage");
    if (config->residue_treatment == READ_RESIDUE_DATA)
      scanclimatefilename(&input,&config->residue_data_filename,config->inputdir,FALSE,"residue_on_field");
    if(grassfix == GRASS_FIXED_PFT)
    {
      scanclimatefilename(&input,&config->grassfix_filename,config->inputdir,FALSE,"grassland_fixed_pft");
    }
    else
      config->grassfix_filename.name = NULL;
    if(grassharvest)
    {
      scanclimatefilename(&input,&config->grassharvest_filename,config->inputdir,FALSE,"Grassland harvest options");
    }
    else
      config->grassharvest_filename.name = NULL;

  }
  else
  {
    config->grassfix_filename.name = NULL;
    config->grassharvest_filename.name = NULL;
  }
  if(config->river_routing)
  {
    scanclimatefilename(&input,&config->lakes_filename,config->inputdir,FALSE,"lakes");
    scanclimatefilename(&input,&config->drainage_filename,config->inputdir,FALSE,"drainage");
    if(config->drainage_filename.fmt==CDF)
    {
      scanclimatefilename(&input,&config->river_filename,config->inputdir,FALSE,"river");
    }
    if(config->withlanduse!=NO_LANDUSE)
    {
      scanclimatefilename(&input,&config->neighb_irrig_filename,config->inputdir,FALSE,"neighb_irrig");
      if(config->reservoir)
      {
        scanclimatefilename(&input,&config->elevation_filename,config->inputdir,FALSE,"elevation");
        scanfilename(&input,&config->reservoir_filename,config->inputdir,"reservoir");
      }
    }
    if(config->sim_id==LPJML_FMS)
    {
      scanfilename(&input,&config->runoff2ocean_filename,config->inputdir,"runoff2ocean_map");
    }
  }
  scanclimatefilename(&input,&config->temp_filename,config->inputdir,config->sim_id==LPJML_FMS,"temp");
  scanclimatefilename(&input,&config->prec_filename,config->inputdir,config->sim_id==LPJML_FMS,"prec");
  switch(config->with_radiation)
  {
    case RADIATION:
      scanclimatefilename(&input,&config->lwnet_filename,config->inputdir,config->sim_id==LPJML_FMS,"lwnet");
      scanclimatefilename(&input,&config->swdown_filename,config->inputdir,config->sim_id==LPJML_FMS,"swdown");
      break;
    case RADIATION_LWDOWN:
      scanclimatefilename(&input,&config->lwnet_filename,config->inputdir,config->sim_id==LPJML_FMS,"lwdown");
      scanclimatefilename(&input,&config->swdown_filename,config->inputdir,config->sim_id==LPJML_FMS,"swdown");
      break;
    case CLOUDINESS:
      scanclimatefilename(&input,&config->cloud_filename,config->inputdir,config->sim_id==LPJML_FMS,"cloud");
      break;
    case RADIATION_SWONLY:
      scanclimatefilename(&input,&config->swdown_filename,config->inputdir,config->sim_id==LPJML_FMS,"swdown");
      break;
    default:
      if(verbose)
        fprintf(stderr,"ERROR213: Invalid setting %d for radiation.\n",config->with_radiation);
      return TRUE;
  }
  if(config->with_nitrogen)
  {
    if(config->with_nitrogen!=UNLIM_NITROGEN)
    {
      scanclimatefilename(&input,&config->no3deposition_filename,config->inputdir,config->sim_id==LPJML_FMS,"no3deposition");
      scanclimatefilename(&input,&config->nh4deposition_filename,config->inputdir,config->sim_id==LPJML_FMS,"nh4deposition");
    }
    scanclimatefilename(&input,&config->soilph_filename,config->inputdir,config->sim_id==LPJML_FMS,"soilpH");
  }
  else
    config->no3deposition_filename.name=config->nh4deposition_filename.name=config->soilph_filename.name=NULL;
  if(config->with_nitrogen || config->fire==SPITFIRE || config->fire==SPITFIRE_TMAX)
  {
    scanclimatefilename(&input,&config->wind_filename,config->inputdir,config->sim_id==LPJML_FMS,"wind");
  }
  if(config->fire==SPITFIRE || config->fire==SPITFIRE_TMAX)
  {
    scanclimatefilename(&input,&config->tamp_filename,config->inputdir,config->sim_id==LPJML_FMS,(config->fire==SPITFIRE_TMAX) ? "tmin" : "tamp");
    if(config->fire==SPITFIRE_TMAX)
    {
      scanclimatefilename(&input,&config->tmax_filename,config->inputdir,config->sim_id==LPJML_FMS,"tmax");
    }
    else
      config->tmax_filename.name=NULL;
    scanclimatefilename(&input,&config->lightning_filename,config->inputdir,FALSE,"lightning");
    scanclimatefilename(&input,&config->human_ignition_filename,
                        config->inputdir,FALSE,"human_ignition");
  }
  if(config->ispopulation)
  {
    scanclimatefilename(&input,&config->popdens_filename,config->inputdir,FALSE,"popdens");
  }
  if(config->prescribe_burntarea)
  {
    scanclimatefilename(&input,&config->burntarea_filename,config->inputdir,FALSE,"burntarea");
  }
  if(config->prescribe_landcover!=NO_LANDCOVER)
  {
    scanclimatefilename(&input,&config->landcover_filename,config->inputdir,FALSE,"landcover");
  }
  if(readfilename(&input,&config->co2_filename,"co2",config->inputdir,FALSE,verbose))
  {
    if(verbose)
      fputs("ERROR209: Cannot read input filename for 'co2'.\n",stderr);
    return TRUE;
  }
  if(config->co2_filename.fmt!=TXT &&  (config->sim_id!=LPJML_FMS || config->co2_filename.fmt!=FMS))
  {
    if(verbose)
      fputs("ERROR197: Only text file is supported for CO2 input in this version of LPJmL.\n",stderr);
    return TRUE;
  }

  if(israndom==RANDOM_PREC)
  {
    scanclimatefilename(&input,&config->wet_filename,config->inputdir,config->sim_id==LPJML_FMS,"wetdays");
  }
#ifdef IMAGE
  else if(config->sim_id==LPJML_IMAGE)
  {
    if(verbose)
      fputs("ERROR175: Number of wet days must be supplied for IMAGE coupler.\n",stderr);
    return TRUE;
  }
#endif
  else
    config->wet_filename.name=NULL;
  if(config->wateruse)
  {
    scanclimatefilename(&input,&config->wateruse_filename,config->inputdir,FALSE,"wateruse");
  }
#ifdef IMAGE
  if(config->sim_id==LPJML_IMAGE)
  {
    /* reading IMAGE-coupling specific information */
    scanclimatefilename(&input,&config->temp_var_filename,config->inputdir,FALSE,"temp_var");
    scanclimatefilename(&input,&config->prec_var_filename,config->inputdir,FALSE,"prec_var");
    scanclimatefilename(&input,&config->prodpool_init_filename,config->inputdir,FALSE,"prodpool_init");
  }
#endif

  /*=================================================================*/
  /* IV. Reading output data section                                 */
  /*=================================================================*/

  if (verbose>=VERB) puts("// IV. output data section");
  if(fscanoutput(file,config,nout_max))
  {
    if(verbose)
      fputs("ERROR230: Cannot read output data.\n",stderr);
    return TRUE;
  }

  /*=================================================================*/
  /* V. Reading run settings section                                 */
  /*=================================================================*/

  if (verbose>=VERB) puts("// V. run settings");
  if(iskeydefined(file,"restartpath"))
  {
    if(fscanstring(file,name,"restartpath",FALSE,verbose))
      return TRUE;
    free(config->restartdir);
    config->restartdir=strdup(name);
  }
  config->startgrid=ALL; /* set default value */
  if(fscanint(file,&config->startgrid,"startgrid",TRUE,verbose))
    return TRUE;
  if(config->startgrid==ALL)
  {
    config->startgrid=0;
    endgrid=getnsoilcode(&config->soil_filename,config->nsoil,isroot(*config));
    if(endgrid==-1)
      return TRUE;
    else if(endgrid==0)
    {
      if(verbose)
        fputs("ERROR135: No soil code found.\n",stderr);
      return TRUE;
    }
    else
      endgrid--;
  }
  else
  {
    endgrid=config->startgrid;
    if(fscanint(file,&endgrid,"endgrid",TRUE,verbose))
      return TRUE;
  }
  if(endgrid<config->startgrid)
  {
    if(verbose)
      fprintf(stderr,"ERROR136: Endgrid=%d less than startgrid=%d.\n",
              endgrid,config->startgrid);
    return TRUE;
  }
  config->nall=endgrid-config->startgrid+1;
  config->firstgrid=config->startgrid;
  if(config->nall<config->ntask)
  {
    if(verbose)
      fprintf(stderr,"ERROR198: Number of cells %d less than number of tasks.\n",config->nall);
    return TRUE;
  }
  if(config->ntask>1) /* parallel mode? */
    divide(&config->startgrid,&endgrid,config->rank,
           config->ntask);
  config->ngridcell=endgrid-config->startgrid+1;
  fscanint2(file,&config->nspinup,"nspinup");
  config->isfirstspinupyear=FALSE;
  if(config->nspinup)
  {
    fscanint2(file,&config->nspinyear,"nspinyear");
    if(iskeydefined(file,"firstspinupyear"))
    {
      fscanint2(file,&config->firstspinupyear,"firstspinupyear");
      config->isfirstspinupyear=TRUE;
    }
  }
  fscanint2(file,&config->firstyear,"firstyear");
  fscanint2(file,&config->lastyear,"lastyear");
#ifdef IMAGE
  if(config->sim_id==LPJML_IMAGE)
  {
    fscanint2(file,&config->start_imagecoupling,"start_imagecoupling");
  }
  else
    config->start_imagecoupling=INT_MAX;
#endif
  if(config->firstyear-config->nspinup>config->lastyear)
  {
    if(verbose)
      fprintf(stderr,"ERROR105: First simulation year=%d greater than last simulation year=%d.\n",
              config->firstyear-config->nspinup,config->lastyear);
    return TRUE;
  }
  if(iskeydefined(file,"outputyear"))
  {
    fscanint2(file,&config->outputyear,"outputyear");
    if(config->outputyear>config->lastyear)
    {
      if(verbose)
        fprintf(stderr,"ERROR230: First year output is written=%d greater than last simulation year=%d.\n",
                config->outputyear,config->lastyear);
      return TRUE;
    }
    else if(config->outputyear<config->firstyear-config->nspinup)
    {
      if(verbose)
        fprintf(stderr,"ERROR230: First year output is written=%d less than first simulation year=%d.\n",
                config->outputyear,config->firstyear-config->nspinup);
      return TRUE;
    }
  }
  else
    config->outputyear=config->firstyear;
  fscanbool2(file,&config->from_restart,"restart");
  if(config->from_restart)
  {
    fscanname(file,name,"restart_filename");
    config->restart_filename=addpath(name,config->restartdir);
  }
  else
  {
    config->restart_filename=NULL;
    if(verbose && config->nspinup<soil_equil_year)
      fprintf(stderr,"WARNING031: Number of spinup years less than %d necessary for soil equilibration.\n",soil_equil_year);
  }
  if(iskeydefined(file,"checkpoint_filename"))
  {
    fscanname(file,name,"checkpoint_filename");
    config->checkpoint_restart_filename=addpath(name,config->restartdir);
  }
  else
    config->checkpoint_restart_filename=NULL;
  fscanbool2(file,&restart,"write_restart");
  if(restart)
  {
    fscanname(file,name,"write_restart_filename");
    config->write_restart_filename=addpath(name,config->restartdir);
    fscanint2(file,&config->restartyear,"restart_year");
    if(config->restartyear>config->lastyear || config->restartyear<config->firstyear-config->nspinup)
    {
      if(verbose)
        fprintf(stderr,"WARNING014: Restart year=%d not in simulation years, no restart file written.\n",config->restartyear);
      free(config->write_restart_filename);
      config->write_restart_filename=NULL;
    }
  }
  else
    config->write_restart_filename=NULL;
  return FALSE;
} /* of 'fscanconfig' */
