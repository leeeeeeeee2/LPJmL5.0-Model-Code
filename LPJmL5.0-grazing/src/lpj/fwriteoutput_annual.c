/**************************************************************************************/
/**                                                                                \n**/
/**      f  w  r  i  t  e  o  u  t  p  u  t  _  a  n  n  u  a  l  .  c             \n**/
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
#include "tree.h"

#define writeoutputvar(index,name) if(isopen(output,index))\
  {\
    count=0;\
    for(cell=0;cell<config->ngridcell;cell++)\
      if(!grid[cell].skip)\
        vec[count++]=(float)(grid[cell].output.name);\
    writeannual(output,index,vec,year,config);\
  }

static void writeannual(Outputfile *output,int index,float data[],int year,
                        const Config *config)
{
  int i;
#ifdef USE_MPI
  MPI_Status status;
#endif
  for(i=0;i<config->count;i++)
    data[i]=(float)(config->outnames[index].scale*data[i]+config->outnames[index].offset);
#ifdef USE_MPI
  switch(output->method)
  {
    case LPJ_MPI2:
      MPI_File_write_at(output->files[index].fp.mpi_file,
         (year-config->outputyear)*config->total+config->offset,data,config->count,
                        MPI_FLOAT,&status);
      break;
    case LPJ_GATHER:
      switch(output->files[index].fmt)
      {
        case RAW: case CLM:
          mpi_write(output->files[index].fp.file,data,MPI_FLOAT,config->total,
                    output->counts,output->offsets,config->rank,config->comm);
          break;
        case TXT:
          mpi_write_txt(output->files[index].fp.file,data,MPI_FLOAT,config->total,
                        output->counts,output->offsets,config->rank,config->comm);
          break;
        case CDF:
          mpi_write_netcdf(&output->files[index].fp.cdf,data,MPI_FLOAT,config->total,
                           output->files[index].oneyear ? NO_TIME : year-config->outputyear,
                           output->counts,output->offsets,config->rank,config->comm);
          break;
      }
      break;
    case LPJ_SOCKET:
      if(isroot(*config))
        writeint_socket(output->socket,&index,1);
      mpi_write_socket(output->socket,data,MPI_FLOAT,config->total,
                       output->counts,output->offsets,config->rank,config->comm);
      break;
  } /* of switch */
#else
  if(output->method==LPJ_FILES)
    switch(output->files[index].fmt)
    {
      case RAW: case CLM:
        if(fwrite(data,sizeof(float),config->count,output->files[index].fp.file)!=config->count)
          fprintf(stderr,"ERROR204: Error writing output: %s.\n",strerror(errno));
        break;
      case TXT:
        for(i=0;i<config->count-1;i++)
          fprintf(output->files[index].fp.file,"%g ",data[i]);
        fprintf(output->files[index].fp.file,"%g\n",data[config->count-1]);
        break;
      case CDF:
        write_float_netcdf(&output->files[index].fp.cdf,data,
                           output->files[index].oneyear ? NO_TIME : year-config->outputyear,
                           config->count);
        break;
    }
  else
  {
    writeint_socket(output->socket,&index,1);
    writefloat_socket(output->socket,data,config->count);
  }
#endif
} /* of 'writeannual' */

static void writeshortannual(Outputfile *output,int index,short data[],int year,
                             const Config *config)
{
#ifdef USE_MPI
  MPI_Status status;
  switch(output->method)
  {
    case LPJ_MPI2:
      MPI_File_write_at(output->files[index].fp.mpi_file,
         (year-config->outputyear)*config->total+config->offset,data,config->count,
                        MPI_SHORT,&status);
      break;
    case LPJ_GATHER:
      switch(output->files[index].fmt)
      {
        case RAW: case CLM:
          mpi_write(output->files[index].fp.file,data,MPI_SHORT,config->total,
                    output->counts,output->offsets,config->rank,config->comm);
          break;
        case TXT:
          mpi_write_txt(output->files[index].fp.file,data,MPI_SHORT,config->total,
                        output->counts,output->offsets,config->rank,config->comm);
          break;
        case CDF:
          mpi_write_netcdf(&output->files[index].fp.cdf,data,MPI_SHORT,config->total,
                           output->files[index].oneyear ? NO_TIME : year-config->outputyear,
                           output->counts,output->offsets,config->rank,config->comm);
          break;
      }
      break;
    case LPJ_SOCKET:
      if(isroot(*config))
        writeint_socket(output->socket,&index,1);
      mpi_write_socket(output->socket,data,MPI_SHORT,config->total,
                       output->counts,output->offsets,config->rank,config->comm);
      break;
  } /* of switch */
#else
  int i;
  if(output->method==LPJ_FILES)
    switch(output->files[index].fmt)
    {
      case RAW: case CLM:
        fwrite(data,sizeof(short),config->count,output->files[index].fp.file);
        break;
      case TXT:
        for(i=0;i<config->count-1;i++)
          fprintf(output->files[index].fp.file,"%d ",data[i]);
        fprintf(output->files[index].fp.file,"%d\n",data[config->count-1]);
        break;
      case CDF:
        write_short_netcdf(&output->files[index].fp.cdf,data,
                           output->files[index].oneyear ? NO_TIME : year-config->outputyear,
                           config->count);
        break;
    }
  else
  {
    writeint_socket(output->socket,&index,1);
    writeshort_socket(output->socket,data,config->count);
  }
#endif
} /* of 'writeshortannual' */

static void writeannualall(Outputfile *output,int index,float data[],int year,
                           const Config *config)
{
  int i;
#ifdef USE_MPI
  int *counts,*offsets;
  MPI_Status status;
#endif
  for(i=0;i<config->ngridcell;i++)
    data[i]=(float)(config->outnames[index].scale*data[i]+config->outnames[index].offset);
#ifdef USE_MPI
  switch(output->method)
  {
    case LPJ_MPI2:
      MPI_File_write_at(output->files[index].fp.mpi_file,
        (year-config->outputyear)*config->nall+config->offset,data,config->ngridcell,
                        MPI_FLOAT,&status);
      break;
    case LPJ_GATHER:
      counts=newvec(int,config->ntask);
      check(counts);
      offsets=newvec(int,config->ntask);
      check(offsets);
      getcounts(counts,offsets,config->nall,1,config->ntask);
      switch(output->files[index].fmt)
      {
        case RAW: case CLM:
          mpi_write(output->files[index].fp.file,data,MPI_FLOAT,config->nall,counts,
                    offsets,config->rank,config->comm);
          break;
        case TXT:
          mpi_write_txt(output->files[index].fp.file,data,MPI_FLOAT,config->nall,counts,
                        offsets,config->rank,config->comm);
          break;
        case CDF:
          mpi_write_netcdf(&output->files[index].fp.cdf,data,MPI_FLOAT,config->nall,
                           output->files[index].oneyear ? NO_TIME : year-config->outputyear,
                           counts,offsets,config->rank,config->comm);
          break;
      }
      free(counts);
      free(offsets);
      break;
    case LPJ_SOCKET:
      counts=newvec(int,config->ntask);
      check(counts);
      offsets=newvec(int,config->ntask);
      check(offsets);
      getcounts(counts,offsets,config->nall,1,config->ntask);
      if(isroot(*config))
        writeint_socket(output->socket,&index,1);
      mpi_write_socket(output->socket,data,MPI_FLOAT,config->nall,counts,
                       offsets,config->rank,config->comm);
      free(counts);
      free(offsets);
      break;
  } /* of switch */
#else
  if(output->method==LPJ_FILES)
    switch(output->files[index].fmt)
    {
      case RAW: case CLM:
        fwrite(data,sizeof(float),config->ngridcell,output->files[index].fp.file);
        break;
      case TXT:
        for(i=0;i<config->ngridcell-1;i++)
          fprintf(output->files[index].fp.file,"%g ",data[i]);
        fprintf(output->files[index].fp.file,"%g\n",data[config->ngridcell-1]);
        break;
      case CDF:
        write_float_netcdf(&output->files[index].fp.cdf,data,
                           output->files[index].oneyear ? NO_TIME : year-config->outputyear,
                           config->ngridcell);
        break;
    }
  else
  {
    writeint_socket(output->socket,&index,1);
    writefloat_socket(output->socket,data,config->ngridcell);
  }
#endif
} /* of 'writeannualall' */

void fwriteoutput_annual(Outputfile *output,  /**< output file array */
                         const Cell grid[],   /**< grid cell array */
                         int year,            /**< simulation year (AD) */
                         const Config *config /**< LPJ configuration */
                        )
{
  int count,s,p,cell,l;
  Stand *stand;
  Pft *pft;
  float *vec;
  short *svec;
  Real fracs;
  Real grassfrac=0;
  if(isopen(output,SEASONALITY))
  {
    count=0;
    svec=newvec(short,config->ngridcell);
    check(svec);
    for(cell=0;cell<config->ngridcell;cell++)
      if(!grid[cell].skip)
        svec[count++]=(short)(grid[cell].ml.seasonality_type);
    writeshortannual(output,SEASONALITY,svec,year,config);
    free(svec);
  }
  vec=newvec(float,config->ngridcell);
  check(vec);
  writeoutputvar(ALITFALLC,alittfall.carbon);
  writeoutputvar(ALITFALLN,alittfall.nitrogen);
  writeoutputvar(FIREC,fire.carbon);
  writeoutputvar(FIREN,fire.nitrogen);
  writeoutputvar(ABURNTAREA,aburntarea);
  writeoutputvar(FLUX_FIREWOOD,flux_firewood.carbon);
  writeoutputvar(FIREF,firef);
  if(isopen(output,VEGC))
  {
    count=0;
    for(cell=0;cell<config->ngridcell;cell++)
      if(!grid[cell].skip)
      {
        vec[count]=0;
        foreachstand(stand,s,grid[cell].standlist)
          /*if(stand->type->landusetype==NATURAL) */
            foreachpft(pft,p,&stand->pftlist)
               vec[count]+=(float)(vegc_sum(pft)*stand->frac);
        count++;
      }
    writeannual(output,VEGC,vec,year,config);
  }
  if(isopen(output,SOILC))
  {
    count=0;
    for(cell=0;cell<config->ngridcell;cell++)
      if(!grid[cell].skip)
      {
        vec[count]=0;
        foreachstand(stand,s,grid[cell].standlist)
        {
          for(p=0;p<stand->soil.litter.n;p++)
            vec[count]+=(float)(stand->soil.litter.item[p].bg.carbon*stand->frac);
          forrootsoillayer(l)
            vec[count]+=(float)((stand->soil.pool[l].slow.carbon+stand->soil.pool[l].fast.carbon)*stand->frac);
        }
        count++;
      }
    writeannual(output,SOILC,vec,year,config);
  }
  if(isopen(output,SOILC_SLOW))
  {
    count=0;
    for(cell=0;cell<config->ngridcell;cell++)
      if(!grid[cell].skip)
      {
        vec[count]=0;
        foreachstand(stand,s,grid[cell].standlist)
        {
          forrootsoillayer(l)
            vec[count]+=(float)((stand->soil.pool[l].slow.carbon)*stand->frac);
          vec[count]+=(float)(stand->soil.YEDOMA*stand->frac);
        }
        count++;
      }
    writeannual(output,SOILC_SLOW,vec,year,config);
  }
  if (output->files[SOILC_AGR].isopen)
  {
    count = 0;
    for (cell = 0; cell<config->ngridcell; cell++)
      if (!grid[cell].skip)
      {
        vec[count] = 0;
        fracs = 0;
        foreachstand(stand, s, grid[cell].standlist)
        {
          if (stand->type->landusetype == SETASIDE_RF || stand->type->landusetype == SETASIDE_IR ||
            stand->type->landusetype == AGRICULTURE)
          {
            fracs += stand->frac;
            for (p = 0; p<stand->soil.litter.n; p++)
              vec[count] += (float)(stand->soil.litter.item[p].bg.carbon*stand->frac);
            forrootsoillayer(l)
              vec[count] += (float)((stand->soil.pool[l].slow.carbon + stand->soil.pool[l].fast.carbon)*stand->frac);
          }
        }
        if (fracs>0) vec[count] /= fracs; /* as fracs don't add up to 1 here */
        count++;
      }
    writeannual(output, SOILC_AGR, vec, year, config);
  }

  if(isopen(output,LITC))
  {
    count=0;
    for(cell=0;cell<config->ngridcell;cell++)
      if(!grid[cell].skip)
      {
        vec[count]=0;
        foreachstand(stand,s,grid[cell].standlist)
          /* if(stand->type->landusetype==NATURAL) */
            vec[count]+=(float)((litter_ag_sum(&stand->soil.litter)+litter_agsub_sum(&stand->soil.litter))*stand->frac);
        count++;
      }
      writeannual(output,LITC,vec,year,config);
  }
  if(output->files[LITC_ALL].isopen)
  {
    count=0;
    for(cell=0;cell<config->ngridcell;cell++)
      if(!grid[cell].skip)
      {
        vec[count]=0;
        foreachstand(stand,s,grid[cell].standlist)
        {
          vec[count]+=(float)((litter_ag_sum(&stand->soil.litter)+litter_agsub_sum(&stand->soil.litter))*stand->frac);
          for(p=0;p<stand->soil.litter.n;p++)
            vec[count]+=(float)(stand->soil.litter.item[p].bg.carbon*stand->frac);
        }
        count++;
      }
      writeannual(output,LITC_ALL,vec,year,config);
  }
  if(output->files[LITC_AG].isopen)
  {
    count=0;
    for(cell=0;cell<config->ngridcell;cell++)
      if(!grid[cell].skip)
      {
        vec[count]=0;
        foreachstand(stand,s,grid[cell].standlist)
          vec[count]+=(float)(litter_ag_sum(&stand->soil.litter)*stand->frac);
        count++;
      }
      writeannual(output,LITC_AG,vec,year,config);
  }
  if (output->files[LITC_AGR].isopen)
  {
    count = 0;
    for (cell = 0; cell<config->ngridcell; cell++)
      if (!grid[cell].skip)
      {
        vec[count] = 0;
        fracs = 0;
        foreachstand(stand, s, grid[cell].standlist)
        {
          if (stand->type->landusetype == SETASIDE_RF || stand->type->landusetype == SETASIDE_IR || stand->type->landusetype == AGRICULTURE)
          {
            fracs += stand->frac;
            vec[count] += (float)((litter_ag_sum(&stand->soil.litter) + litter_agsub_sum(&stand->soil.litter))*stand->frac);
          }
        }
        if (fracs>0) vec[count] /= fracs; /* as fracs don't add up to 1 here */
        count++;
      }
    writeannual(output, LITC_AGR, vec, year, config);
  }

  if(isopen(output,VEGN))
  {
    count=0;
    for(cell=0;cell<config->ngridcell;cell++)
      if(!grid[cell].skip)
      {
        vec[count]=0;
        foreachstand(stand,s,grid[cell].standlist)
          /*if(stand->type->landusetype==NATURAL) */
            foreachpft(pft,p,&stand->pftlist)
               vec[count]+=(float)(vegn_sum(pft)*stand->frac);
        count++;
      }
    writeannual(output,VEGN,vec,year,config);
  }
  if(isopen(output,SOILN))
  {
    count=0;
    for(cell=0;cell<config->ngridcell;cell++)
      if(!grid[cell].skip)
      {
        vec[count]=0;
        foreachstand(stand,s,grid[cell].standlist)
        {
          for(p=0;p<stand->soil.litter.n;p++)
            vec[count]+=(float)(stand->soil.litter.item[p].bg.nitrogen*stand->frac);
          forrootsoillayer(l)
            vec[count]+=(float)((stand->soil.pool[l].slow.nitrogen+stand->soil.pool[l].fast.nitrogen)*stand->frac);
        }
        count++;
      }
    writeannual(output,SOILN,vec,year,config);
  }
  if(isopen(output,SOILN_SLOW))
  {
    count=0;
    for(cell=0;cell<config->ngridcell;cell++)
      if(!grid[cell].skip)
      {
        vec[count]=0;
        foreachstand(stand,s,grid[cell].standlist)
        {
          forrootsoillayer(l)
            vec[count]+=(float)((stand->soil.pool[l].slow.nitrogen)*stand->frac);
          /*vec[count]+=(float)(stand->soil.YEDOMA*stand->frac);*/
        }
        count++;
      }
    writeannual(output,SOILN_SLOW,vec,year,config);
  }
  if(isopen(output,LITN))
  {
    count=0;
    for(cell=0;cell<config->ngridcell;cell++)
      if(!grid[cell].skip)
      {
        vec[count]=0;
        foreachstand(stand,s,grid[cell].standlist)
          /* if(stand->type->landusetype==NATURAL) */
            vec[count]+=(float)((litter_ag_sum_n(&stand->soil.litter)+litter_agsub_sum_n(&stand->soil.litter))*stand->frac);
        count++;
      }
      writeannual(output,LITN,vec,year,config);
  }
  if(isopen(output,SOILNO3))
  {
    count=0;
    for(cell=0;cell<config->ngridcell;cell++)
      if(!grid[cell].skip)
      {
        vec[count]=0;
        foreachstand(stand,s,grid[cell].standlist)
        {
          forrootsoillayer(l)
		  {
            if(stand->soil.mean_maxthaw>=layerbound[l])
        	  vec[count]+=(float)((stand->soil.NO3[l])*stand->frac);
          /*vec[count]+=(float)(stand->soil.YEDOMA*stand->frac);*/
		  }
        }
        count++;
      }
    writeannual(output,SOILNO3,vec,year,config);
  }
  if(isopen(output,SOILNH4))
  {
    count=0;
    for(cell=0;cell<config->ngridcell;cell++)
      if(!grid[cell].skip)
      {
        vec[count]=0;
        foreachstand(stand,s,grid[cell].standlist)
        {
          forrootsoillayer(l)
		  {
		    if(stand->soil.mean_maxthaw>=layerbound[l])
              vec[count]+=(float)((stand->soil.NH4[l])*stand->frac);
            /*vec[count]+=(float)(stand->soil.YEDOMA*stand->frac);*/
		  }
        }
        count++;
      }
      writeannual(output,SOILNH4,vec,year,config);
  }
  if(isopen(output,MAXTHAW_DEPTH))
  {
    count=0;
    for(cell=0;cell<config->ngridcell;cell++)
      if(!grid[cell].skip)
      {
        vec[count]=0;
        foreachstand(stand,s,grid[cell].standlist)
          vec[count]+=(float)(stand->soil.maxthaw_depth*stand->frac*(1.0/(1-stand->cell->lakefrac-stand->cell->ml.reservoirfrac)));
        count++;
      }
      writeannual(output,MAXTHAW_DEPTH,vec,year,config);
  }
  writeoutputvar(RUNOFF_SURF,runoff_surf);
  writeoutputvar(RUNOFF_LAT,runoff_lat);
  writeoutputvar(FLUX_ESTABC,flux_estab.carbon);
  writeoutputvar(FLUX_ESTABN,flux_estab.nitrogen);
  writeoutputvar(HARVESTC,flux_harvest.carbon);
  writeoutputvar(HARVESTN,flux_harvest.nitrogen);
  writeoutputvar(RHARVEST_BURNTC,flux_rharvest_burnt.carbon);
  writeoutputvar(RHARVEST_BURNTN,flux_rharvest_burnt.nitrogen);
  writeoutputvar(RHARVEST_BURNT_IN_FIELDC,flux_rharvest_burnt_in_field.carbon);
  writeoutputvar(RHARVEST_BURNT_IN_FIELDN,flux_rharvest_burnt_in_field.nitrogen);
  writeoutputvar(ANPP,anpp);
  writeoutputvar(ANPP_AGR,anpp_agr);
  writeoutputvar(ARH,arh);
  writeoutputvar(ARH_AGR,arh_agr);
  if(isopen(output,MG_VEGC))
  {
    count=0;
    for(cell=0;cell<config->ngridcell;cell++)
      if(!grid[cell].skip)
      {
        vec[count]=0;
        foreachstand(stand,s,grid[cell].standlist)
          if(stand->type->landusetype!=NATURAL)
            foreachpft(pft,p,&stand->pftlist)
              vec[count]+=(float)(vegc_sum(pft)*stand->frac);
        count++;
      }
      writeannual(output,MG_VEGC,vec,year,config);
  }
  if(isopen(output,MG_SOILC))
  {
      count=0;
      for(cell=0;cell<config->ngridcell;cell++)
        if(!grid[cell].skip)
        {
          vec[count]=0;
          foreachstand(stand,s,grid[cell].standlist)
          {
            if(stand->type->landusetype!=NATURAL)
            {
              for(p=0;p<stand->soil.litter.n;p++)
                vec[count]+=(float)(stand->soil.litter.item[p].bg.carbon*stand->frac);
              forrootsoillayer(l)
                vec[count]+=(float)((stand->soil.pool[l].slow.carbon+stand->soil.pool[l].fast.carbon)*stand->frac);
            }
          }
          count++;
        }
      writeannual(output,MG_SOILC,vec,year,config);
  }
  if(isopen(output,MG_LITC))
  {
    count=0;
    for(cell=0;cell<config->ngridcell;cell++)
      if(!grid[cell].skip)
      {
        vec[count]=0;
        foreachstand(stand,s,grid[cell].standlist)
          if(stand->type->landusetype!=NATURAL)
            vec[count]+=(float)((litter_ag_sum(&stand->soil.litter)+litter_agsub_sum(&stand->soil.litter))*stand->frac);
        count++;
      }
    writeannual(output,MG_LITC,vec,year,config);
  }
  if(isopen(output,APREC))
  {
    count=0;
    for(cell=0;cell<config->ngridcell;cell++)
      if(!grid[cell].skip)
        vec[count++]=(float)grid[cell].balance.aprec;
    writeannual(output,APREC,vec,year,config);
  }
  writeoutputvar(INPUT_LAKE,input_lake*1e-9);
  if(isopen(output,ADISCHARGE))
  {
    for(cell=0;cell<config->ngridcell;cell++)
      vec[cell]=(float)(grid[cell].output.adischarge*1e-9);
    writeannualall(output,ADISCHARGE,vec,year,config);
  }
  writeoutputvar(DEFOREST_EMIS,deforest_emissions.carbon);
  writeoutputvar(TRAD_BIOFUEL,trad_biofuel);
  writeoutputvar(AIRRIG,airrig);
  writeoutputvar(FBURN,fburn);
  writeoutputvar(FTIMBER,ftimber);
  writeoutputvar(TIMBER_HARVESTC,timber_harvest.carbon);
#ifdef IMAGE
  writeoutputvar(PRODUCT_POOL_FAST,product_pool_fast);
  writeoutputvar(PRODUCT_POOL_SLOW,product_pool_slow);
  writeoutputvar(PROD_TURNOVER,prod_turnover);
#else
  writeoutputvar(PRODUCT_POOL_FAST,product_pool_fast.carbon);
  writeoutputvar(PRODUCT_POOL_SLOW,product_pool_slow.carbon);
  writeoutputvar(PROD_TURNOVER,prod_turnover.carbon);
#endif
  if(isopen(output,AFRAC_WD_UNSUST))
  {
    count=0;
    for(cell=0;cell<config->ngridcell;cell++)
      if(!grid[cell].skip)
         vec[count++]=(float)(grid[cell].output.awd_unsustainable/((grid[cell].output.airrig+
            grid[cell].output.aconv_loss_evap + grid[cell].output.aconv_loss_drain)*grid[cell].coord.area));
    writeannual(output,AFRAC_WD_UNSUST,vec,year,config);
  }
  writeoutputvar(ACONV_LOSS_EVAP,aconv_loss_evap);
  writeoutputvar(ACONV_LOSS_DRAIN,aconv_loss_drain);
  writeoutputvar(AWATERUSE_HIL,awateruse_hil);
  if(isopen(output,AGB))
  {
    count=0;
    for(cell=0;cell<config->ngridcell;cell++)
      if(!grid[cell].skip)
      {
        vec[count]=0;
        foreachstand(stand,s,grid[cell].standlist)
          foreachpft(pft,p,&stand->pftlist)
            vec[count]+=(float)(agb(pft)*stand->frac);
        count++;
      }
    writeannual(output,AGB,vec,year,config);
  }
  if(isopen(output,AGB_TREE))
  {
    count=0;
    for(cell=0;cell<config->ngridcell;cell++)
      if(!grid[cell].skip)
      {
        vec[count]=0;
        foreachstand(stand,s,grid[cell].standlist)
          foreachpft(pft,p,&stand->pftlist)
            if(istree(pft))
              vec[count]+=(float)(agb_tree(pft)*stand->frac);
        count++;
      }
    writeannual(output,AGB_TREE,vec,year,config);
  }
  if(isopen(output,MGRASS_SOILC))
  {
    count=0;
    for(cell=0;cell<config->ngridcell;cell++)
      if(!grid[cell].skip)
      {
        vec[count]=0;
        grassfrac=0;
        foreachstand(stand,s,grid[cell].standlist)
        {
          if(stand->type->landusetype==GRASSLAND){
            for(p=0;p<stand->soil.litter.n;p++)
              vec[count]+=(float)(stand->soil.litter.item[p].bg.carbon*stand->frac);
            forrootsoillayer(l)
              vec[count]+=(float)((stand->soil.pool[l].slow.carbon+stand->soil.pool[l].fast.carbon)*stand->frac);
            grassfrac+=stand->frac;
          }
        }
        if(grassfrac>epsilon)
          vec[count]/=grassfrac;
        count++;
      }
    writeannual(output,MGRASS_SOILC,vec,year,config);
  }
  if(isopen(output,MGRASS_LITC))
  {
    count=0;
    for(cell=0;cell<config->ngridcell;cell++)
      if(!grid[cell].skip)
      {
        vec[count]=0;
        grassfrac=0;
        foreachstand(stand,s,grid[cell].standlist)
        {
          if(stand->type->landusetype==GRASSLAND){
            vec[count]+=(float)(litter_ag_sum(&stand->soil.litter)*stand->frac);
            grassfrac+=stand->frac;
          }
        }
        if(grassfrac>epsilon)
          vec[count]/=grassfrac;
        count++;
      }
    writeannual(output,MGRASS_LITC,vec,year,config);
  }
  if(isopen(output,MGRASS_SOILN))
  {
    count=0;
    for(cell=0;cell<config->ngridcell;cell++)
      if(!grid[cell].skip)
      {
        vec[count]=0;
        grassfrac=0;
        foreachstand(stand,s,grid[cell].standlist)
        {
          if(stand->type->landusetype==GRASSLAND){
            for(p=0;p<stand->soil.litter.n;p++)
              vec[count]+=(float)(stand->soil.litter.item[p].bg.nitrogen*stand->frac);
            forrootsoillayer(l)
              vec[count]+=(float)((stand->soil.pool[l].slow.nitrogen+stand->soil.pool[l].fast.nitrogen)*stand->frac);
            grassfrac+=stand->frac;
          }
        }
        if(grassfrac>epsilon)
          vec[count]/=grassfrac;
        count++;
      }
    writeannual(output,MGRASS_SOILN,vec,year,config);
  }
  if(isopen(output,MGRASS_LITN))
  {
    count=0;
    for(cell=0;cell<config->ngridcell;cell++)
      if(!grid[cell].skip)
      {
        vec[count]=0;
        grassfrac=0;
        foreachstand(stand,s,grid[cell].standlist)
        {
          if(stand->type->landusetype==GRASSLAND){
            vec[count]+=(float)(litter_ag_sum_n(&stand->soil.litter)*stand->frac);
            grassfrac+=stand->frac;
          }
        }
        if(grassfrac>epsilon)
          vec[count]/=grassfrac;
        count++;
      }
    writeannual(output,MGRASS_LITN,vec,year,config);
  }

  writeoutputvar(NEGC_FLUXES,neg_fluxes.carbon);
  writeoutputvar(NEGN_FLUXES,neg_fluxes.nitrogen);
  writeoutputvar(MEANVEGCMANGRASS,mean_vegc_mangrass);
  writeoutputvar(ABNF_AGR,abnf_agr);
  writeoutputvar(ANFERT_AGR,anfert_agr);
  writeoutputvar(ANFERT_MGRASS,anfert_mgrass);
  writeoutputvar(ANMANURE_AGR,anmanure_agr);
  writeoutputvar(ANDEPO_AGR,andepo_agr);
  writeoutputvar(ANDEPO_MGRASS,andepo_mgrass);
  writeoutputvar(ANMINERALIZATION_AGR,anmineralization_agr);
  writeoutputvar(ANMINERALIZATION_MGRASS,anmineralization_mgrass);
  writeoutputvar(ANIMMOBILIZATION_AGR,animmobilization_agr);
  writeoutputvar(ANIMMOBILIZATION_MGRASS,animmobilization_mgrass);
  writeoutputvar(ANUPTAKE_AGR,anuptake_agr);
  writeoutputvar(ANUPTAKE_MGRASS,anuptake_mgrass);
  writeoutputvar(ANLEACHING_AGR,anleaching_agr);
  writeoutputvar(AN2O_DENIT_AGR,an2o_denit_agr);
  writeoutputvar(AN2O_NIT_AGR,an2o_nit_agr);
  writeoutputvar(ANH3_AGR,anh3_agr);
  writeoutputvar(AN2_AGR,an2_agr);
  writeoutputvar(ALITFALLN_AGR,alitfalln_agr);
  writeoutputvar(ALITFALLN_MGRASS,alitfalln_mgrass);
  writeoutputvar(AHARVESTN_AGR,aharvestn_agr);
  writeoutputvar(AHARVESTN_MGRASS,aharvestn_mgrass);
  writeoutputvar(ASEEDN_AGR,aseedn_agr);
  writeoutputvar(ASEEDN_MGRASS,aseedn_mgrass);
  writeoutputvar(ADELTA_NORG_SOIL_AGR,adelta_norg_soil_agr);
  writeoutputvar(ADELTA_NMIN_SOIL_AGR,adelta_nmin_soil_agr);
  writeoutputvar(ADELTA_NVEG_SOIL_AGR,adelta_nveg_soil_agr);
  writeoutputvar(ADELTA_NORG_SOIL_MGRASS,adelta_norg_soil_mgrass);
  writeoutputvar(ADELTA_NMIN_SOIL_MGRASS,adelta_nmin_soil_mgrass);
  writeoutputvar(ADELTA_NVEG_SOIL_MGRASS,adelta_nveg_soil_mgrass);
  writeoutputvar(CELLFRAC_AGR,cellfrac_agr);
  writeoutputvar(ALITFALLC_WOOD,alittfall_wood.carbon);
  writeoutputvar(ALITFALLN_WOOD,alittfall_wood.nitrogen);
  writeoutputvar(DECAY_WOOD_AGR,decay_wood_agr);
  writeoutputvar(DECAY_WOOD_NV,decay_wood_nv);
  writeoutputvar(DECAY_LEAF_AGR,decay_leaf_agr);
  writeoutputvar(DECAY_LEAF_NV,decay_leaf_nv);
  writeoutputvar(ALITBURNC,alitburnc);
  writeoutputvar(ALITBURNC_WOOD,alitburnc_wood);
  writeoutputvar(MGRASS_N2O_DENIT,mgrass_n2o_denit);
  writeoutputvar(MGRASS_N2O_NIT,mgrass_n2o_nit);
  writeoutputvar(MGRASS_N2_EMIS,mgrass_n2_emis);
  writeoutputvar(MGRASS_LEACHING,mgrass_leaching);
  writeoutputvar(MGRASS_BNF,mgrass_bnf);
  writeoutputvar(MGRASS_VOLATILIZATION,mgrass_volatilization);
  writeoutputvar(YIELDC_MGRASS,yieldc_mgrass);
  writeoutputvar(MANUREC_MGRASS,manurec_mgrass);
  writeoutputvar(URINEC_MGRASS,urinec_mgrass);
  writeoutputvar(METHANEC_MGRASS,methanec_mgrass);
  writeoutputvar(RESPC_MGRASS,respc_mgrass);
  writeoutputvar(YIELDN_MGRASS,yieldn_mgrass);
  writeoutputvar(MANUREN_MGRASS,manuren_mgrass);
  writeoutputvar(URINEN_MGRASS,urinen_mgrass);
  writeoutputvar(UPTAKEC_MGRASS,uptakec_mgrass);
  writeoutputvar(UPTAKEN_MGRASS,uptaken_mgrass);
  writeoutputvar(DEFICIT_LSU_MP,deficit_lsu_mp);
  writeoutputvar(DEFICIT_LSU_NE,deficit_lsu_ne);
  writeoutputvar(MGRASS_GPP,mgrass_gpp);
  writeoutputvar(MGRASS_RA,mgrass_ra);
  writeoutputvar(MGRASS_RH,mgrass_rh);
  writeoutputvar(MGRASS_DELTAC,mgrass_deltac);
  writeoutputvar(MGRASS_SEEDC,mgrass_seedc);
  writeoutputvar(MGRASS_AGBMEAN,mgrass_agbmean);

  free(vec);
} /* of 'fwriteoutput_annual' */
