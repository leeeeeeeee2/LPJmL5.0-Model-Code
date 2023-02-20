/**************************************************************************************/
/**                                                                                \n**/
/**                           l  a  n  d  u  s  e  .  c                            \n**/
/**                                                                                \n**/
/**     C implementation of LPJmL                                                  \n**/
/**                                                                                \n**/
/**     Function initializes landuse datatype                                      \n**/
/**     - opens the landuse input file (see also cfts26_lu2clm.c)                  \n**/
/**     - sets the landuse variables (see also landuse.h)                          \n**/
/**     - order of landuse input data:                                             \n**/
/**        0-10   CFT                                                              \n**/
/**          11   OTHERS                                                           \n**/
/**          12   PASTURES                                                         \n**/
/**          13   BIOMASS_GRASS                                                    \n**/
/**          14   BIOMASS_TREE                                                     \n**/
/**       15-25   CFT_irr                                                          \n**/
/**          26   others_irr                                                       \n**/
/**          27   PASTURE_irr                                                      \n**/
/**          28   BIOMASS_GRASS IRR                                                \n**/
/**          29   BIOMASS_TREE IRR                                                 \n**/
/**     - called in iterate()                                                      \n**/
/**     - reads every year the fractions of the bands for all cells from           \n**/
/**       the input file                                                           \n**/
/**     - checks if sum of fraction is not greater 1                               \n**/
/**       -> if sum of fraction is greater 1: subtraction from fraction            \n**/
/**          of managed grass if possible                                          \n**/
/**       -> else fail incorrect input file                                        \n**/
/**                                                                                \n**/
/**                                                                                \n**/
/** (C) Potsdam Institute for Climate Impact Research (PIK), see COPYRIGHT file    \n**/
/** authors, and contributors see AUTHORS file                                     \n**/
/** This file is part of LPJmL and licensed under GNU AGPL Version 3               \n**/
/** or later. See LICENSE file or go to http://www.gnu.org/licenses/               \n**/
/** Contact: https://github.com/PIK-LPJmL/LPJmL                                    \n**/
/**                                                                                \n**/
/**************************************************************************************/

#include "lpj.h"
#include "natural.h"
#include "agriculture.h"
#include "grassland.h"

/* define a tiny fraction for allcrops that is always at least 10x epsilon */

Real tinyfrac=max(epsilon*10,1e-6);

struct landuse
{
  Bool intercrop;      /**< intercropping possible (TRUE/FALSE) */
  Bool allcrops;       /**< all crops establish (TRUE/FALSE) */
  Bool onlycrops;      /* scale crop/grass shares to add to one (no natural vegetation) */
  Bool onlygrass;      /**< sumulate only managed grasslands (TRUE/FALSE) */
  int nbands;          /**< number of data items per cell */
  int nbands_sdate;    /**< number of data items per cell for sowing dates */
  Climatefile landuse; /**< file pointer */
  Climatefile fertilizer_nr; /**< file pointer to nitrogen fertilizer file */
  Climatefile manure_nr; /* file pointer to manure fertilizer file */
  Climatefile with_tillage; /* file pointer to tillage file */
  Climatefile residue_on_field; /* file pointer to residue extraction file */
  Climatefile sdate;   /**< file pointer to prescribed sdates */
  Climatefile crop_phu;   /**< file pointer to prescribed crop phus */
};                     /**< definition of opaque datatype Landuse */

Landuse initlanduse(int ncft,            /**< number of crop PFTs */
  const Config *config /**< LPJ configuration */
)                     /** \return allocated landuse or NULL */
{
  Header header;
  String headername;
  int version;
  Landuse landuse;
  size_t offset;
  landuse=new(struct landuse);
  if(landuse==NULL)
  {
    printallocerr("landuse");
    return NULL;
  }
  landuse->allcrops=(config->withlanduse==ALL_CROPS);
  landuse->onlycrops=(config->withlanduse==ONLY_CROPS);
  landuse->onlygrass=(config->withlanduse==ONLY_GRASS);
  landuse->landuse.fmt=config->landuse_filename.fmt;
  if(config->landuse_filename.fmt==CDF)
  {
    if(opendata_netcdf(&landuse->landuse,config->landuse_filename.name,config->landuse_filename.var,NULL,config))
    {
      free(landuse);
      return NULL;
    }
    landuse->nbands=landuse->landuse.var_len;
    if(landuse->landuse.var_len!=2*(ncft+NGRASS+NBIOMASSTYPE)&&landuse->landuse.var_len!=4*(ncft+NGRASS+NBIOMASSTYPE))
    {
      if(landuse->landuse.var_len!=2*(ncft+NGRASS))
      {
        closeclimate_netcdf(&landuse->landuse,isroot(*config));
        if(isroot(*config))
          fprintf(stderr,
            "ERROR147: Invalid number of bands=%d in landuse data file.\n",
            (int)landuse->landuse.var_len);
        free(landuse);
        return NULL;
      }
      else
      {
        if(isroot(*config))
          fputs("WARNING022: No landuse for biomass defined.\n",stderr);
      }
    }
    if(landuse->landuse.var_len!=4*(ncft+NGRASS+NBIOMASSTYPE))
    {
      if(isroot(*config))
        fputs("WARNING024: Land-use input does not include irrigation systems, suboptimal country values are used.\n",stderr);
    }
  }
  else
  {
    if((landuse->landuse.file=openinputfile(&header,&landuse->landuse.swap,
      &config->landuse_filename,
      headername,
      &version,&offset,config))==NULL)
    {
      free(landuse);
      return NULL;
    }
    if(config->landuse_filename.fmt==RAW)
    {
      header.nbands=2*(ncft+NGRASS+NBIOMASSTYPE);
      landuse->landuse.datatype=LPJ_SHORT;
      landuse->landuse.offset=config->startgrid*header.nbands*sizeof(short);
    }
    else
    {
      if(header.nbands!=2*(ncft+NGRASS+NBIOMASSTYPE)&&header.nbands!=4*(ncft+NGRASS+NBIOMASSTYPE))
      {
        if(header.nbands!=2*(ncft+NGRASS))
        {
          fclose(landuse->landuse.file);
          if(isroot(*config))
            fprintf(stderr,
              "ERROR147: Invalid number of bands=%d in landuse data file.\n",
              header.nbands);
          free(landuse);
          return NULL;
        }
        else
        {
          if(isroot(*config))
            fputs("WARNING022: No landuse for biomass defined.\n",stderr);
        }
      }
      landuse->landuse.datatype=header.datatype;
      if(header.nbands!=4*(ncft+NGRASS+NBIOMASSTYPE))
      {
        if(isroot(*config))
          fputs("WARNING024: Land-use input does not include irrigation systems, suboptimal country values are used.\n",stderr);
      }
      landuse->landuse.offset=(config->startgrid-header.firstcell)*header.nbands*typesizes[header.datatype]+headersize(headername,version)+offset;
    }
    landuse->landuse.firstyear=header.firstyear;
    landuse->landuse.nyear=header.nyear;
    landuse->landuse.size=(long long)header.ncell*(long long)header.nbands*typesizes[landuse->landuse.datatype];
    landuse->landuse.n=config->ngridcell*header.nbands;
    landuse->nbands=header.nbands;
    landuse->landuse.scalar=(version==1) ? 0.001 : header.scalar;
  }

    /* Multiple-years PRESCRIBED_SDATE */
  if(config->sdate_option==PRESCRIBED_SDATE)
  {
    /* read sdate input metadata */
    landuse->sdate.fmt=config->sdate_filename.fmt;
    if(config->sdate_filename.fmt==CDF)
    {
      if(opendata_netcdf(&landuse->sdate,config->sdate_filename.name,config->sdate_filename.var,NULL,config))
      {
        if(landuse->landuse.fmt==CDF)
          closeclimate_netcdf(&landuse->landuse,isroot(*config));
        else
          fclose(landuse->landuse.file);
        free(landuse);
        return NULL;
      }
      if(landuse->sdate.var_len!=2*ncft)
      {
        closeclimate_netcdf(&landuse->sdate,isroot(*config));
        if(isroot(*config))
          fprintf(stderr,
            "ERROR147: Invalid number of bands=%d in sowing date Nr data file.\n",
            (int)landuse->sdate.var_len);
        if(landuse->landuse.fmt==CDF)
          closeclimate_netcdf(&landuse->landuse,isroot(*config));
        else
          fclose(landuse->landuse.file);
        free(landuse);
        return NULL;
      }
    }
    else
    {
      if((landuse->sdate.file=openinputfile(&header,&landuse->sdate.swap,
        &config->sdate_filename,headername,
        &version,&offset,config))==NULL)
      {
        if(landuse->landuse.fmt==CDF)
          closeclimate_netcdf(&landuse->landuse,isroot(*config));
        else
          fclose(landuse->landuse.file);
        free(landuse);
        return NULL;
      }
      if(config->sdate_filename.fmt==RAW)
      {
        header.nbands=2*ncft;
        landuse->sdate.offset=(long long)config->startgrid*header.nbands*sizeof(short);
      }
      else
      {
        if(header.nbands!=2*ncft)
        {
          if(landuse->landuse.fmt==CDF)
            closeclimate_netcdf(&landuse->landuse,isroot(*config));
          else
            fclose(landuse->landuse.file);
          fclose(landuse->sdate.file);
          if(isroot(*config))
            fprintf(stderr,
              "ERROR147: Invalid number of bands=%d in sowing date data file.\n",
              header.nbands);
          free(landuse);
          return(NULL);
        }
        landuse->sdate.datatype=header.datatype;
        landuse->sdate.offset=((long long)config->startgrid-(long long)header.firstcell)*header.nbands*typesizes[landuse->sdate.datatype]+headersize(headername,version)+offset;
      }
      landuse->sdate.firstyear=header.firstyear;
      landuse->sdate.nyear=header.nyear;
      landuse->sdate.size=(long long)header.ncell*(long long)header.nbands*typesizes[landuse->sdate.datatype];
      landuse->sdate.n=config->ngridcell*header.nbands;
      landuse->nbands_sdate=header.nbands;
      landuse->sdate.scalar=header.scalar;
    }
  }
  else
  {
    landuse->sdate.file=NULL;
  } /* End sdate */

    /* Multiple-years PRESCRIBED_CROP_PHU */
  if(config->crop_phu_option==PRESCRIBED_CROP_PHU)
  {
    /* read sdate input metadata */
    landuse->crop_phu.fmt=config->crop_phu_filename.fmt;
    if(config->crop_phu_filename.fmt==CDF)
    {
      if(opendata_netcdf(&landuse->crop_phu,config->crop_phu_filename.name,config->crop_phu_filename.var,NULL,config))
      {
        if(landuse->landuse.fmt==CDF)
          closeclimate_netcdf(&landuse->landuse,isroot(*config));
        else
          fclose(landuse->landuse.file);
        if(landuse->sdate.file!=NULL)
        {
          if(landuse->sdate.fmt==CDF)
            closeclimate_netcdf(&landuse->sdate,isroot(*config));
          else
            fclose(landuse->sdate.file);
        }
        free(landuse);
        return NULL;
      }

      if(landuse->crop_phu.var_len!=2*ncft)
      {
        closeclimate_netcdf(&landuse->crop_phu,isroot(*config));
        if(isroot(*config))
          fprintf(stderr,
            "ERROR147: Invalid number of bands=%d in sowing date Nr data file.\n",
            (int)landuse->crop_phu.var_len);
        if(landuse->landuse.fmt==CDF)
          closeclimate_netcdf(&landuse->landuse,isroot(*config));
        else
          fclose(landuse->landuse.file);
        if(landuse->sdate.file!=NULL)
        {
          if(landuse->sdate.fmt==CDF)
            closeclimate_netcdf(&landuse->sdate,isroot(*config));
          else
            fclose(landuse->sdate.file);
        }
        free(landuse);
        return NULL;
      }
    }
    else
    {
      if((landuse->crop_phu.file=openinputfile(&header,&landuse->crop_phu.swap,
        &config->crop_phu_filename,headername,
        &version,&offset,config))==NULL)
      {
        if(landuse->landuse.fmt==CDF)
          closeclimate_netcdf(&landuse->landuse,isroot(*config));
        else
          fclose(landuse->landuse.file);
        if(landuse->sdate.file!=NULL)
        {
          if(landuse->sdate.fmt==CDF)
            closeclimate_netcdf(&landuse->sdate,isroot(*config));
          else
            fclose(landuse->sdate.file);
        }
        free(landuse);
        return NULL;
      }
      if(config->crop_phu_filename.fmt==RAW)
      {
        header.nbands=2*ncft;
        landuse->crop_phu.offset=(long long)config->startgrid*header.nbands*sizeof(short);
      }
      else
      {
        if(header.nbands!=2*ncft)
        {
          if(landuse->landuse.fmt==CDF)
            closeclimate_netcdf(&landuse->landuse,isroot(*config));
          else
            fclose(landuse->landuse.file);
          if(landuse->sdate.file!=NULL)
          {
            if(landuse->sdate.fmt==CDF)
              closeclimate_netcdf(&landuse->sdate,isroot(*config));
            else
              fclose(landuse->sdate.file);
          }
          fclose(landuse->crop_phu.file);
          if(isroot(*config))
            fprintf(stderr,
              "ERROR147: Invalid number of bands=%d in crop phu data file.\n",
              header.nbands);
          free(landuse);
          return(NULL);
        }
        landuse->crop_phu.datatype=header.datatype;
        landuse->crop_phu.offset=((long long)config->startgrid-(long long)header.firstcell)*header.nbands*typesizes[landuse->crop_phu.datatype]+headersize(headername,version)+offset;
      }
      landuse->crop_phu.firstyear=header.firstyear;
      landuse->crop_phu.nyear=header.nyear;
      landuse->crop_phu.size=(long long)header.ncell*(long long)header.nbands*typesizes[landuse->crop_phu.datatype];
      landuse->crop_phu.n=config->ngridcell*header.nbands;
      landuse->crop_phu.var_len=header.nbands;
      landuse->crop_phu.scalar=header.scalar;
    }
  }
  else
  {
    landuse->crop_phu.file=NULL;
  } /* End crop_phu */


  if(config->with_nitrogen&&config->fertilizer_input&&!config->fix_fertilization)
  {
    /* read fertilizer data */
    landuse->fertilizer_nr.fmt=config->fertilizer_nr_filename.fmt;
    if(config->fertilizer_nr_filename.fmt==CDF)
    {
      if(opendata_netcdf(&landuse->fertilizer_nr,config->fertilizer_nr_filename.name,config->fertilizer_nr_filename.var,NULL,config))
      {
        if(landuse->landuse.fmt==CDF)
          closeclimate_netcdf(&landuse->landuse,isroot(*config));
        else
          fclose(landuse->landuse.file);
        if(landuse->sdate.file!=NULL)
        {
          if(landuse->sdate.fmt==CDF)
            closeclimate_netcdf(&landuse->sdate,isroot(*config));
          else
            fclose(landuse->sdate.file);
        }
        if(landuse->crop_phu.file!=NULL)
        {
          if(landuse->crop_phu.fmt==CDF)
            closeclimate_netcdf(&landuse->crop_phu,isroot(*config));
          else
            fclose(landuse->crop_phu.file);
        }
        free(landuse);
        return NULL;
      }
      if(landuse->fertilizer_nr.var_len!=2*(ncft+NGRASS+NBIOMASSTYPE))
      {
        closeclimate_netcdf(&landuse->fertilizer_nr,isroot(*config));
        if(isroot(*config))
          fprintf(stderr,
                  "ERROR147: Invalid number of bands=%d in fertilizer Nr data file.\n",
                  (int)landuse->fertilizer_nr.var_len);
        if(landuse->landuse.fmt==CDF)
          closeclimate_netcdf(&landuse->landuse,isroot(*config));
        else
          fclose(landuse->landuse.file);
        if(landuse->sdate.file!=NULL)
        {
          if(landuse->sdate.fmt==CDF)
            closeclimate_netcdf(&landuse->sdate,isroot(*config));
          else
            fclose(landuse->sdate.file);
        }
        if(landuse->crop_phu.file!=NULL)
        {
          if(landuse->crop_phu.fmt==CDF)
            closeclimate_netcdf(&landuse->crop_phu,isroot(*config));
          else
            fclose(landuse->crop_phu.file);
        }
        free(landuse);
        return NULL;
      }
    }
    else
    {
      if((landuse->fertilizer_nr.file=openinputfile(&header,&landuse->fertilizer_nr.swap,
        &config->fertilizer_nr_filename,headername,
        &version,&offset,config))==NULL)
      {
        if(landuse->landuse.fmt==CDF)
          closeclimate_netcdf(&landuse->landuse,isroot(*config));
        else
          fclose(landuse->landuse.file);
        if(landuse->sdate.file!=NULL)
        {
          if(landuse->sdate.fmt==CDF)
            closeclimate_netcdf(&landuse->sdate,isroot(*config));
          else
            fclose(landuse->sdate.file);
        }
        if(landuse->crop_phu.file!=NULL)
        {
          if(landuse->crop_phu.fmt==CDF)
            closeclimate_netcdf(&landuse->crop_phu,isroot(*config));
          else
            fclose(landuse->crop_phu.file);
        }
        free(landuse);
        return NULL;
      }
      if(config->fertilizer_nr_filename.fmt==RAW)
      {
        header.nbands=2*(ncft+NGRASS+NBIOMASSTYPE);
        landuse->fertilizer_nr.offset=config->startgrid*header.nbands*sizeof(short);
      }
      else
      {
        if(header.nbands!=2*(ncft+NGRASS+NBIOMASSTYPE))
        {
          if(landuse->landuse.fmt==CDF)
            closeclimate_netcdf(&landuse->landuse,isroot(*config));
          else
            fclose(landuse->landuse.file);
          if(landuse->sdate.file!=NULL)
          {
            if(landuse->sdate.fmt==CDF)
              closeclimate_netcdf(&landuse->sdate,isroot(*config));
            else
              fclose(landuse->sdate.file);
          }
          if(landuse->crop_phu.file!=NULL)
          {
            if(landuse->crop_phu.fmt==CDF)
              closeclimate_netcdf(&landuse->crop_phu,isroot(*config));
            else
              fclose(landuse->crop_phu.file);
          }
          fclose(landuse->fertilizer_nr.file);
          if(isroot(*config))
            fprintf(stderr,
              "ERROR147: Invalid number of bands=%d in fertilizer Nr data file.\n",
              header.nbands);
          free(landuse);
          return(NULL);
        }
        landuse->fertilizer_nr.datatype=header.datatype;
        landuse->fertilizer_nr.offset=(config->startgrid-header.firstcell)*header.nbands*typesizes[header.datatype]+headersize(headername,version);
      }
      landuse->fertilizer_nr.firstyear=header.firstyear;
      landuse->fertilizer_nr.nyear=header.nyear;
      landuse->fertilizer_nr.size=header.ncell*header.nbands*typesizes[header.datatype];
      landuse->fertilizer_nr.n=config->ngridcell*header.nbands;
      landuse->fertilizer_nr.var_len=header.nbands;
      landuse->fertilizer_nr.scalar=header.scalar;
    }

    if(config->with_nitrogen&&config->manure_input&&!config->fix_fertilization)
    {
      /* read manure fertilizer data */
      landuse->manure_nr.fmt=config->manure_nr_filename.fmt;
      if(config->manure_nr_filename.fmt==CDF)
      {
        if(opendata_netcdf(&landuse->manure_nr,config->manure_nr_filename.name,config->manure_nr_filename.var,NULL,config))
        {
          if(landuse->landuse.fmt==CDF)
            closeclimate_netcdf(&landuse->landuse,isroot(*config));
          else
            fclose(landuse->landuse.file);
          if(landuse->sdate.file!=NULL)
          {
            if(landuse->sdate.fmt==CDF)
              closeclimate_netcdf(&landuse->sdate,isroot(*config));
            else
              fclose(landuse->sdate.file);
          }
          if(landuse->crop_phu.file!=NULL)
          {
            if(landuse->crop_phu.fmt==CDF)
              closeclimate_netcdf(&landuse->crop_phu,isroot(*config));
            else
              fclose(landuse->crop_phu.file);
          }
          if(landuse->fertilizer_nr.file!=NULL)
          {
            if(landuse->fertilizer_nr.fmt==CDF)
              closeclimate_netcdf(&landuse->fertilizer_nr,isroot(*config));
            else
              fclose(landuse->fertilizer_nr.file);
          }

          free(landuse);
          return NULL;
        }
        if(landuse->manure_nr.var_len!=2*(ncft+NGRASS+NBIOMASSTYPE))
        {
          closeclimate_netcdf(&landuse->manure_nr,isroot(*config));
          if(isroot(*config))
            fprintf(stderr,
              "ERROR147: Invalid number of bands=%d in fertilizer Nr data file.\n",
              (int)landuse->manure_nr.var_len);
          if(landuse->landuse.fmt==CDF)
            closeclimate_netcdf(&landuse->landuse,isroot(*config));
          else
            fclose(landuse->landuse.file);
          if(landuse->sdate.file!=NULL)
          {
            if(landuse->sdate.fmt==CDF)
              closeclimate_netcdf(&landuse->sdate,isroot(*config));
            else
              fclose(landuse->sdate.file);
            if(landuse->crop_phu.file!=NULL)
            {
              if(landuse->crop_phu.fmt==CDF)
                closeclimate_netcdf(&landuse->crop_phu,isroot(*config));
              else
                fclose(landuse->crop_phu.file);
            }
            if(landuse->fertilizer_nr.file!=NULL)
            {
              if(landuse->fertilizer_nr.fmt==CDF)
                closeclimate_netcdf(&landuse->fertilizer_nr,isroot(*config));
              else
                fclose(landuse->fertilizer_nr.file);
            }
          }
          free(landuse);
          return NULL;
        }
      }
      else
      {
        if((landuse->manure_nr.file=openinputfile(&header,&landuse->manure_nr.swap,
          &config->manure_nr_filename,headername,
          &version,&offset,config))==NULL)
        {
          if(landuse->landuse.fmt==CDF)
            closeclimate_netcdf(&landuse->landuse,isroot(*config));
          else
            fclose(landuse->landuse.file);
          if(landuse->sdate.file!=NULL)
          {
            if(landuse->sdate.fmt==CDF)
              closeclimate_netcdf(&landuse->sdate,isroot(*config));
            else
              fclose(landuse->sdate.file);
          }
          if(landuse->crop_phu.file!=NULL)
          {
            if(landuse->crop_phu.fmt==CDF)
              closeclimate_netcdf(&landuse->crop_phu,isroot(*config));
            else
              fclose(landuse->crop_phu.file);
          }
          if(landuse->fertilizer_nr.file!=NULL)
          {
            if(landuse->fertilizer_nr.fmt==CDF)
              closeclimate_netcdf(&landuse->fertilizer_nr,isroot(*config));
            else
              fclose(landuse->fertilizer_nr.file);
          }
          free(landuse);
          return NULL;
        }
        if(config->manure_nr_filename.fmt==RAW)
        {
          header.nbands=2*(ncft+NGRASS+NBIOMASSTYPE);
          landuse->manure_nr.offset=config->startgrid*header.nbands*sizeof(short);
        }
        else
        {
          if(header.nbands!=2*(ncft+NGRASS+NBIOMASSTYPE))
          {
            if(landuse->landuse.fmt==CDF)
              closeclimate_netcdf(&landuse->landuse,isroot(*config));
            else
              fclose(landuse->landuse.file);
            if(landuse->sdate.file!=NULL)
            {
              if(landuse->sdate.fmt==CDF)
                closeclimate_netcdf(&landuse->sdate,isroot(*config));
              else
                fclose(landuse->sdate.file);
            }
            if(landuse->crop_phu.file!=NULL)
            {
              if(landuse->crop_phu.fmt==CDF)
                closeclimate_netcdf(&landuse->crop_phu,isroot(*config));
              else
                fclose(landuse->crop_phu.file);
            }
            if(landuse->fertilizer_nr.file!=NULL)
            {
              if(landuse->fertilizer_nr.fmt==CDF)
                closeclimate_netcdf(&landuse->fertilizer_nr,isroot(*config));
              else
                fclose(landuse->fertilizer_nr.file);
            }
            fclose(landuse->manure_nr.file);
            if(isroot(*config))
              fprintf(stderr,
                "ERROR147: Invalid number of bands=%d in manure data file.\n",
                header.nbands);
            free(landuse);
            return(NULL);
          }
          landuse->manure_nr.datatype=header.datatype;
          landuse->manure_nr.offset=(config->startgrid-header.firstcell)*header.nbands*typesizes[header.datatype]+headersize(headername,version);
        }
        landuse->manure_nr.firstyear=header.firstyear;
        landuse->manure_nr.nyear=header.nyear;
        landuse->manure_nr.size=header.ncell*header.nbands*typesizes[header.datatype];
        landuse->manure_nr.n=config->ngridcell*header.nbands;
        landuse->manure_nr.var_len=header.nbands;
        landuse->manure_nr.scalar=header.scalar;
      }
    }
  }
  else
  {
    landuse->fertilizer_nr.file=NULL;
    landuse->manure_nr.file=NULL;
  }

  if(config->tillage_type==READ_TILLAGE)
  {
    landuse->with_tillage.fmt=config->with_tillage_filename.fmt;
    if(config->with_tillage_filename.fmt==CDF)
    {
      if(opendata_netcdf(&landuse->with_tillage,config->with_tillage_filename.name,config->with_tillage_filename.var,NULL,config))
      {
        if(landuse->landuse.fmt==CDF)
          closeclimate_netcdf(&landuse->landuse,isroot(*config));
        else
          fclose(landuse->landuse.file);
        if(landuse->sdate.file!=NULL)
        {
          if(landuse->sdate.fmt==CDF)
            closeclimate_netcdf(&landuse->sdate,isroot(*config));
          else
            fclose(landuse->sdate.file);
        }
        if(landuse->crop_phu.file!=NULL)
        {
          if(landuse->crop_phu.fmt==CDF)
            closeclimate_netcdf(&landuse->crop_phu,isroot(*config));
          else
            fclose(landuse->crop_phu.file);
        }
        if(landuse->fertilizer_nr.file!=NULL)
        {
          if(landuse->fertilizer_nr.fmt==CDF)
            closeclimate_netcdf(&landuse->fertilizer_nr,isroot(*config));
          else
            fclose(landuse->fertilizer_nr.file);
        }
        if(landuse->manure_nr.file!=NULL)
        {
          if(landuse->manure_nr.fmt==CDF)
            closeclimate_netcdf(&landuse->manure_nr,isroot(*config));
          else
            fclose(landuse->manure_nr.file);
        }
        free(landuse);
        return NULL;
      }
    }
    else
    {
      if((landuse->with_tillage.file=openinputfile(&header,&landuse->with_tillage.swap,
        &config->with_tillage_filename,headername,
        &version,&offset,config))==NULL)
      {
        if(landuse->landuse.fmt==CDF)
          closeclimate_netcdf(&landuse->landuse,isroot(*config));
        else
          fclose(landuse->landuse.file);
        fclose(landuse->with_tillage.file);
        if(landuse->sdate.file!=NULL)
        {
          if(landuse->sdate.fmt==CDF)
            closeclimate_netcdf(&landuse->sdate,isroot(*config));
          else
            fclose(landuse->sdate.file);
        }
        if(landuse->crop_phu.file!=NULL)
        {
          if(landuse->crop_phu.fmt==CDF)
            closeclimate_netcdf(&landuse->crop_phu,isroot(*config));
          else
            fclose(landuse->crop_phu.file);
        }
        if(landuse->fertilizer_nr.file!=NULL)
        {
          if(landuse->fertilizer_nr.fmt==CDF)
            closeclimate_netcdf(&landuse->fertilizer_nr,isroot(*config));
          else
            fclose(landuse->fertilizer_nr.file);
        }
        if(landuse->manure_nr.file!=NULL)
        {
          if(landuse->manure_nr.fmt==CDF)
            closeclimate_netcdf(&landuse->manure_nr,isroot(*config));
          else
            fclose(landuse->manure_nr.file);
        }
        free(landuse);
        return NULL;
      }
      if(config->with_tillage_filename.fmt==RAW)
      {
        landuse->with_tillage.offset=config->startgrid*header.nbands*sizeof(short);
      }
      else
      {
        if(header.nbands!=1)
        {
          if(landuse->landuse.fmt==CDF)
            closeclimate_netcdf(&landuse->landuse,isroot(*config));
          else
            fclose(landuse->landuse.file);
          fclose(landuse->with_tillage.file);
          if(landuse->sdate.file!=NULL)
          {
            if(landuse->sdate.fmt==CDF)
              closeclimate_netcdf(&landuse->sdate,isroot(*config));
            else
              fclose(landuse->sdate.file);
          }
          if(landuse->crop_phu.file!=NULL)
          {
            if(landuse->crop_phu.fmt==CDF)
              closeclimate_netcdf(&landuse->crop_phu,isroot(*config));
            else
              fclose(landuse->crop_phu.file);
          }
          if(landuse->fertilizer_nr.file!=NULL)
          {
            if(landuse->fertilizer_nr.fmt==CDF)
              closeclimate_netcdf(&landuse->fertilizer_nr,isroot(*config));
            else
              fclose(landuse->fertilizer_nr.file);
          }
          if(landuse->manure_nr.file!=NULL)
          {
            if(landuse->manure_nr.fmt==CDF)
              closeclimate_netcdf(&landuse->manure_nr,isroot(*config));
            else
              fclose(landuse->manure_nr.file);
          }
          if(isroot(*config))
            fprintf(stderr,
              "ERROR147: Invalid number of bands=%d in tillage type file.\n",
              header.nbands);
          free(landuse);
          return(NULL);
        }
        landuse->with_tillage.offset=(config->startgrid-header.firstcell)*header.nbands*sizeof(short)+headersize(headername,version);
      }
      landuse->with_tillage.firstyear=header.firstyear;
      landuse->with_tillage.nyear=header.nyear;
      landuse->with_tillage.size=header.ncell*header.nbands*sizeof(short);
      landuse->with_tillage.n=config->ngridcell*header.nbands;
      landuse->with_tillage.var_len=header.nbands;
      landuse->with_tillage.scalar=header.scalar;
    }
  }
  else
    landuse->with_tillage.file=NULL;

  if(config->residue_treatment==READ_RESIDUE_DATA)
  {
    /* read residue data */
    landuse->residue_on_field.fmt=config->residue_data_filename.fmt;
    if(config->residue_data_filename.fmt==CDF)
    {
      if(opendata_netcdf(&landuse->residue_on_field,config->residue_data_filename.name,config->residue_data_filename.var,NULL,config))
      {
        if(landuse->landuse.fmt==CDF)
          closeclimate_netcdf(&landuse->landuse,isroot(*config));
        else
          fclose(landuse->landuse.file);
        if(landuse->sdate.file!=NULL)
        {
          if(landuse->sdate.fmt==CDF)
            closeclimate_netcdf(&landuse->sdate,isroot(*config));
          else
            fclose(landuse->sdate.file);
        }
        if(landuse->crop_phu.file!=NULL)
        {
          if(landuse->crop_phu.fmt==CDF)
            closeclimate_netcdf(&landuse->crop_phu,isroot(*config));
          else
            fclose(landuse->crop_phu.file);
        }
        if(landuse->fertilizer_nr.file!=NULL)
        {
          if(landuse->fertilizer_nr.fmt==CDF)
            closeclimate_netcdf(&landuse->fertilizer_nr,isroot(*config));
          else
            fclose(landuse->fertilizer_nr.file);
        }
        if(landuse->manure_nr.file!=NULL)
        {
          if(landuse->manure_nr.fmt==CDF)
            closeclimate_netcdf(&landuse->manure_nr,isroot(*config));
          else
            fclose(landuse->manure_nr.file);
        }
        if(landuse->with_tillage.file!=NULL)
        {
          if(landuse->with_tillage.fmt==CDF)
            closeclimate_netcdf(&landuse->with_tillage,isroot(*config));
          else
            fclose(landuse->with_tillage.file);
        }
        free(landuse);
        return NULL;
      }
      if(landuse->residue_on_field.var_len!=ncft+NGRASS+NBIOMASSTYPE)
      {
        closeclimate_netcdf(&landuse->residue_on_field,isroot(*config));
        if(isroot(*config))
          fprintf(stderr,
            "ERROR147: Invalid number of bands=%d in residue extraction data file.\n",
            (int)landuse->residue_on_field.var_len);
        if(landuse->landuse.fmt==CDF)
          closeclimate_netcdf(&landuse->landuse,isroot(*config));
        else
          fclose(landuse->landuse.file);
        if(landuse->sdate.file!=NULL)
        {
          if(landuse->sdate.fmt==CDF)
            closeclimate_netcdf(&landuse->sdate,isroot(*config));
          else
            fclose(landuse->sdate.file);
        }
        if(landuse->crop_phu.file!=NULL)
        {
          if(landuse->crop_phu.fmt==CDF)
            closeclimate_netcdf(&landuse->crop_phu,isroot(*config));
          else
            fclose(landuse->crop_phu.file);
        }
        if(landuse->fertilizer_nr.file!=NULL)
        {
          if(landuse->fertilizer_nr.fmt==CDF)
            closeclimate_netcdf(&landuse->fertilizer_nr,isroot(*config));
          else
            fclose(landuse->fertilizer_nr.file);
        }
        if(landuse->manure_nr.file!=NULL)
        {
          if(landuse->manure_nr.fmt==CDF)
            closeclimate_netcdf(&landuse->manure_nr,isroot(*config));
          else
            fclose(landuse->manure_nr.file);
        }
        if(landuse->with_tillage.file!=NULL)
        {
          if(landuse->with_tillage.fmt==CDF)
            closeclimate_netcdf(&landuse->with_tillage,isroot(*config));
          else
            fclose(landuse->with_tillage.file);
        }
        free(landuse);
        return NULL;
      }
    }
    else
    {
      if((landuse->residue_on_field.file=openinputfile(&header,&landuse->residue_on_field.swap,
        &config->residue_data_filename,headername,
        &version,&offset,config))==NULL)
      {
        if(landuse->landuse.fmt==CDF)
          closeclimate_netcdf(&landuse->landuse,isroot(*config));
        else
          fclose(landuse->landuse.file);
        if(landuse->sdate.file!=NULL)
        {
          if(landuse->sdate.fmt==CDF)
            closeclimate_netcdf(&landuse->sdate,isroot(*config));
          else
            fclose(landuse->sdate.file);
        }
        if(landuse->crop_phu.file!=NULL)
        {
          if(landuse->crop_phu.fmt==CDF)
            closeclimate_netcdf(&landuse->crop_phu,isroot(*config));
          else
            fclose(landuse->crop_phu.file);
        }
        if(landuse->fertilizer_nr.file!=NULL)
        {
          if(landuse->fertilizer_nr.fmt==CDF)
            closeclimate_netcdf(&landuse->fertilizer_nr,isroot(*config));
          else
            fclose(landuse->fertilizer_nr.file);
        }
        if(landuse->manure_nr.file!=NULL)
        {
          if(landuse->manure_nr.fmt==CDF)
            closeclimate_netcdf(&landuse->manure_nr,isroot(*config));
          else
            fclose(landuse->manure_nr.file);
        }
        if(landuse->with_tillage.file!=NULL)
        {
          if(landuse->with_tillage.fmt==CDF)
            closeclimate_netcdf(&landuse->with_tillage,isroot(*config));
          else
            fclose(landuse->with_tillage.file);
        }
        free(landuse);
        return NULL;
      }
      if(config->residue_data_filename.fmt==RAW)
      {
        header.nbands=ncft+NGRASS+NBIOMASSTYPE;
        landuse->residue_on_field.offset=config->startgrid*header.nbands*sizeof(short);
      }
      else
      {
        if(header.nbands!=ncft+NGRASS+NBIOMASSTYPE)
        {
          if(landuse->landuse.fmt==CDF)
            closeclimate_netcdf(&landuse->landuse,isroot(*config));
          else
            fclose(landuse->landuse.file);
          if(landuse->sdate.file!=NULL)
          {
            if(landuse->sdate.fmt==CDF)
              closeclimate_netcdf(&landuse->sdate,isroot(*config));
            else
              fclose(landuse->sdate.file);
          }
          if(landuse->crop_phu.file!=NULL)
          {
            if(landuse->crop_phu.fmt==CDF)
              closeclimate_netcdf(&landuse->crop_phu,isroot(*config));
            else
              fclose(landuse->crop_phu.file);
          }
          if(landuse->fertilizer_nr.file!=NULL)
          {
            if(landuse->fertilizer_nr.fmt==CDF)
              closeclimate_netcdf(&landuse->fertilizer_nr,isroot(*config));
            else
              fclose(landuse->fertilizer_nr.file);
          }
          if(landuse->manure_nr.file!=NULL)
          {
            if(landuse->manure_nr.fmt==CDF)
              closeclimate_netcdf(&landuse->manure_nr,isroot(*config));
            else
              fclose(landuse->manure_nr.file);
          }
          if(landuse->with_tillage.file!=NULL)
          {
            if(landuse->with_tillage.fmt==CDF)
              closeclimate_netcdf(&landuse->with_tillage,isroot(*config));
            else
              fclose(landuse->with_tillage.file);
          }
          fclose(landuse->residue_on_field.file);
          if(isroot(*config))
            fprintf(stderr,
              "ERROR147: Invalid number of bands=%d in residue extract data file.\n",
              header.nbands);
          free(landuse);
          return(NULL);
        }
        landuse->residue_on_field.datatype=header.datatype;
        landuse->residue_on_field.offset=(config->startgrid-header.firstcell)*header.nbands*typesizes[header.datatype]+headersize(headername,version);
      }
      landuse->residue_on_field.firstyear=header.firstyear;
      landuse->residue_on_field.nyear=header.nyear;
      landuse->residue_on_field.size=header.ncell*header.nbands*typesizes[header.datatype];
      landuse->residue_on_field.n=config->ngridcell*header.nbands;
      landuse->residue_on_field.var_len=header.nbands;
      landuse->residue_on_field.scalar=header.scalar;
    }
  }
  else
    landuse->residue_on_field.file=NULL;

  landuse->intercrop=config->intercrop;
  return landuse;
} /* of 'initlanduse' */


/**************************************************************************************/
/**                                                                                \n**/
/**     C implementation of LPJmL                                                  \n**/
/**                                                                                \n**/
/**     Function reads landuse data for a specific year                            \n**/
/**                                                                                \n**/
/**     - order of landuse input data:                                             \n**/
/**        0-10   CFT                                                              \n**/
/**          11   OTHERS                                                           \n**/
/**          12   PASTURES                                                         \n**/
/**          13   BIOMASS_GRASS                                                    \n**/
/**          14   BIOMASS_TREE                                                     \n**/
/**       15-25   CFT_irr                                                          \n**/
/**          26   others_irr                                                       \n**/
/**          27   PASTURE_irr                                                      \n**/
/**          28   BIOMASS_GRASS IRR                                                \n**/
/**          29   BIOMASS_TREE IRR                                                 \n**/
/**     - called in iterate()                                                      \n**/
/**     - reads every year the fractions of the bands for all cells from           \n**/
/**       the input file                                                           \n**/
/**     - checks if sum of fraction is not greater 1                               \n**/
/**       -> if sum of fraction is greater 1: subtraction from fraction            \n**/
/**          of managed grass if possible                                          \n**/
/**       -> else fail incorrect input file                                        \n**/
/**                                                                                \n**/
/**************************************************************************************/

static Real reducelanduse(Cell *cell,Real sum,int ncft)
{
  int i,j;
  if(cell->ml.landfrac[0].grass[1]>sum)
  {
    cell->ml.landfrac[0].grass[1]-=sum;
    return 0.0;
  }
  if(cell->ml.landfrac[1].grass[1]>sum)
  {
    cell->ml.landfrac[1].grass[1]-=sum;
    return 0.0;
  }
  for(j=0; j<2; j++)
  {
    for(i=0; i<ncft; i++)
      if(cell->ml.landfrac[j].crop[i]>sum)
      {
        cell->ml.landfrac[j].crop[i]-=sum;
        return 0;
      }
    for(i=0; i<NGRASS; i++)
      if(cell->ml.landfrac[j].grass[i]>sum)
      {
        cell->ml.landfrac[j].grass[i]-=sum;
        return 0;
      }
    if(cell->ml.landfrac[j].biomass_tree>sum)
    {
      cell->ml.landfrac[j].biomass_tree-=sum;
      return 0;
    }
    if(cell->ml.landfrac[j].biomass_grass>sum)
    {
      cell->ml.landfrac[j].biomass_grass-=sum;
      return 0;
    }
  }
  return sum;
} /* of 'reducelanduse' */

Bool getlanduse(Landuse landuse,     /**< Pointer to landuse data */
  Cell grid[],         /**< LPJ cell array */
  int year,            /**< year (AD) */
  int actual_year,     /**< year (AD) but not the static in case of CONST_LANDUSE */
  int ncft,            /**< number of crop PFTs */
  const Config *config /**< LPJ configuration */
)                     /** \return TRUE on error */
{
  short *vec;
  int i,j,p,count,cell;
  Real sum,*data,*fert_nr,*manu_nr,*res_on_field;
  int *dates;
  int *phus;
  Bool *tilltypes;
  int yearsdate=actual_year;     /*sdate year*/
  int yearphu=actual_year;       /*crop phu year*/
  int yearf=year;
  int yearm=year;
  int yeart=year;
  int yearr=year;
  /* for testing soil type to avoid all crops on ROCK and ICE cells */
  Stand *stand;
  int soiltype=-1;
  int yearl=year;

  /* Initialize yearly prescribed sdate */
  if(config->sdate_option==PRESCRIBED_SDATE)
  {
    /* assigning sdate data */
    yearsdate-=landuse->sdate.firstyear;
    if(yearsdate>=landuse->sdate.nyear)
      yearsdate=landuse->sdate.nyear-1; /* use last year sdate */
    else if(yearsdate<0)
      yearsdate=0;                        /* use first year sdate */

    if(landuse->sdate.fmt==CDF)
    {
      dates=newvec(int,config->ngridcell*landuse->sdate.var_len);
      if(dates==NULL)
      {
        printallocerr("dates");
        return TRUE;
      }
      if(readintdata_netcdf(&landuse->sdate,dates,grid,yearsdate,config))
      {
        fprintf(stderr,
          "ERROR149: Cannot read sowing dates of year %d in getlanduse().\n",
          yearsdate+landuse->sdate.firstyear);
        fflush(stderr);
        return TRUE;
      }
      count=0;
      for(cell=0; cell<config->ngridcell; cell++)
        if(!grid[cell].skip)
          for(j=0; j<2*ncft; j++)
            grid[cell].ml.sdate_fixed[j]=dates[count++];
        else
          count+=2*ncft;
      free(dates);
    }
    else
    {
      if(fseek(landuse->sdate.file,(long long)yearsdate*landuse->sdate.size+landuse->sdate.offset,SEEK_SET))
      {
        fprintf(stderr,
          "ERROR148: Cannot seek sowing dates to year %d in getlanduse().\n",
          yearsdate);
        return TRUE;
      }
      dates=newvec(int,landuse->sdate.n);
      if(dates==NULL)
      {
        printallocerr("dates");
        return TRUE;
      }
      if(readintvec(landuse->sdate.file,dates,landuse->sdate.n,landuse->sdate.swap,landuse->sdate.datatype))
      {
        fprintf(stderr,
          "ERROR149: Cannot read sowing dates of year %d in getlanduse().\n",
          yearsdate);
        free(dates);
        return TRUE;
      }
      count=0;
      for(cell=0; cell<config->ngridcell; cell++)
        if(!grid[cell].skip)
        {
          for(j=0; j<2*ncft; j++)
          {
            grid[cell].ml.sdate_fixed[j]=dates[count++];
          }
        }
        else
          count+=2*ncft;
      free(dates);
    }
  }

  /* Initialize yearly prescribed crop_phu */
  if(config->crop_phu_option==PRESCRIBED_CROP_PHU)
  {
    /* assigning crop phus data */
    yearphu-=landuse->crop_phu.firstyear;
    if(yearphu>=landuse->crop_phu.nyear)
      yearphu=landuse->crop_phu.nyear-1; /* use last year sdate */
    else if(yearphu<0)
      yearphu=0;                        /* use first year sdate */

    if(landuse->crop_phu.fmt==CDF)
    {
      phus=newvec(int,config->ngridcell*landuse->crop_phu.var_len);
      if(phus==NULL)
      {
        printallocerr("dates");
        return TRUE;
      }
      if(readintdata_netcdf(&landuse->crop_phu,phus,grid,yearphu,config))
      {
        fprintf(stderr,
          "ERROR149: Cannot read crop phus of year %d in getlanduse().\n",
          yearphu+landuse->crop_phu.firstyear);
        fflush(stderr);
        return TRUE;
      }
      count=0;
      for(cell=0; cell<config->ngridcell; cell++)
        if(!grid[cell].skip)
          for(j=0; j<2*ncft; j++)
            grid[cell].ml.crop_phu_fixed[j]=(Real)phus[count++];
        else
          count+=2*ncft;
      free(phus);
    }
    else
    {
      if(fseek(landuse->crop_phu.file,(long long)yearphu*landuse->crop_phu.size+landuse->crop_phu.offset,SEEK_SET))
      {
        fprintf(stderr,
          "ERROR148: Cannot seek crop phus to year %d in getlanduse().\n",
          yearphu);
        return TRUE;
      }
      phus=newvec(int,landuse->crop_phu.n);
      if(phus==NULL)
      {
        printallocerr("dates");
        return TRUE;
      }
      if(readintvec(landuse->crop_phu.file,phus,landuse->crop_phu.n,landuse->crop_phu.swap,landuse->crop_phu.datatype))
      {
        fprintf(stderr,
          "ERROR149: Cannot read crop phus of year %d in getlanduse().\n",
          yearphu);
        free(phus);
        return TRUE;
      }
      count=0;
      for(cell=0; cell<config->ngridcell; cell++)
        if(!grid[cell].skip)
        {
          for(j=0; j<2*ncft; j++)
          {
            grid[cell].ml.crop_phu_fixed[j]=(Real)phus[count++];
          }
        }
        else
          count+=2*ncft;
      free(phus);
    }
  } /* end crop_phu*/

  yearl-=landuse->landuse.firstyear;
  if(yearl>=landuse->landuse.nyear)
    yearl=landuse->landuse.nyear-1;
  else if(yearl<0)
    yearl=0;
  if(landuse->landuse.fmt==CDF)
  {
    data=newvec(Real,config->ngridcell*landuse->landuse.var_len);
    if(data==NULL)
    {
      printallocerr("data");
      return TRUE;
    }
    if(readdata_netcdf(&landuse->landuse,data,grid,yearl,config))
    {
      fprintf(stderr,
        "ERROR149: Cannot read landuse of year %d in getlanduse().\n",
        yearl+landuse->landuse.firstyear);
      fflush(stderr);
      return TRUE;
    }
  }
  else
  {
    if(fseek(landuse->landuse.file,(long long)yearl*landuse->landuse.size+landuse->landuse.offset,SEEK_SET))
    {
      fprintf(stderr,
        "ERROR148: Cannot seek landuse to year %d in getlanduse().\n",
        yearl+landuse->landuse.firstyear);
      fflush(stderr);
      return TRUE;
    }
    data=newvec(Real,landuse->landuse.n);
    if(data==NULL)
    {
      printallocerr("data");
      return TRUE;
    }
    if(readrealvec(landuse->landuse.file,data,0,landuse->landuse.scalar,landuse->landuse.n,landuse->landuse.swap,landuse->landuse.datatype))
    {
      fprintf(stderr,
        "ERROR149: Cannot read landuse of year %d in getlanduse().\n",
        yearl+landuse->landuse.firstyear);
      fflush(stderr);
      free(data);
      return TRUE;
    }
  }
  count=0;


  for(cell=0; cell<config->ngridcell; cell++)
  {
  //  if (!grid[cell].skip){
    //  fprintf(stdout,"in landuse 1317 per cell %i\n",cell);
    //  fflush(stdout);
      /* get soiltype of first stand (not stored in cell structure) */
      if(!grid[cell].skip && grid[cell].standlist->n>0)
      {
        stand=getstand(grid[cell].standlist,0);
        soiltype=stand->soil.par->type;
      }
      else
      {
        soiltype=-1;
      }
    //  fprintf(stdout,"in landuse 1328 per cell soil type %i\n",soiltype);
    //  fflush(stdout);
      for(i=0; i<WIRRIG; i++)
      {
        /* read cropfrac from 32 bands or rain-fed cropfrac from 64 bands input */
        if(landuse->nbands!=4*(ncft+NGRASS+NBIOMASSTYPE)||i<1)
        {
          for(j=0; j<ncft; j++)
          {
            grid[cell].ml.landfrac[i].crop[j]=data[count++];
            if(i>0&&!grid[cell].skip)
              grid[cell].ml.irrig_system->crop[j]=grid[cell].ml.manage.par->default_irrig_system; /*default national irrigation system (Rohwer & Gerten 2007)*/
          }
          for(j=0; j<NGRASS; j++)
          {
            grid[cell].ml.landfrac[i].grass[j]=data[count++];
            if(i>0&&!grid[cell].skip)
              grid[cell].ml.irrig_system->grass[j]=grid[cell].ml.manage.par->default_irrig_system;
          }
          if(landuse->nbands!=2*(ncft+NGRASS))
          {
            grid[cell].ml.landfrac[i].biomass_grass=data[count++];
            if(i>0&&!grid[cell].skip)
              grid[cell].ml.irrig_system->biomass_grass=grid[cell].ml.manage.par->default_irrig_system;
            grid[cell].ml.landfrac[i].biomass_tree=data[count++];
            if(i>0&&!grid[cell].skip)
              grid[cell].ml.irrig_system->biomass_tree=grid[cell].ml.manage.par->default_irrig_system;
          }
          else
            grid[cell].ml.landfrac[i].biomass_grass=grid[cell].ml.landfrac[i].biomass_tree=0;
        }
        else /* read irrigated cropfrac and irrigation systems from 64 bands input */
        {
          for(j=0; j<ncft; j++) /* initialization needed */
            grid[cell].ml.landfrac[i].crop[j]=0;
          for(j=0; j<NGRASS; j++)
            grid[cell].ml.landfrac[i].grass[j]=0;
          grid[cell].ml.landfrac[i].biomass_grass=0;
          grid[cell].ml.landfrac[i].biomass_tree=0;

          for(p=1; p<4; p++) /* irrigation system loop; 1: surface, 2: sprinkler, 3: drip */
          {
            for(j=0; j<ncft; j++)
            {
              if(data[count]>0)
              {
                grid[cell].ml.landfrac[i].crop[j]=data[count++];
                grid[cell].ml.irrig_system->crop[j]=p;
              }
              else
                count++;
            }
            for(j=0; j<NGRASS; j++)
            {
              if(data[count]>0)
              {
                grid[cell].ml.landfrac[i].grass[j]=data[count++];
                grid[cell].ml.irrig_system->grass[j]=p;
              }
              else
                count++;
            }
            if(data[count]>0)
            {
              grid[cell].ml.landfrac[i].biomass_grass=data[count++];
              grid[cell].ml.irrig_system->biomass_grass=p;
            }
            else
              count++;
            if(data[count]>0)
            {
              grid[cell].ml.landfrac[i].biomass_tree=data[count++];
              grid[cell].ml.irrig_system->biomass_tree=p;
            }
            else
              count++;
          }
        }
      }
      switch(config->irrig_scenario)
      {
      case NO_IRRIGATION:
        for(j=0; j<ncft; j++)
        {
          grid[cell].ml.landfrac[0].crop[j]+=grid[cell].ml.landfrac[1].crop[j];
          grid[cell].ml.landfrac[1].crop[j]=0;
          grid[cell].ml.irrig_system->crop[j]=NOIRRIG;
        }
        for(j=0; j<NGRASS; j++)
        {
          grid[cell].ml.landfrac[0].grass[j]+=grid[cell].ml.landfrac[1].grass[j];
          grid[cell].ml.landfrac[1].grass[j]=0;
          grid[cell].ml.irrig_system->grass[j]=NOIRRIG;
        }
        grid[cell].ml.landfrac[0].biomass_grass+=grid[cell].ml.landfrac[1].biomass_grass;
        grid[cell].ml.landfrac[1].biomass_grass=0;
        grid[cell].ml.irrig_system->biomass_grass=NOIRRIG;
        grid[cell].ml.landfrac[0].biomass_tree+=grid[cell].ml.landfrac[1].biomass_tree;
        grid[cell].ml.landfrac[1].biomass_tree=0;
        grid[cell].ml.irrig_system->biomass_tree=NOIRRIG;
        break;
      case ALL_IRRIGATION:
        for(j=0; j<ncft; j++)
        {
          grid[cell].ml.landfrac[1].crop[j]+=grid[cell].ml.landfrac[0].crop[j];
          grid[cell].ml.landfrac[0].crop[j]=0;
          if(!grid[cell].skip)
            grid[cell].ml.irrig_system->crop[j]=grid[cell].ml.manage.par->default_irrig_system; /*default national irrigation system (Rohwer & Gerten 2007)*/
        }
        for(j=0; j<NGRASS; j++)
        {
          grid[cell].ml.landfrac[1].grass[j]+=grid[cell].ml.landfrac[0].grass[j];
          grid[cell].ml.landfrac[0].grass[j]=0;
          if(!grid[cell].skip)
            grid[cell].ml.irrig_system->grass[j]=grid[cell].ml.manage.par->default_irrig_system;
        }
        grid[cell].ml.landfrac[1].biomass_grass+=grid[cell].ml.landfrac[0].biomass_grass;
        grid[cell].ml.landfrac[0].biomass_grass=0;
        if(!grid[cell].skip)
          grid[cell].ml.irrig_system->biomass_grass=grid[cell].ml.manage.par->default_irrig_system;
        grid[cell].ml.landfrac[1].biomass_tree+=grid[cell].ml.landfrac[0].biomass_tree;
        grid[cell].ml.landfrac[0].biomass_tree=0;
        grid[cell].ml.irrig_system->biomass_tree=grid[cell].ml.manage.par->default_irrig_system;
        break;
      }   /* of switch(...) */

      /* DEBUG: here you can set land-use fractions manually, it overwrites the land-use input, in all cells */
      /* the irrigation system is set to the default country value, but you can set a number from 1-3 manually */
      /* 1: surface, 2: sprinkler, 3: drip irrigation */
  
/*
      sum=landfrac_sum(grid[cell].ml.landfrac,ncft,FALSE)+landfrac_sum(grid[cell].ml.landfrac,ncft,TRUE);

      for(j=0;j<ncft;j++)
      {
        grid[cell].ml.landfrac[1].crop[j]=0.0;
        grid[cell].ml.landfrac[0].crop[j]=0.0;
        grid[cell].ml.irrig_system->crop[j]=grid[cell].ml.manage.par->default_irrig_system;
      }
      for(j=0;j<NGRASS;j++)
      {
        grid[cell].ml.landfrac[1].grass[j]=0.0;
        grid[cell].ml.landfrac[0].grass[j]=0.0;
        grid[cell].ml.irrig_system->grass[j]=grid[cell].ml.manage.par->default_irrig_system;
      }

      grid[cell].ml.landfrac[1].biomass_grass=0.0;
      grid[cell].ml.landfrac[0].biomass_grass=0.0;
      grid[cell].ml.irrig_system->biomass_grass=grid[cell].ml.manage.par->default_irrig_system;
      grid[cell].ml.landfrac[1].biomass_tree=0.0;
      grid[cell].ml.landfrac[0].biomass_tree=0.0;
      grid[cell].ml.irrig_system->biomass_tree=grid[cell].ml.manage.par->default_irrig_system;

      grid[cell].ml.landfrac[1].grass[1]=0.0;
      grid[cell].ml.landfrac[0].grass[0]=1.0;
      //if (sum>1.00001) grid[cell].ml.landfrac[0].grass[0]=1.0;
*/
/* END DEBUG */

    if(config->others_to_crop)
    {
      if(grid[cell].coord.lat>30||grid[cell].coord.lat<-30)
      {
        grid[cell].ml.landfrac[0].crop[0]+=grid[cell].ml.landfrac[0].grass[0];
        grid[cell].ml.landfrac[1].crop[0]+=grid[cell].ml.landfrac[1].grass[0];
        grid[cell].ml.landfrac[0].grass[0]=grid[cell].ml.landfrac[1].grass[0]=0;
      }
      else
      {
        grid[cell].ml.landfrac[0].crop[2]+=grid[cell].ml.landfrac[0].grass[0];
        grid[cell].ml.landfrac[1].crop[2]+=grid[cell].ml.landfrac[1].grass[0];
        grid[cell].ml.landfrac[0].grass[0]=grid[cell].ml.landfrac[1].grass[0]=0;
      }
    }
    if(config->grassonly)
    {
      for(j=0; j<ncft; j++)
        grid[cell].ml.landfrac[0].crop[j]=grid[cell].ml.landfrac[1].crop[j]=0;
      grid[cell].ml.landfrac[0].grass[0]=grid[cell].ml.landfrac[1].grass[0]=0;
      grid[cell].ml.landfrac[0].biomass_grass=grid[cell].ml.landfrac[1].biomass_grass=
        grid[cell].ml.landfrac[0].biomass_tree=grid[cell].ml.landfrac[1].biomass_tree=0;
    }
    
    /* force tinyfrac for all crops only on pixels with valid soil */
    if (landuse->allcrops && !grid[cell].skip && soiltype!=ROCK && soiltype!=ICE && soiltype >= 0)
    {
      for(j=0; j<ncft; j++)
      {
        if (grid[cell].ml.landfrac[1].crop[j] < tinyfrac) 
        {
          grid[cell].ml.irrig_system->crop[j] = grid[cell].ml.manage.par->default_irrig_system;
          grid[cell].ml.landfrac[1].crop[j] = tinyfrac;
        }
        if (grid[cell].ml.landfrac[0].crop[j] < tinyfrac) grid[cell].ml.landfrac[0].crop[j] = tinyfrac;
      }
      for(j=0; j<NGRASS; j++)
      {
        if (grid[cell].ml.landfrac[0].grass[j] < tinyfrac) grid[cell].ml.landfrac[0].grass[j] = tinyfrac;
        if (grid[cell].ml.landfrac[1].grass[j] < tinyfrac) 
        {
          grid[cell].ml.landfrac[1].grass[j] = tinyfrac;
          grid[cell].ml.irrig_system->grass[j] = grid[cell].ml.manage.par->default_irrig_system;
        }
      }
      if (grid[cell].ml.landfrac[1].biomass_tree < tinyfrac) 
      {
        grid[cell].ml.landfrac[1].biomass_tree = tinyfrac;
        grid[cell].ml.irrig_system->biomass_tree = grid[cell].ml.manage.par->default_irrig_system;
      }
      if (grid[cell].ml.landfrac[0].biomass_tree < tinyfrac) grid[cell].ml.landfrac[0].biomass_tree = tinyfrac;
      if (grid[cell].ml.landfrac[1].biomass_grass < tinyfrac) 
      {
        grid[cell].ml.landfrac[1].biomass_grass = tinyfrac;
        grid[cell].ml.irrig_system->biomass_grass = grid[cell].ml.manage.par->default_irrig_system;
      }
      if (grid[cell].ml.landfrac[0].biomass_grass < tinyfrac) grid[cell].ml.landfrac[0].biomass_grass = tinyfrac;
    }

    /* force tinyfrac for all crops only on pixels with valid soil */
    if (landuse->onlygrass && !grid[cell].skip && soiltype!=ROCK && soiltype!=ICE && soiltype >= 0)
    {
      for(j=0; j<ncft; j++)
      {
        grid[cell].ml.landfrac[0].crop[j]=0;
        grid[cell].ml.landfrac[1].crop[j]=0;
        grid[cell].ml.irrig_system->crop[j]=NOIRRIG;
      }
      for(j=0; j<NGRASS; j++)
      {
        grid[cell].ml.landfrac[0].grass[j]=0;
        grid[cell].ml.landfrac[1].grass[j]=0;
        grid[cell].ml.irrig_system->grass[j]=NOIRRIG;
      }
      grid[cell].ml.landfrac[0].biomass_grass=0;
      grid[cell].ml.landfrac[1].biomass_grass=0;
      grid[cell].ml.irrig_system->biomass_grass=NOIRRIG;
      grid[cell].ml.landfrac[0].biomass_tree=0;
      grid[cell].ml.landfrac[1].biomass_tree=0;
      grid[cell].ml.irrig_system->biomass_tree=NOIRRIG;

      if(config->irrig_scenario==NO_IRRIGATION)
      {
        grid[cell].ml.landfrac[0].grass[1]=1;
      }
      else
      {
        grid[cell].ml.landfrac[1].grass[1]=1;
        grid[cell].ml.irrig_system->grass[j]=SPRINK;
      }
    }


    sum=landfrac_sum(grid[cell].ml.landfrac,ncft,FALSE)+landfrac_sum(grid[cell].ml.landfrac,ncft,TRUE);

    /* set landuse to zero if no valid soil */
    if ((grid[cell].skip || soiltype==ROCK || soiltype==ICE || soiltype < 0) && sum>0)
      {
        //fprintf(stderr,"WARNING!! setting LU (sum:%g) to zero, because of invalid soil type %d (%g/%g) in cell %d at year %d\n",
        //        sum,soiltype,grid[cell].coord.lon,grid[cell].coord.lat, cell+config->startgrid,yearl+landuse->landuse.firstyear);
        for(j=0; j<ncft; j++)
        {
          grid[cell].ml.landfrac[0].crop[j]/=sum;
          grid[cell].ml.landfrac[1].crop[j]/=sum;
        }
        for(j=0; j<NGRASS; j++)
        {
          grid[cell].ml.landfrac[0].grass[j]=0;
          grid[cell].ml.landfrac[1].grass[j]=0;
        }
        grid[cell].ml.landfrac[0].biomass_grass=0;
        grid[cell].ml.landfrac[1].biomass_grass=0;
        grid[cell].ml.landfrac[0].biomass_tree=0;
        grid[cell].ml.landfrac[1].biomass_tree=0;
        sum /= sum;
      }

      if(sum>1.00001)
      {
        if(yearl>0&&sum>1.01)
        {
          fprintf(stderr,"WARNING013: in cell %d at year %d: sum of crop fractions greater 1: %f\n",
            cell+config->startgrid,yearl+landuse->landuse.firstyear,sum);
          fflush(stderr);
        }
        sum=reducelanduse(grid+cell,sum-1,ncft);
        if(sum>0.00001)
          fail(CROP_FRACTION_ERR,FALSE,
            "crop fraction greater 1: %f cell: %d, managed grass is 0",
            sum+1,cell+config->startgrid);
      }
      if(landuse->onlycrops) {
        sum=0;
        for(j=0; j<ncft; j++)
        {
          sum+=grid[cell].ml.landfrac[0].crop[j];
          sum+=grid[cell].ml.landfrac[1].crop[j];
        }
        if(sum < 1&&sum > epsilon) {
          for(j=0; j<ncft; j++)
          {
            grid[cell].ml.landfrac[0].crop[j]/=sum;
            grid[cell].ml.landfrac[1].crop[j]/=sum;
          }
          for(j=0; j<NGRASS; j++)
          {
            grid[cell].ml.landfrac[0].grass[j]=0;
            grid[cell].ml.landfrac[1].grass[j]=0;
          }
          grid[cell].ml.landfrac[0].biomass_grass=0;
          grid[cell].ml.landfrac[1].biomass_grass=0;
          grid[cell].ml.landfrac[0].biomass_tree=0;
          grid[cell].ml.landfrac[1].biomass_tree=0;
        }
      }
/** temporary set everything to irrigated maize */
/*        for(j=0; j<ncft; j++)
        {
          grid[cell].ml.landfrac[0].crop[j]=0;
          grid[cell].ml.landfrac[1].crop[j]=0;
        }
        for(j=0; j<NGRASS; j++)
        {
          grid[cell].ml.landfrac[0].grass[j]=0;
          grid[cell].ml.landfrac[1].grass[j]=0;
        }
        grid[cell].ml.landfrac[0].biomass_grass=0;
        grid[cell].ml.landfrac[1].biomass_grass=0;
        grid[cell].ml.landfrac[0].biomass_tree=0;
        grid[cell].ml.landfrac[1].biomass_tree=0;
        grid[cell].ml.landfrac[1].crop[0]=0.5;
        grid[cell].ml.landfrac[1].crop[2]=0.5;
        sum=1;*/
   // }
  } /* for(cell=0;...) */
  free(data);
  if(config->with_nitrogen)
  {

    for(cell=0; cell<config->ngridcell; cell++)
      for(i=0; i<WIRRIG; i++)
      {
        for(j=0; j<ncft; j++)
        {
          grid[cell].ml.fertilizer_nr[i].crop[j]=0;
          grid[cell].ml.manure_nr[i].crop[j]=0;
        }
        for(j=0; j<NGRASS; j++)
        {
          grid[cell].ml.fertilizer_nr[i].grass[j]=0;
          grid[cell].ml.manure_nr[i].grass[j]=0;
        }
        grid[cell].ml.fertilizer_nr[i].biomass_grass=0;
        grid[cell].ml.fertilizer_nr[i].biomass_grass=0;
        grid[cell].ml.manure_nr[i].biomass_tree=0;
        grid[cell].ml.manure_nr[i].biomass_tree=0;
      }

    if(config->fertilizer_input&&!config->fix_fertilization)
    {
      /* assigning fertilizer Nr data */
      yearf-=landuse->fertilizer_nr.firstyear;
      if(yearf>=landuse->fertilizer_nr.nyear)
        yearf=landuse->fertilizer_nr.nyear-1;
      else if(yearf<0)
        yearf=0;
      if(landuse->fertilizer_nr.fmt==CDF)
      {
        fert_nr=newvec(Real,config->ngridcell*landuse->fertilizer_nr.var_len);
        if(fert_nr==NULL)
        {
          printallocerr("fert_nr");
          return TRUE;
        }
        if(readdata_netcdf(&landuse->fertilizer_nr,fert_nr,grid,yearf,config))
        {
          fprintf(stderr,
            "ERROR149: Cannot read fertilizer of year %d in getlanduse().\n",
            yearf+landuse->fertilizer_nr.firstyear);
          fflush(stderr);
          return TRUE;
        }
      }
      else
      {
        if(fseek(landuse->fertilizer_nr.file,(long long)yearf*landuse->fertilizer_nr.size+landuse->fertilizer_nr.offset,SEEK_SET))
        {
          fprintf(stderr,
            "ERROR148: Cannot seek fertilizer Nr to year %d in getlanduse().\n",
            yearf+landuse->fertilizer_nr.firstyear);
          fflush(stderr);
          return TRUE;
        }
        fert_nr=newvec(Real,landuse->fertilizer_nr.n);
        //vec=newvec(short,landuse->fertilizer_nr.n);
        if(fert_nr==NULL)
        {
          printallocerr("fert_nr");
          return TRUE;
        }
        if(readrealvec(landuse->fertilizer_nr.file,fert_nr,0,landuse->fertilizer_nr.scalar,landuse->fertilizer_nr.n,
          landuse->fertilizer_nr.swap,landuse->fertilizer_nr.datatype))
        {
          fprintf(stderr,
            "ERROR149: Cannot read fertilizer Nr of year %d in getlanduse().\n",
            yearf+landuse->fertilizer_nr.firstyear);
          fflush(stderr);
          free(fert_nr);
          return TRUE;
        }
      }
      count=0;

      /* do changes here for the fertilization*/
      for(cell=0; cell<config->ngridcell; cell++)
      {
        for(i=0; i<WIRRIG; i++)
        {
          for(j=0; j<ncft; j++)
            grid[cell].ml.fertilizer_nr[i].crop[j]=fert_nr[count++];
          for(j=0; j<NGRASS; j++)
            grid[cell].ml.fertilizer_nr[i].grass[j]=fert_nr[count++];

          if(landuse->fertilizer_nr.var_len!=2*(ncft+NGRASS))
          {
            grid[cell].ml.fertilizer_nr[i].biomass_grass=fert_nr[count++];
            grid[cell].ml.fertilizer_nr[i].biomass_tree=fert_nr[count++];
          }
          else
            grid[cell].ml.fertilizer_nr[i].biomass_grass=grid[cell].ml.fertilizer_nr[i].biomass_tree=0;
        }
      } /* for(cell=0;...) */
      free(fert_nr);
    }

    if(config->manure_input&&!config->fix_fertilization)
    {
      /* assigning manure fertilizer nr data */
      yearm-=landuse->manure_nr.firstyear;
      if(yearm>=landuse->manure_nr.nyear)
        yearm=landuse->manure_nr.nyear-1;
      else if(yearm<0)
        yearm=0;
      if(landuse->manure_nr.fmt==CDF)
      {
        manu_nr=newvec(Real,config->ngridcell*landuse->manure_nr.var_len);
        if(manu_nr==NULL)
        {
          printallocerr("manu_nr");
          return TRUE;
        }
        if(readdata_netcdf(&landuse->manure_nr,manu_nr,grid,yearm,config))
        {
          fprintf(stderr,
            "ERROR149: Cannot read manure fertilizer of year %d in getlanduse().\n",
            yearm+landuse->manure_nr.firstyear);
          fflush(stderr);
          return TRUE;
        }
      }
      else
      {
        if(fseek(landuse->manure_nr.file,(long long)yearm*landuse->manure_nr.size+landuse->manure_nr.offset,SEEK_SET))
        {
          fprintf(stderr,
            "ERROR148: Cannot seek manure fertilizer to year %d in getlanduse().\n",
            yearm+landuse->manure_nr.firstyear);
          fflush(stderr);
          return TRUE;
        }
        manu_nr=newvec(Real,landuse->manure_nr.n);
        if(manu_nr==NULL)
        {
          printallocerr("manu_nr");
          return TRUE;
        }
        if(readrealvec(landuse->manure_nr.file,manu_nr,0,landuse->manure_nr.scalar,landuse->manure_nr.n,
          landuse->manure_nr.swap,landuse->manure_nr.datatype))
        //if(fread(vec,sizeof(short),landuse->manure_nr.n,landuse->manure_nr.file)!=landuse->manure_nr.n)
        {
          fprintf(stderr,
            "ERROR149: Cannot read manure fertilizer of year %d in getlanduse().\n",
            yearm+landuse->manure_nr.firstyear);
          fflush(stderr);
          free(manu_nr);
          return TRUE;
        }
      }
      count=0;

      /* do changes here for the manure*/
      for(cell=0; cell<config->ngridcell; cell++)
      {
        for(i=0; i<WIRRIG; i++)
        {
          for(j=0; j<ncft; j++)
            grid[cell].ml.manure_nr[i].crop[j]=manu_nr[count++];
          for(j=0; j<NGRASS; j++)
            grid[cell].ml.manure_nr[i].grass[j]=manu_nr[count++];

          if(landuse->manure_nr.var_len!=2*(ncft+NGRASS))
          {
            grid[cell].ml.manure_nr[i].biomass_grass=manu_nr[count++];
            grid[cell].ml.manure_nr[i].biomass_tree=manu_nr[count++];
          }
          else
            grid[cell].ml.manure_nr[i].biomass_grass=grid[cell].ml.manure_nr[i].biomass_tree=0;
        }
      } /* for(cell=0;...) */
      free(manu_nr);
    }

    if(config->fix_fertilization)
    {
      for(cell=0; cell<config->ngridcell; cell++)
      {
        for(i=0; i<WIRRIG; i++)
        {
          for(j=0; j<ncft; j++){
            grid[cell].ml.fertilizer_nr[i].crop[j]=param.fertilizer_rate;
            grid[cell].ml.manure_nr[i].crop[j]=param.manure_rate;
          }
          for(j=0; j<NGRASS; j++){
            grid[cell].ml.fertilizer_nr[i].grass[j]=param.fertilizer_rate;
            grid[cell].ml.manure_nr[i].grass[j]=param.manure_rate;
          }
          grid[cell].ml.fertilizer_nr[i].biomass_grass=param.fertilizer_rate;
          grid[cell].ml.fertilizer_nr[i].biomass_tree=param.fertilizer_rate;
          grid[cell].ml.manure_nr[i].biomass_grass=param.manure_rate;
          grid[cell].ml.manure_nr[i].biomass_tree=param.manure_rate;
        }
      }
    }

  }

  if(config->tillage_type==READ_TILLAGE)
  {
    /* read in tillage data */
    yeart-=landuse->with_tillage.firstyear;
    if(yeart>=landuse->with_tillage.nyear)
      yeart=landuse->with_tillage.nyear-1;
    else if(yeart<0)
      yeart=0;
    if(landuse->with_tillage.fmt==CDF)
    {
      tilltypes=newvec(int,config->ngridcell*landuse->with_tillage.var_len);
      if(tilltypes==NULL)
      {
        printallocerr("tilltypes");
        return TRUE;
      }
      if(readintdata_netcdf(&landuse->with_tillage,tilltypes,grid,yeart,config))
      {
        fprintf(stderr,
          "ERROR149: Cannot read tillage types of year %d in getlanduse().\n",
          year+landuse->with_tillage.firstyear);
        fflush(stderr);
        return TRUE;
      }
      count=0;
      for(cell=0; cell<config->ngridcell; cell++)
        if(!grid[cell].skip)
          grid[cell].ml.with_tillage=tilltypes[count++];
        else
          count+=1;
      free(tilltypes);
    }
    else
    {
      if(fseek(landuse->with_tillage.file,(long long)yeart*landuse->with_tillage.size+landuse->with_tillage.offset,SEEK_SET))
      {
        fprintf(stderr,
          "ERROR148: Cannot seek tillage types to year %d in getlanduse().\n",
          yeart+landuse->with_tillage.firstyear);
        return TRUE;
      }
      vec=newvec(short,landuse->with_tillage.n);
      if(vec==NULL)
      {
        printallocerr("vec2");
        return TRUE;
      }
      if(fread(vec,sizeof(short),landuse->with_tillage.n,landuse->with_tillage.file)!=landuse->with_tillage.n)
      {
        fprintf(stderr,
          "ERROR149: Cannot read tillage types of year %d in getlanduse().\n",
          yeart+landuse->with_tillage.firstyear);
        free(vec);
        return TRUE;
      }
      if(landuse->with_tillage.swap)
        for(i=0; i<landuse->with_tillage.n; i++)
          vec[i]=swapshort(vec[i]);
      count=0;
      for(cell=0; cell<config->ngridcell; cell++)
        if(!grid[cell].skip)
          grid[cell].ml.with_tillage=(Bool)vec[count++]*landuse->with_tillage.scalar;
        else
          count+=1;
      free(vec);
    }
  }
  if(config->tillage_type!=READ_TILLAGE)
  {
    for(cell=0; cell<config->ngridcell; cell++)
    {
      grid[cell].ml.with_tillage=config->tillage_type==NO_TILLAGE ? FALSE : TRUE;
    }
  }

  if(config->residue_treatment==READ_RESIDUE_DATA)
  {
    /* assigning residue extraction data */
    yearr-=landuse->residue_on_field.firstyear;
    if(yearr>=landuse->residue_on_field.nyear)
      yearr=landuse->residue_on_field.nyear-1;
    else if(yearr<0)
      yearr=0;
    if(landuse->residue_on_field.fmt==CDF)
    {
      res_on_field=newvec(Real,config->ngridcell*landuse->residue_on_field.var_len);
      if(res_on_field==NULL)
      {
        printallocerr("residue_on_field");
        return TRUE;
      }
      if(readdata_netcdf(&landuse->residue_on_field,res_on_field,grid,yearr,config))
      {
        fprintf(stderr,
          "ERROR149: Cannot read residue extract of year %d in getlanduse().\n",
          yearr+landuse->residue_on_field.firstyear);
        fflush(stderr);
        return TRUE;
      }
    }
    else
    {
      if(fseek(landuse->residue_on_field.file,(long long)yearr*landuse->residue_on_field.size+landuse->residue_on_field.offset,SEEK_SET))
      {
        fprintf(stderr,
          "ERROR148: Cannot seek residue extraction to year %d in getlanduse().\n",
          yearr+landuse->residue_on_field.firstyear);
        fflush(stderr);
        return TRUE;
      }
      vec=newvec(short,landuse->residue_on_field.n);
      if(vec==NULL)
      {
        printallocerr("vec");
        return TRUE;
      }
      if(fread(vec,sizeof(short),landuse->residue_on_field.n,landuse->residue_on_field.file)!=landuse->residue_on_field.n)
      {
        fprintf(stderr,
          "ERROR149: Cannot read residue extract of year %d in getlanduse().\n",
          yearr+landuse->residue_on_field.firstyear);
        fflush(stderr);
        free(vec);
        return TRUE;
      }
      if(landuse->residue_on_field.swap)
        for(i=0; i<landuse->residue_on_field.n; i++)
          vec[i]=swapshort(vec[i]);
    }
    count=0;

    /* do changes for residue rate left on field*/
    for(cell=0; cell<config->ngridcell; cell++)
    {
      if(landuse->residue_on_field.fmt==CDF)
      {
        for(i=0; i<WIRRIG; i++)
        {
          for(j=0; j<ncft; j++)
            grid[cell].ml.residue_on_field[i].crop[j]=res_on_field[count++];
          for(j=0; j<NGRASS; j++)
            grid[cell].ml.residue_on_field[i].grass[j]=res_on_field[count++];

          if(landuse->residue_on_field.var_len!=ncft+NGRASS)
          {
            grid[cell].ml.residue_on_field[i].biomass_grass=res_on_field[count++];
            grid[cell].ml.residue_on_field[i].biomass_tree=res_on_field[count++];
          }
          else
            grid[cell].ml.residue_on_field[i].biomass_grass=grid[cell].ml.residue_on_field[i].biomass_tree=0;
        }
      }
      else
      {
        for(j=0; j<ncft; j++)
        {
          for(i=0; i<WIRRIG; i++)
            grid[cell].ml.residue_on_field[i].crop[j]=vec[count]*landuse->residue_on_field.scalar;
          count++;
        }
        for(j=0; j<NGRASS; j++)
        {
          for(i=0; i<WIRRIG; i++)
            grid[cell].ml.residue_on_field[i].grass[j]=vec[count]*landuse->residue_on_field.scalar;
          count++;
        }
        if(landuse->residue_on_field.var_len!=ncft+NGRASS)
        {
          grid[cell].ml.residue_on_field[0].biomass_grass=vec[count]*landuse->residue_on_field.scalar;
          grid[cell].ml.residue_on_field[1].biomass_grass=vec[count++]*landuse->residue_on_field.scalar;
          grid[cell].ml.residue_on_field[0].biomass_tree=vec[count]*landuse->residue_on_field.scalar;
          grid[cell].ml.residue_on_field[1].biomass_tree=vec[count++]*landuse->residue_on_field.scalar;
        }
        else
          for(i=0; i<WIRRIG; i++)
            grid[cell].ml.residue_on_field[i].biomass_grass=grid[cell].ml.residue_on_field[i].biomass_tree=0;
      }
    } /* for(cell=0;...) */
    if(landuse->residue_on_field.fmt==CDF)
      free(res_on_field);
    else
      free(vec);
  }
  return FALSE;
} /* of 'getlanduse' */

Bool getintercrop(const Landuse landuse /**< pointer to landuse data */
)                      /** \return intercropping enabled? (TRUE/FALSE) */
{
  return (landuse==NULL) ? FALSE : landuse->intercrop;
} /* of 'getintercrop' */

void freelanduse(Landuse landuse, /**< pointer to landuse data */
  Bool isroot      /**< task is root task */
)
{
  if(landuse!=NULL)
  {
    if(landuse->landuse.fmt==CDF)
      closeclimate_netcdf(&landuse->landuse,isroot);
    else
      fclose(landuse->landuse.file);
    if(landuse->sdate.file!=NULL)
    {
      if(landuse->sdate.fmt==CDF)
        closeclimate_netcdf(&landuse->sdate,isroot);
      else
        fclose(landuse->sdate.file);
    }
    if(landuse->fertilizer_nr.file!=NULL)
    {
      if(landuse->fertilizer_nr.fmt==CDF)
        closeclimate_netcdf(&landuse->fertilizer_nr,isroot);
      else
        fclose(landuse->fertilizer_nr.file);
    }
    free(landuse);
  }
} /* of 'freelanduse' */
