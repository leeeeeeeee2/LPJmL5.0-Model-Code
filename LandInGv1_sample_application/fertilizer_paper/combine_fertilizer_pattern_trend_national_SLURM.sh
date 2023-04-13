#!/bin/bash
#SBATCH --ntasks=32
#SBATCH --qos=short
#SBATCH --export=ALL
#SBATCH --account=lpjml
#SBATCH --comment="combine_fertilizer_pattern_trend_national"
#SBATCH --output=combine_fertilizer_pattern_trend_national_%j.out
#SBATCH --error=combine_fertilizer_pattern_trend_national_%j.err
#SBATCH --mail-type=END
#SBATCH --job-name=combine_fertilizer_pattern_trend_national

################################################################################
## Copyright (C) 2022 Potsdam Institute for Climate Impact Research (PIK),    ##
## see COPYRIGHT file.                                                        ##
##                                                                            ##
## This file is part of LandInG and licensed under GNU AGPL Version 3 or      ##
## later. See LICENSE file or go to http://www.gnu.org/licenses/              ##
## Contact: https://github.com/PIK-LPJmL/LandInG/                             ##
################################################################################

ulimit -c unlimited
export I_MPI_PMI_LIBRARY=/p/system/slurm/lib/libpmi.so

module load R/3.6.2
module load intel/2018.1
module load geos/3.6.1
module load udunits/2.2.19
module load proj4/5.2.0
module load netcdf-c/4.2.1.1/serial
module load gdal/2.4.0
module load cdo/1.9.6/gnu-threadsafe
module load curl/7.58.0

export R_LIBS=/p/projects/lpjml/R.3.6.2/library # Must match with R module selected above

srun --propagate R --no-save --file=combine_fertilizer_pattern_trend_national.R --silent --slave
