#!/bin/bash

################################################################################
## Copyright (C) 2022 Potsdam Institute for Climate Impact Research (PIK),    ##
## see COPYRIGHT file.                                                        ##
##                                                                            ##
## This file is part of LandInG and licensed under GNU AGPL Version 3 or      ##
## later. See LICENSE file or go to http://www.gnu.org/licenses/              ##
## Contact: https://github.com/PIK-LPJmL/LandInG/                             ##
################################################################################

################################################################################
## This script creates a river routing input file for LPJmL from a drainage   ##
## direction map.                                                             ##
## Parameters:                                                                ##
##   - DRAINAGEBIN: Path to drainage executable from LPJmL utilities.         ##
##     Drainage utility is included in LPJmL source code repository and       ##
##     compiled with "make all".                                              ##
##   - LPJGRID: Path to LPJmL grid input file for which a river routing input ##
##     will be created.                                                       ##
##   - DDMASCII: Path to a drainage direction map as an ASCII grid. This must ##
##     have the same spatial resolution as LPJGRID and cover at least the     ##
##     full spatial extent of LPJGRID. For encoding of drainage directions    ##
##     see /src/utils/drainage.c in LPJmL source code reposity.               ##
##   - RIVERROUTINGFILE: River routing input file to be created.              ##
################################################################################

DRAINAGEBIN="./drainage"
# 5 minute versions
#LPJGRID="../gadm_paper/grid_gadm_5arcmin.bin"
# 5 minute HydroSHEDS v.1.1
#DDMASCII="HydroSHEDS_v1.1/HydroSheds_v1.1.asc"
# 5 minute HydroSHEDS v.1.1
#RIVERROUTINGFILE="drainage_hydrosheds_v1.1_5arcmin.bin"
# 5 minute MERIT
#DDMASCII="MERIT_Hydro_IHU/05min_flwdir.asc"
# 5 minute MERIT
#RIVERROUTINGFILE="drainage_merit_5arcmin.bin"

# 30 minute versions
LPJGRID="../gadm_paper/grid_gadm_30arcmin_predefined_grid.bin"
# 30 minute DDM30
#DDMASCII="DDM30/DDM30.asc"
# 30 minute DDM30
#RIVERROUTINGFILE="drainage_ddm30_30arcmin.bin"
# 30 minute STN
DDMASCII="STN-30/global_30_minute_potential_network_v601_asc/g_network.asc"
# 30 minute STN
RIVERROUTINGFILE="drainage_stn_30arcmin.bin"

## Convert drainage direction map into LPJmL river routing input
$DRAINAGEBIN "$LPJGRID" "$DDMASCII" "$RIVERROUTINGFILE"
