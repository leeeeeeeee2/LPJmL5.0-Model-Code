#################################################################################
##                                                                             ##
##          l  p  j  _  p  a  t  h  s  .  s  h                                 ##
##                                                                             ##
##    sh script to set environment variables for LPJmL. A call to this         ##
##    script has to be put in ~/.profile in the following way:                 ##
##                                                                             ##
##    . $LPJROOT/bin/lpj_paths.sh                                              ##
##                                                                             ##
##    LPJROOT has to be set to your root directory of LPJmL                    ##
##                                                                             ##
## (C) Potsdam Institute for Climate Impact Research (PIK), see COPYRIGHT file ##
## authors, and contributors see AUTHORS file                                  ##
## This file is part of LPJmL and licensed under GNU AGPL Version 3            ##
## or later. See LICENSE file or go to http://www.gnu.org/licenses/            ##
## Contact: https://github.com/PIK-LPJmL/LPJmL                                 ##
##                                                                             ##
##    Created: 23.06.2020                                                      ##
##                                                                             ##
#################################################################################

export LPJROOT=/home/herzfeld/lpj_git/tillage_carbon/LPJmL_internal # change here to your directory

# set search path for LPJmL commands

export PATH=$LPJROOT/bin:$PATH

# set path for input files

export LPJINPATH=/p/projects/lpjml/input/historical

# include manpages of LPJmL

export MANPATH=/home/herzfeld/lpj_git/tillage_carbon/LPJmL_internal/man:$MANPATH

# define alias

alias printheader="printclm -data"
alias lpjml='lpjml.sh'
