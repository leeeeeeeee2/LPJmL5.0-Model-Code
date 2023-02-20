#################################################################################
##                                                                             ##
##          l  p  j  _  p  a  t  h  s  .  c  s  h                              ##
##                                                                             ##
##    csh script to set environment variables for LPJmL. A call to this        ##
##    script has to be put in ~/.cshrc in the following way:                   ##
##                                                                             ##
##    source $LPJROOT/bin/lpj_paths.csh                                        ##
##                                                                             ##
##    LPJROOT has to be set to your root directory of LPJmL                    ##
## (C) Potsdam Institute for Climate Impact Research (PIK), see COPYRIGHT file ##
## authors, and contributors see AUTHORS file                                  ##
## This file is part of LPJmL and licensed under GNU AGPL Version 3            ##
## or later. See LICENSE file or go to http://www.gnu.org/licenses/            ##
## Contact: https://github.com/PIK-LPJmL/LPJmL                                 ##
##                                                                             ##
##    Created: 23.06.2020                                                      ##
##                                                                             ##
#################################################################################

setenv LPJROOT /home/herzfeld/lpj_git/tillage_carbon/LPJmL_internal # change here to your directory

# set search path for LPJmL commands

setenv PATH $LPJROOT/bin\:$PATH

# set path for input files

setenv LPJINPATH /p/projects/lpjml/input/historical

# include manpages of LPJmL

setenv MANPATH $LPJROOT/man\:$MANPATH

# define alias

alias printheader "printclm -data"
alias lpjml 'lpjml.sh'
