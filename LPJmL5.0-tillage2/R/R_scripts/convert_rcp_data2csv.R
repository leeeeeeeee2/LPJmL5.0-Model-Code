###############################################################################################
#                                                                                             #
#  Script by Tobias Herzfeld (04.2020) for:                                                   #
# "SOC sequestration potentials for agricultural management practices under climate change"   #
#                                                                                             #
#  Reads in the data created by LPJmL5.0-tillage2 model runs and calculates global sums.      #
#  Output data is stored in csv files to be used further.                                     #
#                                                                                             #
###############################################################################################

rm(list=ls(all=T));gc()

library(oce)
library(Hmisc)
library(weights)
library(plotrix)

##define simulations:
sims<-c("rcp26_hadgem",
        "rcp85_hadgem",
        "rcp26_gfdl",
        "rcp85_gfdl",
        "rcp26_ipsl",
        "rcp85_ipsl",
        "rcp26_miroc",
        "rcp85_miroc")

##define paths:
output_path<-"path/.../"           #path of model output data for rcps
grid_path<-"path/.../R/grid/"      #path for lpj grid
landuse_path<-"path/.../"          #path for landuse input data
plots_path<-"path/.../"            #path for output boxplots and barplots created

##define scens:
scens <- c("T_NR","T_R","NT_NR","NT_R")

##define NCELL and years:
NCELL <- 67420
refyear <- 2006
syear <- 2006
eyear <- 2099

nyears <- length(syear:eyear)
years<-eyear-syear+1

dirnames <- array(dim=length(scens),dimnames=list(scens))

dirnames["T_R"] <- "OUTPUT_T_R"
dirnames["T_NR"] <- "OUTPUT_T_NR"
dirnames["NT_R"] <- "OUTPUT_NT_R"
dirnames["NT_NR"] <- "OUTPUT_NT_NR"

bnames_cft <- c("wheat","rice","maize","millet","pulses","sugarbeet","cassava","sunflower","soybeans","groundnuts","rapeseed","sugarcane","others","grasses","begrasses","betrees")
bbnames_cft<- c(paste0(bnames_cft,"_rf"),paste0(bnames_cft,"_surf"),paste0(bnames_cft,"_sprink"), paste0(bnames_cft,"_drip"))
bnames_cft <- c(paste0(bnames_cft,"_rf"),paste0(bnames_cft,"_ir"))
nbands_cft <- length(bnames_cft)
bnames_pft <- c("TrBEG","TrBRG","TeNEG","TeBEG","TeBSG","BorNEG","BorBSG","BorNSG","C3_GRASS","C4_GRASS",bnames_cft)
nbands_pft <- length(bnames_pft)

cnames_cft <- c("wheat","rice","maize","millet","pulses","sugarbeet","cassava","sunflower","soybeans","groundnuts","rapeseed","sugarcane","others","grasses")
cnames_cft <- c(paste0(cnames_cft,"_rf"),paste0(cnames_cft,"_ir"))
cbands_cft <- length(cnames_cft)

bnames_soil <- c("1","2","3","4","5")
nbands_soil <- length(bnames_soil)

plotcrops <- c("wheat_rf","maize_rf","pulses_rf","rapeseed_rf")
npc <- length(plotcrops)
plotcrops_names <- array(c("rainfed wheat","rainfed maize","rainfed pulses","rainfed rapeseed"),dimnames=list(plotcrops))

################################################################################################################

###read grid + calculate cellara:
zz <- file(paste0(grid_path,"grid.bin"),"rb")
seek(zz,where=43,origin = "start")
NCELL<-67420
x <- readBin(zz, integer(), n=2*NCELL, size=2) / 100
lon <- x[c(1:NCELL)*2-1]
lat <- x[c(1:NCELL)*2]
cellarea <- (111e3*0.5)*(111e3*0.5)*cos(lat/180*pi)
ilon <- as.integer((lon+180)/0.5 + 1.01)
ilat <- as.integer((lat+90)/0.5 + 1.01)
mat.ilonilat <- matrix(c(ilon,ilat),nrow=NCELL,ncol=2,byrow=F)
close(zz);rm(zz,x);gc()

######read landuse data functions:

read_data_landuse_single_year <- function(filename="",startyear=2005, refyear=1700, NCELL=67420, lu=64)
{
  data.out <- array(data=0,dim=c(NCELL,lu)) #
  f <- file(filename,"rb")
  seek(f,where=43+NCELL*(startyear-refyear)*2*lu,origin="start")
  data.out <- matrix(readBin(f,integer(),n=NCELL*lu,size=2), ncol=lu,byrow=T)
  #colnames(data.out) <- bbnames_cft #for 64bands
  colnames(data.out) <- bnames_cft #for 32 bands
  close(f)
  data.out/1000
}

read_landuse_cropland_area_to_years_32 <- function(filename="",startyear=1700, refyear=850, header=43, NCELL=67420, lu=32, NYEARS=years, area=cellarea)
{
  data.out <- array(data=0,dim=c(NCELL,NYEARS))
  f <- file(filename,"rb")
  dummy<-matrix(1:NCELL,ncol=1)
  for (iy in 1:NYEARS){
    dummy2 <- array(data=0,dim=c(NCELL,lu))
    seek(f,where=header+NCELL*((startyear-1+iy)-refyear)*2*lu,origin="start")
    dummy2 <- matrix(readBin(f,integer(),n=NCELL*lu,size=2), ncol=lu,byrow=T)
    dummy2 <- rowSums(dummy2[,-c(rep(14:16,2)+rep(16*(0:1),each=3))]) # sum up cropland bands including others
    dummy <- cbind(dummy,dummy2)
  }
  data.out<-dummy[,2:(NYEARS+1)]
  colnames(data.out)<-c(startyear:(startyear+NYEARS-1))
  close(f);rm(f,dummy,dummy2)
  data.out <- data.out*area/1000      #scaling 1000 for LU
}

###in case of transient LANDUSE:
landuse32bands<-paste0(landuse_path,"lu_madrat_850-2015_32bands.clm")
data.landuse_cropland_area_32<-read_landuse_cropland_area_to_years_32(landuse32bands)

###read landuse for single year + calculate cropfrac:
data_landuse2005_32 <- read_data_landuse_single_year(landuse32bands,startyear=2005,refyear=850,lu=32)
frac_cropland2005_32 <- rowSums(data_landuse2005_32[,-c(rep(14:16,2)+rep(16*(0:1),each=3))])
is_cropland2005_32 <- ifelse(frac_cropland2005_32>0,1,NA)
area_cropland2005_32 <- frac_cropland2005_32*cellarea # in m2

################################################################

###function to read outputs:

read_output_global_years <- function(filename="", area=cellarea, NCELL=67420,syear=1700, refyear=1700, NYEARS=years,header=0,size=4, scaling=1e+15)
{
  data.out <- array(data=0,dim=c(NCELL,NYEARS))
  f <- file(filename,"rb")
  seek(f,where=header+NCELL*(syear-refyear)*size*NYEARS,origin="start")
  data.out<-matrix(readBin(f,double(),n=NCELL*NYEARS,size), ncol=NYEARS,byrow=F)
  colnames(data.out) <- c(syear:(syear+NYEARS-1))
  close(f)
  data.out<-data.out*area/scaling
  colSums(data.out)
}

###create csv:

for (sim in sims) {
  
  outpath<-paste0(output_path,sim,"/")

  for (scen in scens)
  {
    npp_agr<-read_output_global_years(paste0(outpath,"OUTPUT_",scen,"/anpp_agr.bin"),area=area_cropland2005_32,syear=2006,refyear=2006)
    soilc_agr<-read_output_global_years(paste0(outpath,"OUTPUT_",scen,"/soilc_agr.bin"),area=area_cropland2005_32,syear=2006,refyear=2006)
    litc_agr<-read_output_global_years(paste0(outpath,"OUTPUT_",scen,"/litc_agr.bin"),area=area_cropland2005_32,syear=2006,refyear=2006)
    allc_agr<-soilc_agr+litc_agr
    rh_agr<-read_output_global_years(paste0(outpath,"OUTPUT_",scen,"/arh_agr.bin"),area=area_cropland2005_32,syear=2006,refyear=2006)
    litf_agr<-read_output_global_years(paste0(outpath,"OUTPUT_",scen,"/alittfallc_agr.bin"),area=area_cropland2005_32,syear=2006,refyear=2006)
    mrt_agr<-ifelse(rh_agr>0,allc_agr/rh_agr,NA)
    
    npp<-read_output_global_years(paste0(outpath,"OUTPUT_",scen,"/anpp.bin"),syear=2006,refyear=2006)
    soilc<-read_output_global_years(paste0(outpath,"OUTPUT_",scen,"/soilc.bin"),syear=2006,refyear=2006)
    litc<-read_output_global_years(paste0(outpath,"OUTPUT_",scen,"/litc.bin"),syear=2006,refyear=2006)
    vegc<-read_output_global_years(paste0(outpath,"OUTPUT_",scen,"/vegc.bin"),syear=2006,refyear=2006)
    
    allc<-soilc+litc
    allc_veg<-soilc+litc+vegc
    rh<-read_output_global_years(paste0(outpath,"OUTPUT_",scen,"/arh.bin"),syear=2006,refyear=2006)
    df<-data.frame(npp_agr,soilc_agr,litc_agr,allc_agr,rh_agr,mrt_agr,litf_agr,npp,soilc,litc,allc,allc_veg,rh)
    write.csv(df,file=paste0(plots_path,"OUTPUT_",scen,"_df_",sim,".csv"))
    #write.csv(df,file=paste0(plots_path,"OUTPUT_",scen,"_df_",sim,"_CONST_CO2.csv"))
  }
}



