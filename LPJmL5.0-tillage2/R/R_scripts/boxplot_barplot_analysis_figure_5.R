###############################################################################################
#                                                                                             #
#  Script by Tobias Herzfeld (04.2020) for:                                                   #
# "SOC sequestration potentials for agricultural management practices under climate change"   #
#                                                                                             #
#  Reads in the data created by LPJmL5.0-tillage2 model runs and creates the boxplot          #
#  and and barplot in figure 5.                                                               #
#                                                                                             #
###############################################################################################

rm(list=ls(all=T));gc()

library(oce)
library(Hmisc)
library(weights)
library(plotrix)

##define paths:
output_path<-"path/.../"           #path of model output data for rcps
lpjinput_path<-"path/.../"         #path for lpj input data
grid_path<-"path/.../R/grid/"      #path for lpj grid
landuse_path<-"path/.../"          #path for landuse input data
mpet_path<-"path/.../"             #path for mpet input
plots_path<-"path/.../"            #path for output boxplots and barplots created

##path of model output data for rcps:
outpath_h1<-paste0(output_path,"rcp26_hadgem/")
outpath_h2<-paste0(output_path,"rcp85_hadgem/")
outpath_g1<-paste0(output_path,"rcp26_gfdl/")
outpath_g2<-paste0(output_path,"rcp85_gfdl/")
outpath_i1<-paste0(output_path,"rcp26_ipsl/")
outpath_i2<-paste0(output_path,"rcp85_ipsl/")
outpath_m1<-paste0(output_path,"rcp26_miroc/")
outpath_m2<-paste0(output_path,"rcp85_miroc/")

#define scens:
scens <- c("T_R","NT_R","T_NR","NT_NR")

#define NCELL and years:
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

##read grid + calculate cellara:

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

##read landuse data functions:

read_data_landuse_single_year <- function(filename="",startyear=2005, refyear=1700, NCELL=67420, lu=64)
{
  data.out <- array(data=0,dim=c(NCELL,lu)) #
  f <- file(filename,"rb")
  seek(f,where=43+NCELL*(startyear-refyear)*2*lu,origin="start")
  data.out <- matrix(readBin(f,integer(),n=NCELL*lu,size=2), ncol=lu,byrow=T)
  colnames(data.out) <- bnames_cft #for 32 bands
  close(f)
  data.out/1000
}

##read landuse for single year + calculate cropfrac:

landuse32bands<-paste0(landuse_path,"lu_madrat_850-2015_32bands.clm")
data_landuse2005_32 <- read_data_landuse_single_year(landuse32bands,startyear=2005,refyear=850,lu=32)
frac_cropland2005_32 <- rowSums(data_landuse2005_32[,-c(rep(14:16,2)+rep(16*(0:1),each=3))])
is_cropland2005_32 <- ifelse(frac_cropland2005_32>0,1,0)
area_cropland2005_32 <- frac_cropland2005_32*cellarea # in m2

################################################################

###funtions read outputs:

read_data_annual <- function(filename="",outp=outpath,syr=syear,eyr=eyear,ryr=refyear)
{
  nyears=length(syr:eyr)
  data.out <- array(data=0,dim=c(NCELL,length(scens)),dimnames=list(NULL,scens))
  for(scen in scens)
  {
    f <- file(paste0(outp,"/",dirnames[scen],"/",filename),"rb")
    seek(f,where=NCELL*(syr-ryr)*4,origin="start")
    for(iy in 1:nyears) data.out[,scen] <- data.out[,scen] + readBin(f,double(),n=NCELL,size=4)
    close(f)
  }
  data.out/nyears
}

read_data_monthly <- function(filename="",outp=outpath,syr=syear,ryr=refyear,nyrs=nyears)
{
  data.out <- array(dim=c(NCELL,length(scens)),dimnames=list(NULL,scens))
  for(scen in scens)
  {
    f <- file(paste0(outp,dirnames[scen],"/",filename),"rb")
    seek(f,where=NCELL*(syr-ryr)*12*4,origin="start")
    dummy <- numeric(NCELL*12)
    for(iy in 1:nyrs) dummy <- dummy + readBin(f,double(),n=NCELL*12,size=4)
    data.out[,scen] <- apply(array(dummy/nyrs,dim=c(NCELL,12)),1,sum)
    close(f);rm(f,dummy);gc()
  }
  data.out
}

read_data_monthly_onepath <- function(filename="",syr=syear,ryr=refyear,nyrs=nyears)
{
  data.out <- array(dim=c(NCELL,length(scens)),dimnames=list(NULL,scens))
  f <- file(filename,"rb")
  seek(f,where=NCELL*(syr-ryr)*12*4,origin="start")
  dummy <- numeric(NCELL*12)
  for(iy in 1:nyrs) dummy <- dummy + readBin(f,double(),n=NCELL*12,size=4)
  data.out <- apply(array(dummy/nyrs,dim=c(NCELL,12)),1,sum)
  close(f);rm(f,dummy);gc()
  data.out
}


read_prec <- function(filename="", header.size=43, scaling=1, data.size=2, cols=12, startyr=syear, refyr=refyear, ncell=NCELL, years=nyears)
{
  data.out <- array(data=0,dim=c(ncell,cols))
  f <- file(filename,"rb")
  seek(f,where=header.size+ncell*(startyr-refyr)*data.size*cols,origin="start")
  for (iy in 1:years) {
    datayear <- matrix(readBin(f,integer(),n=ncell*cols,size=data.size), ncol=cols,byrow=T)
    data.out <- data.out + datayear
  }
  data.out <- rowSums(data.out)
  close(f)
  data.out/scaling/years
}

read_climate_monthly <- function(filename="", header.size=43, scaling=1, data.size=2, cols=12, startyr=syear, refyr=refyear, ncell=NCELL, years=nyears)
{
  data.out <- array(data=0,dim=c(ncell,cols))
  f <- file(filename,"rb")
  seek(f,where=header.size+ncell*(startyr-refyr)*data.size*cols,origin="start")
  for (iy in 1:years) {
    datayear <- matrix(readBin(f,integer(),n=ncell*cols,size=data.size), ncol=cols,byrow=T)
    data.out <- data.out + datayear
  }
  data.out <- rowSums(data.out)/cols
  close(f)
  data.out/scaling/years
}

#Climate zones
#categorising data into different climates after IPCC 2005

#CRU input:
temp<- "cru_ts3.21.1901.2012.tmp.clm" # use if readinputmonthly
temp<- paste0(lpjinput_path, temp) #use if read clim

temperature<- read_climate_monthly(temp,startyr=2000,refyr=1901,years=10,header=43,scaling=10)

#read in precipitation
prec<- "gpcc_v6_cruts3_21_precip_1901_2010.clm" # use if readinputmonthly
prec<- paste0(lpjinput_path, prec)#use if read clim
precipitation<-read_prec(prec, startyr=2000,refyr=1901,years=10, header=43, scaling=1)
prec_trop_dry<-prec_trop_moist<-prec_temperate_moist<-prec_trop_dry<-precipitation

#read in PET
PETdata<-paste0(mpet_path,"mpet.bin")
PET<- read_data_monthly_onepath(PETdata,syr=2000,ryr=1700,nyrs=10) #num [1:67420, 1:4]

index<-precipitation/PET

#index as in IPCC (2006)
index_arid<-index
index_arid[index<1]<-1
index_arid[index>=1]<-0

index_humid<-index
index_humid[index<1]<-0
index_humid[index>=1]<-1

tropical_moist<-tropical_dry<-tropical_wet<-precipitation
tropical_moist[precipitation>1000]<-1
tropical_moist[precipitation<=1000]<-0
tropical_moist[precipitation>2000]<-0
tropical_dry[precipitation<=1000]<-1
tropical_dry[precipitation>1000]<-0
tropical_wet[precipitation>2000]<-1
tropical_wet[precipitation<=2000]<-0

#temperatures
t_cold_temperate<-t_warm_temperate<-t_tropical<-t_boreal<-temperature
t_tropical[temperature<19]<-0
t_tropical[temperature>=19]<-1

t_cold_temperate[temperature<11]<-1
t_cold_temperate[temperature>=11]<-0
t_cold_temperate[temperature<0]<-0
t_warm_temperate[temperature<=10]<-0
t_warm_temperate[temperature>10]<-1
t_warm_temperate[temperature>18]<-0
t_boreal[temperature<=0]<-1
t_boreal[temperature>0]<-0

is_cold_temperate_dry<-is_cropland2005_32*t_cold_temperate*index_arid
is_cold_temperate_moist<-is_cropland2005_32*t_cold_temperate*index_humid
is_warm_temperate_dry<-is_cropland2005_32*t_warm_temperate*index_arid
is_warm_temperate_moist<-is_cropland2005_32*t_warm_temperate*index_humid

is_tropical_dry<-is_cropland2005_32*t_tropical*tropical_dry
is_tropical_moist<-is_cropland2005_32*t_tropical*tropical_moist
is_tropical_wet<-is_cropland2005_32*t_tropical*tropical_wet

is_boreal_dry<-is_cropland2005_32*t_boreal*index_arid
is_boreal_moist<-is_cropland2005_32*t_boreal*index_humid

area_tropical_dry<-(sum(area_cropland2005_32*is_tropical_dry))/1000000            #2852093  #3
area_tropical_moist<-(sum(area_cropland2005_32*is_tropical_moist))/1000000        #2805863  #4
area_tropical_wet<-(sum(area_cropland2005_32*is_tropical_wet))/1000000            #861803   #5
area_warmtemp_dry<-(sum(area_cropland2005_32*is_warm_temperate_dry))/1000000      #3439034  #2
area_warmtemp_moist<-(sum(area_cropland2005_32*is_warm_temperate_moist))/1000000  #657880   #6
area_coldtemp_dry<-(sum(area_cropland2005_32*is_cold_temperate_dry))/1000000      #3857076  #1
area_coldtemp_moist<-(sum(area_cropland2005_32*is_cold_temperate_moist))/1000000  #494375   #7
area_boreal_dry<-(sum(area_cropland2005_32*is_boreal_dry))/1000000                #85102    #8
area_boreal_moist<-(sum(area_cropland2005_32*is_boreal_moist))/1000000            #6353     #9

#area in km2
# area_arid<-8787269
# area_tropical<-6289651
# area_humid<-6157857
# area_coldtemp<-4754397
# area_warmtemp<-3901134

#######################################################

d_allc99_06_h1<-((read_data_annual("soilc_agr.bin",outp=outpath_h1,syr=2099,eyr=2099,ryr=2006)+
                    read_data_annual("litc_agr.bin",outp=outpath_h1,syr=2099,eyr=2099,ryr=2006)+
                    read_data_annual("soilc_agr.bin",outp=outpath_g1,syr=2099,eyr=2099,ryr=2006)+
                    read_data_annual("litc_agr.bin",outp=outpath_g1,syr=2099,eyr=2099,ryr=2006)+
                    read_data_annual("soilc_agr.bin",outp=outpath_i1,syr=2099,eyr=2099,ryr=2006)+
                    read_data_annual("litc_agr.bin",outp=outpath_i1,syr=2099,eyr=2099,ryr=2006)+
                    read_data_annual("soilc_agr.bin",outp=outpath_m1,syr=2099,eyr=2099,ryr=2006)+
                    read_data_annual("litc_agr.bin",outp=outpath_m1,syr=2099,eyr=2099,ryr=2006))/4-
                   (read_data_annual("soilc_agr.bin",outp=outpath_h1,syr=2006,eyr=2006,ryr=2006)+
                      read_data_annual("litc_agr.bin",outp=outpath_h1,syr=2006,eyr=2006,ryr=2006)+
                      read_data_annual("soilc_agr.bin",outp=outpath_g1,syr=2006,eyr=2006,ryr=2006)+
                      read_data_annual("litc_agr.bin",outp=outpath_g1,syr=2006,eyr=2006,ryr=2006)+
                      read_data_annual("soilc_agr.bin",outp=outpath_i1,syr=2006,eyr=2006,ryr=2006)+
                      read_data_annual("litc_agr.bin",outp=outpath_i1,syr=2006,eyr=2006,ryr=2006)+
                      read_data_annual("soilc_agr.bin",outp=outpath_m1,syr=2006,eyr=2006,ryr=2006)+
                      read_data_annual("litc_agr.bin",outp=outpath_m1,syr=2006,eyr=2006,ryr=2006))/4)/1000 #convert to kg/m2

d_trop_wet_h1<-ifelse(d_allc99_06_h1*is_tropical_wet==0,NA,d_allc99_06_h1*is_tropical_wet)
d_trop_moist_h1<-ifelse(d_allc99_06_h1*is_tropical_moist==0,NA,d_allc99_06_h1*is_tropical_moist)
d_trop_dry_h1<-ifelse(d_allc99_06_h1*is_tropical_dry==0,NA,d_allc99_06_h1*is_tropical_dry)
d_warmt_moist_h1<-ifelse(d_allc99_06_h1*is_warm_temperate_moist==0,NA,d_allc99_06_h1*is_warm_temperate_moist)
d_warmt_dry_h1<-ifelse(d_allc99_06_h1*is_warm_temperate_dry==0,NA,d_allc99_06_h1*is_warm_temperate_dry)
d_coldt_moist_h1<-ifelse(d_allc99_06_h1*is_cold_temperate_moist==0,NA,d_allc99_06_h1*is_cold_temperate_moist)
d_coldt_dry_h1<-ifelse(d_allc99_06_h1*is_cold_temperate_dry==0,NA,d_allc99_06_h1*is_cold_temperate_dry)
d_boreal_moist_h1<-ifelse(d_allc99_06_h1*is_boreal_moist==0,NA,d_allc99_06_h1*is_boreal_moist)
d_boreal_dry_h1<-ifelse(d_allc99_06_h1*is_boreal_dry==0,NA,d_allc99_06_h1*is_boreal_dry)

df_d_average_26<-cbind(d_trop_wet_h1,d_trop_moist_h1,d_trop_dry_h1,
                       d_warmt_moist_h1,d_warmt_dry_h1,
                       d_coldt_moist_h1,d_coldt_dry_h1,
                       d_boreal_moist_h1,d_boreal_dry_h1)

df_names<-c("Tropical wet","Tropical moist","Tropical dry",
            "Warmtemperate moist","Warmtemperate dry",
            "Coldtemperate moist","Coldtemperate dry",
            "Boreal moist","Boreal dry")

df_names1<-c("","","Tropical wet","","","","Tropical moist","","","","Tropical dry","",
             "","","Warm temperate moist","","","","Warm temperate dry","",
             "","","Cold temperate moist","","","","Cold temperate dry","",
             "","","Boreal moist","","","","Boreal dry","")

n_names<-c("","n=2088","","","","n=6098","","","","n=5686",
           "","","","n=1514","","","","n=4893","",
           "","","n=2923","","","","n=7068","",
           "","","n=328","","","","n=989","")

myCol<-rep(c("#80cdc1","#018571","#dfc27d","#a6611a"))

pdf(paste0(plots_path,"boxplot_SOC_change_regions_gcm_average_rcp26_notitle_A4.pdf"))
boxplot(df_d_average_26,col=myCol,outline=FALSE,boxwex=0.7,xaxt="n",space=0,las=2,ylim=c(-15,7),whisklty=1,#main="GCM average for RCP2.6",cex.main=0.8,
        ylab="Cropland SOC density change 2006-2099 (kg "~m^-2~")",cex.lab=0.8)
axis(side=1,at=c(4.5,8.5,12.5,16.5,20.5,24.5,28.5,32.5),labels=FALSE)
abline(0,0,lty=1)
abline(v=4.5,lty=2)
abline(v=8.5,lty=2)
abline(v=12.5,lty=2)
abline(v=16.5,lty=2)
abline(v=20.5,lty=2)
abline(v=24.5,lty=2)
abline(v=28.5,lty=2)
abline(v=32.5,lty=2)
text(x=1:36+1.5,y=par("usr")[3]-0.5,labels=df_names1,xpd=NA,srt=25,adj=1.1,cex=0.9)
text(x=1:36+1.5,y=6.6,labels=n_names,xpd=NA,srt=0,adj=1,cex=0.75)
text(x=1,y=8.5,labels="(A)",xpd=NA)
legend("bottomleft",legend=c("T_R", "NT_R", "T_NR","NT_NR"),pch=c(15,15,15,15),col=c("#80cdc1","#018571","#dfc27d","#a6611a"),bty="n",cex=0.85)
dev.off()

d_allc99_06_h2<-((read_data_annual("soilc_agr.bin",outp=outpath_h2,syr=2099,eyr=2099,ryr=2006)+
                    read_data_annual("litc_agr.bin",outp=outpath_h2,syr=2099,eyr=2099,ryr=2006)+
                    read_data_annual("soilc_agr.bin",outp=outpath_g2,syr=2099,eyr=2099,ryr=2006)+
                    read_data_annual("litc_agr.bin",outp=outpath_g2,syr=2099,eyr=2099,ryr=2006)+
                    read_data_annual("soilc_agr.bin",outp=outpath_i2,syr=2099,eyr=2099,ryr=2006)+
                    read_data_annual("litc_agr.bin",outp=outpath_i2,syr=2099,eyr=2099,ryr=2006)+
                    read_data_annual("soilc_agr.bin",outp=outpath_m2,syr=2099,eyr=2099,ryr=2006)+
                    read_data_annual("litc_agr.bin",outp=outpath_m2,syr=2099,eyr=2099,ryr=2006))/4-
                   (read_data_annual("soilc_agr.bin",outp=outpath_h2,syr=2006,eyr=2006,ryr=2006)+
                      read_data_annual("litc_agr.bin",outp=outpath_h2,syr=2006,eyr=2006,ryr=2006)+
                      read_data_annual("soilc_agr.bin",outp=outpath_g2,syr=2006,eyr=2006,ryr=2006)+
                      read_data_annual("litc_agr.bin",outp=outpath_g2,syr=2006,eyr=2006,ryr=2006)+
                      read_data_annual("soilc_agr.bin",outp=outpath_i2,syr=2006,eyr=2006,ryr=2006)+
                      read_data_annual("litc_agr.bin",outp=outpath_i2,syr=2006,eyr=2006,ryr=2006)+
                      read_data_annual("soilc_agr.bin",outp=outpath_m2,syr=2006,eyr=2006,ryr=2006)+
                      read_data_annual("litc_agr.bin",outp=outpath_m2,syr=2006,eyr=2006,ryr=2006))/4)/1000 #convert to kg/m2

d_trop_wet_h2<-ifelse(d_allc99_06_h2*is_tropical_wet==0,NA,d_allc99_06_h2*is_tropical_wet)
d_trop_moist_h2<-ifelse(d_allc99_06_h2*is_tropical_moist==0,NA,d_allc99_06_h2*is_tropical_moist)
d_trop_dry_h2<-ifelse(d_allc99_06_h2*is_tropical_dry==0,NA,d_allc99_06_h2*is_tropical_dry)
d_warmt_moist_h2<-ifelse(d_allc99_06_h2*is_warm_temperate_moist==0,NA,d_allc99_06_h2*is_warm_temperate_moist)
d_warmt_dry_h2<-ifelse(d_allc99_06_h2*is_warm_temperate_dry==0,NA,d_allc99_06_h2*is_warm_temperate_dry)
d_coldt_moist_h2<-ifelse(d_allc99_06_h2*is_cold_temperate_moist==0,NA,d_allc99_06_h2*is_cold_temperate_moist)
d_coldt_dry_h2<-ifelse(d_allc99_06_h2*is_cold_temperate_dry==0,NA,d_allc99_06_h2*is_cold_temperate_dry)
d_boreal_moist_h2<-ifelse(d_allc99_06_h2*is_boreal_moist==0,NA,d_allc99_06_h2*is_boreal_moist)
d_boreal_dry_h2<-ifelse(d_allc99_06_h2*is_boreal_dry==0,NA,d_allc99_06_h2*is_boreal_dry)

df_d_average_85<-cbind(d_trop_wet_h2,d_trop_moist_h2,d_trop_dry_h2,
                       d_warmt_moist_h2,d_warmt_dry_h2,
                       d_coldt_moist_h2,d_coldt_dry_h2,
                       d_boreal_moist_h2,d_boreal_dry_h2)

df_names<-c("Tropical wet","Tropical moist","Tropical dry",
            "Warmtemperate moist","Warmtemperate dry",
            "Coldtemperate moist","Coldtemperate dry",
            "Boreal moist","Boreal dry")

df_names1<-c("","","Tropical wet","","","","Tropical moist","","","","Tropical dry","",
             "","","Warm temperate moist","","","","Warm temperate dry","",
             "","","Cold temperate moist","","","","Cold temperate dry","",
             "","","Boreal moist","","","","Boreal dry","")

n_names<-c("","n=2088","","","","n=6098","","","","n=5686",
           "","","","n=1514","","","","n=4893","",
           "","","n=2923","","","","n=7068","",
           "","","n=328","","","","n=989","")

myCol<-rep(c("#80cdc1","#018571","#dfc27d","#a6611a"))

pdf(paste0(plots_path,"boxplot_SOC_change_regions_gcm_average_rcp85_notitle_C4.pdf"))
boxplot(df_d_average_85,col=myCol,outline=FALSE,boxwex=0.7,xaxt="n",space=0,las=2,ylim=c(-15,7),whisklty=1,#main="GCM average for RCP8.5",cex.main=0.8,
        ylab="Cropland SOC density change 2006-2099 (kg "~m^-2~")",cex.lab=0.8)
axis(side=1,at=c(4.5,8.5,12.5,16.5,20.5,24.5,28.5,32.5),labels=FALSE)
abline(0,0,lty=1)
abline(v=4.5,lty=2)
abline(v=8.5,lty=2)
abline(v=12.5,lty=2)
abline(v=16.5,lty=2)
abline(v=20.5,lty=2)
abline(v=24.5,lty=2)
abline(v=28.5,lty=2)
abline(v=32.5,lty=2)
text(x=1:36+1.5,y=par("usr")[3]-0.5,labels=df_names1,xpd=NA,srt=25,adj=1.1,cex=0.9)
text(x=1:36+1.5,y=6.6,labels=n_names,xpd=NA,srt=0,adj=1,cex=0.75)
text(x=1,y=8.5,labels="(C)",xpd=NA)
legend("bottomleft",legend=c("T_R", "NT_R", "T_NR","NT_NR"),pch=c(15,15,15,15),col=c("#80cdc1","#018571","#dfc27d","#a6611a"),bty="n",cex=0.85)
dev.off()

#######Total stock change:

#hadgem26:
st2<-(read_data_annual("soilc_agr.bin",outp=outpath_h1,syr=2099,eyr=2099,ryr=2006)+
        read_data_annual("litc_agr.bin",outp=outpath_h1,syr=2099,eyr=2099,ryr=2006))*
  (area_cropland2005_32*is_tropical_wet)/1e+15 #convert to Petagramm (Pg)
st1<-(read_data_annual("soilc_agr.bin",outp=outpath_h1,syr=2006,eyr=2006,ryr=2006)+
        read_data_annual("litc_agr.bin",outp=outpath_h1,syr=2006,eyr=2006,ryr=2006))*
  (area_cropland2005_32*is_tropical_wet)/1e+15
s_trop_wet_h1<-st2-st1
rm(st1,st2)

st2<-(read_data_annual("soilc_agr.bin",outp=outpath_h1,syr=2099,eyr=2099,ryr=2006)+
        read_data_annual("litc_agr.bin",outp=outpath_h1,syr=2099,eyr=2099,ryr=2006))*
  (area_cropland2005_32*is_tropical_moist)/1e+15 #convert to Petagramm (Pg)
st1<-(read_data_annual("soilc_agr.bin",outp=outpath_h1,syr=2006,eyr=2006,ryr=2006)+
        read_data_annual("litc_agr.bin",outp=outpath_h1,syr=2006,eyr=2006,ryr=2006))*
  (area_cropland2005_32*is_tropical_moist)/1e+15
s_trop_moist_h1<-st2-st1
rm(st1,st2)

st2<-(read_data_annual("soilc_agr.bin",outp=outpath_h1,syr=2099,eyr=2099,ryr=2006)+
        read_data_annual("litc_agr.bin",outp=outpath_h1,syr=2099,eyr=2099,ryr=2006))*
  (area_cropland2005_32*is_tropical_dry)/1e+15 #convert to Petagramm (Pg)
st1<-(read_data_annual("soilc_agr.bin",outp=outpath_h1,syr=2006,eyr=2006,ryr=2006)+
        read_data_annual("litc_agr.bin",outp=outpath_h1,syr=2006,eyr=2006,ryr=2006))*
  (area_cropland2005_32*is_tropical_dry)/1e+15
s_trop_dry_h1<-st2-st1
rm(st1,st2)

st2<-(read_data_annual("soilc_agr.bin",outp=outpath_h1,syr=2099,eyr=2099,ryr=2006)+
        read_data_annual("litc_agr.bin",outp=outpath_h1,syr=2099,eyr=2099,ryr=2006))*
  (area_cropland2005_32*is_warm_temperate_moist)/1e+15 #convert to Petagramm (Pg)
st1<-(read_data_annual("soilc_agr.bin",outp=outpath_h1,syr=2006,eyr=2006,ryr=2006)+
        read_data_annual("litc_agr.bin",outp=outpath_h1,syr=2006,eyr=2006,ryr=2006))*
  (area_cropland2005_32*is_warm_temperate_moist)/1e+15
s_warmt_moist_h1<-st2-st1
rm(st1,st2)

st2<-(read_data_annual("soilc_agr.bin",outp=outpath_h1,syr=2099,eyr=2099,ryr=2006)+
        read_data_annual("litc_agr.bin",outp=outpath_h1,syr=2099,eyr=2099,ryr=2006))*
  (area_cropland2005_32*is_warm_temperate_dry)/1e+15 #convert to Petagramm (Pg)
st1<-(read_data_annual("soilc_agr.bin",outp=outpath_h1,syr=2006,eyr=2006,ryr=2006)+
        read_data_annual("litc_agr.bin",outp=outpath_h1,syr=2006,eyr=2006,ryr=2006))*
  (area_cropland2005_32*is_warm_temperate_dry)/1e+15
s_warmt_dry_h1<-st2-st1
rm(st1,st2)

st2<-(read_data_annual("soilc_agr.bin",outp=outpath_h1,syr=2099,eyr=2099,ryr=2006)+
        read_data_annual("litc_agr.bin",outp=outpath_h1,syr=2099,eyr=2099,ryr=2006))*
  (area_cropland2005_32*is_cold_temperate_moist)/1e+15 #convert to Petagramm (Pg)
st1<-(read_data_annual("soilc_agr.bin",outp=outpath_h1,syr=2006,eyr=2006,ryr=2006)+
        read_data_annual("litc_agr.bin",outp=outpath_h1,syr=2006,eyr=2006,ryr=2006))*
  (area_cropland2005_32*is_cold_temperate_moist)/1e+15
s_coldt_moist_h1<-st2-st1
rm(st1,st2)

st2<-(read_data_annual("soilc_agr.bin",outp=outpath_h1,syr=2099,eyr=2099,ryr=2006)+
        read_data_annual("litc_agr.bin",outp=outpath_h1,syr=2099,eyr=2099,ryr=2006))*
  (area_cropland2005_32*is_cold_temperate_dry)/1e+15 #convert to Petagramm (Pg)
st1<-(read_data_annual("soilc_agr.bin",outp=outpath_h1,syr=2006,eyr=2006,ryr=2006)+
        read_data_annual("litc_agr.bin",outp=outpath_h1,syr=2006,eyr=2006,ryr=2006))*
  (area_cropland2005_32*is_cold_temperate_dry)/1e+15
s_coldt_dry_h1<-st2-st1
rm(st1,st2)

st2<-(read_data_annual("soilc_agr.bin",outp=outpath_h1,syr=2099,eyr=2099,ryr=2006)+
        read_data_annual("litc_agr.bin",outp=outpath_h1,syr=2099,eyr=2099,ryr=2006))*
  (area_cropland2005_32*is_boreal_moist)/1e+15 #convert to Petagramm (Pg)
st1<-(read_data_annual("soilc_agr.bin",outp=outpath_h1,syr=2006,eyr=2006,ryr=2006)+
        read_data_annual("litc_agr.bin",outp=outpath_h1,syr=2006,eyr=2006,ryr=2006))*
  (area_cropland2005_32*is_boreal_moist)/1e+15
s_boreal_moist_h1<-st2-st1
rm(st1,st2)

st2<-(read_data_annual("soilc_agr.bin",outp=outpath_h1,syr=2099,eyr=2099,ryr=2006)+
        read_data_annual("litc_agr.bin",outp=outpath_h1,syr=2099,eyr=2099,ryr=2006))*
  (area_cropland2005_32*is_boreal_dry)/1e+15 #convert to Petagramm (Pg)
st1<-(read_data_annual("soilc_agr.bin",outp=outpath_h1,syr=2006,eyr=2006,ryr=2006)+
        read_data_annual("litc_agr.bin",outp=outpath_h1,syr=2006,eyr=2006,ryr=2006))*
  (area_cropland2005_32*is_boreal_dry)/1e+15
s_boreal_dry_h1<-st2-st1
rm(st1,st2)

st2<-(read_data_annual("soilc_agr.bin",outp=outpath_h1,syr=2099,eyr=2099,ryr=2006)+
        read_data_annual("litc_agr.bin",outp=outpath_h1,syr=2099,eyr=2099,ryr=2006))*
  (area_cropland2005_32)/1e+15 #convert to Petagramm (Pg)
st1<-(read_data_annual("soilc_agr.bin",outp=outpath_h1,syr=2006,eyr=2006,ryr=2006)+
        read_data_annual("litc_agr.bin",outp=outpath_h1,syr=2006,eyr=2006,ryr=2006))*
  (area_cropland2005_32)/1e+15
s_allc_h1<-st2-st1
rm(st1,st2)

df_s_h1<-c(sum(s_trop_wet_h1[,1]),sum(s_trop_wet_h1[,2]),sum(s_trop_wet_h1[,3]),sum(s_trop_wet_h1[,4]),
           sum(s_trop_moist_h1[,1]),sum(s_trop_moist_h1[,2]),sum(s_trop_moist_h1[,3]),sum(s_trop_moist_h1[,4]),
           sum(s_trop_dry_h1[,1]),sum(s_trop_dry_h1[,2]),sum(s_trop_dry_h1[,3]),sum(s_trop_dry_h1[,4]),
           sum(s_warmt_moist_h1[,1]),sum(s_warmt_moist_h1[,2]),sum(s_warmt_moist_h1[,3]),sum(s_warmt_moist_h1[,4]),
           sum(s_warmt_dry_h1[,1]),sum(s_warmt_dry_h1[,2]),sum(s_warmt_dry_h1[,3]),sum(s_warmt_dry_h1[,4]),
           sum(s_coldt_moist_h1[,1]),sum(s_coldt_moist_h1[,2]),sum(s_coldt_moist_h1[,3]),sum(s_coldt_moist_h1[,4]),
           sum(s_coldt_dry_h1[,1]),sum(s_coldt_dry_h1[,2]),sum(s_coldt_dry_h1[,3]),sum(s_coldt_dry_h1[,4]),
           sum(s_boreal_moist_h1[,1]),sum(s_boreal_moist_h1[,2]),sum(s_boreal_moist_h1[,3]),sum(s_boreal_moist_h1[,4]),
           sum(s_boreal_dry_h1[,1]),sum(s_boreal_dry_h1[,2]),sum(s_boreal_dry_h1[,3]),sum(s_boreal_dry_h1[,4]))

#gfdl26:
st2<-(read_data_annual("soilc_agr.bin",outp=outpath_g1,syr=2099,eyr=2099,ryr=2006)+
        read_data_annual("litc_agr.bin",outp=outpath_g1,syr=2099,eyr=2099,ryr=2006))*
  (area_cropland2005_32*is_tropical_wet)/1e+15 #convert to Petagramm (Pg)
st1<-(read_data_annual("soilc_agr.bin",outp=outpath_g1,syr=2006,eyr=2006,ryr=2006)+
        read_data_annual("litc_agr.bin",outp=outpath_g1,syr=2006,eyr=2006,ryr=2006))*
  (area_cropland2005_32*is_tropical_wet)/1e+15
s_trop_wet_g1<-st2-st1
rm(st1,st2)

st2<-(read_data_annual("soilc_agr.bin",outp=outpath_g1,syr=2099,eyr=2099,ryr=2006)+
        read_data_annual("litc_agr.bin",outp=outpath_g1,syr=2099,eyr=2099,ryr=2006))*
  (area_cropland2005_32*is_tropical_moist)/1e+15 #convert to Petagramm (Pg)
st1<-(read_data_annual("soilc_agr.bin",outp=outpath_g1,syr=2006,eyr=2006,ryr=2006)+
        read_data_annual("litc_agr.bin",outp=outpath_g1,syr=2006,eyr=2006,ryr=2006))*
  (area_cropland2005_32*is_tropical_moist)/1e+15
s_trop_moist_g1<-st2-st1
rm(st1,st2)

st2<-(read_data_annual("soilc_agr.bin",outp=outpath_g1,syr=2099,eyr=2099,ryr=2006)+
        read_data_annual("litc_agr.bin",outp=outpath_g1,syr=2099,eyr=2099,ryr=2006))*
  (area_cropland2005_32*is_tropical_dry)/1e+15 #convert to Petagramm (Pg)
st1<-(read_data_annual("soilc_agr.bin",outp=outpath_g1,syr=2006,eyr=2006,ryr=2006)+
        read_data_annual("litc_agr.bin",outp=outpath_g1,syr=2006,eyr=2006,ryr=2006))*
  (area_cropland2005_32*is_tropical_dry)/1e+15
s_trop_dry_g1<-st2-st1
rm(st1,st2)

st2<-(read_data_annual("soilc_agr.bin",outp=outpath_g1,syr=2099,eyr=2099,ryr=2006)+
        read_data_annual("litc_agr.bin",outp=outpath_g1,syr=2099,eyr=2099,ryr=2006))*
  (area_cropland2005_32*is_warm_temperate_moist)/1e+15 #convert to Petagramm (Pg)
st1<-(read_data_annual("soilc_agr.bin",outp=outpath_g1,syr=2006,eyr=2006,ryr=2006)+
        read_data_annual("litc_agr.bin",outp=outpath_g1,syr=2006,eyr=2006,ryr=2006))*
  (area_cropland2005_32*is_warm_temperate_moist)/1e+15
s_warmt_moist_g1<-st2-st1
rm(st1,st2)

st2<-(read_data_annual("soilc_agr.bin",outp=outpath_g1,syr=2099,eyr=2099,ryr=2006)+
        read_data_annual("litc_agr.bin",outp=outpath_g1,syr=2099,eyr=2099,ryr=2006))*
  (area_cropland2005_32*is_warm_temperate_dry)/1e+15 #convert to Petagramm (Pg)
st1<-(read_data_annual("soilc_agr.bin",outp=outpath_g1,syr=2006,eyr=2006,ryr=2006)+
        read_data_annual("litc_agr.bin",outp=outpath_g1,syr=2006,eyr=2006,ryr=2006))*
  (area_cropland2005_32*is_warm_temperate_dry)/1e+15
s_warmt_dry_g1<-st2-st1
rm(st1,st2)

st2<-(read_data_annual("soilc_agr.bin",outp=outpath_g1,syr=2099,eyr=2099,ryr=2006)+
        read_data_annual("litc_agr.bin",outp=outpath_g1,syr=2099,eyr=2099,ryr=2006))*
  (area_cropland2005_32*is_cold_temperate_moist)/1e+15 #convert to Petagramm (Pg)
st1<-(read_data_annual("soilc_agr.bin",outp=outpath_g1,syr=2006,eyr=2006,ryr=2006)+
        read_data_annual("litc_agr.bin",outp=outpath_g1,syr=2006,eyr=2006,ryr=2006))*
  (area_cropland2005_32*is_cold_temperate_moist)/1e+15
s_coldt_moist_g1<-st2-st1
rm(st1,st2)

st2<-(read_data_annual("soilc_agr.bin",outp=outpath_g1,syr=2099,eyr=2099,ryr=2006)+
        read_data_annual("litc_agr.bin",outp=outpath_g1,syr=2099,eyr=2099,ryr=2006))*
  (area_cropland2005_32*is_cold_temperate_dry)/1e+15 #convert to Petagramm (Pg)
st1<-(read_data_annual("soilc_agr.bin",outp=outpath_g1,syr=2006,eyr=2006,ryr=2006)+
        read_data_annual("litc_agr.bin",outp=outpath_g1,syr=2006,eyr=2006,ryr=2006))*
  (area_cropland2005_32*is_cold_temperate_dry)/1e+15
s_coldt_dry_g1<-st2-st1
rm(st1,st2)

st2<-(read_data_annual("soilc_agr.bin",outp=outpath_g1,syr=2099,eyr=2099,ryr=2006)+
        read_data_annual("litc_agr.bin",outp=outpath_g1,syr=2099,eyr=2099,ryr=2006))*
  (area_cropland2005_32*is_boreal_moist)/1e+15 #convert to Petagramm (Pg)
st1<-(read_data_annual("soilc_agr.bin",outp=outpath_g1,syr=2006,eyr=2006,ryr=2006)+
        read_data_annual("litc_agr.bin",outp=outpath_g1,syr=2006,eyr=2006,ryr=2006))*
  (area_cropland2005_32*is_boreal_moist)/1e+15
s_boreal_moist_g1<-st2-st1
rm(st1,st2)

st2<-(read_data_annual("soilc_agr.bin",outp=outpath_g1,syr=2099,eyr=2099,ryr=2006)+
        read_data_annual("litc_agr.bin",outp=outpath_g1,syr=2099,eyr=2099,ryr=2006))*
  (area_cropland2005_32*is_boreal_dry)/1e+15 #convert to Petagramm (Pg)
st1<-(read_data_annual("soilc_agr.bin",outp=outpath_g1,syr=2006,eyr=2006,ryr=2006)+
        read_data_annual("litc_agr.bin",outp=outpath_g1,syr=2006,eyr=2006,ryr=2006))*
  (area_cropland2005_32*is_boreal_dry)/1e+15
s_boreal_dry_g1<-st2-st1
rm(st1,st2)

st2<-(read_data_annual("soilc_agr.bin",outp=outpath_g1,syr=2099,eyr=2099,ryr=2006)+
        read_data_annual("litc_agr.bin",outp=outpath_g1,syr=2099,eyr=2099,ryr=2006))*
  (area_cropland2005_32)/1e+15 #convert to Petagramm (Pg)
st1<-(read_data_annual("soilc_agr.bin",outp=outpath_g1,syr=2006,eyr=2006,ryr=2006)+
        read_data_annual("litc_agr.bin",outp=outpath_g1,syr=2006,eyr=2006,ryr=2006))*
  (area_cropland2005_32)/1e+15
s_allc_g1<-st2-st1
rm(st1,st2)

df_s_g1<-c(sum(s_trop_wet_g1[,1]),sum(s_trop_wet_g1[,2]),sum(s_trop_wet_g1[,3]),sum(s_trop_wet_g1[,4]),
           sum(s_trop_moist_g1[,1]),sum(s_trop_moist_g1[,2]),sum(s_trop_moist_g1[,3]),sum(s_trop_moist_g1[,4]),
           sum(s_trop_dry_g1[,1]),sum(s_trop_dry_g1[,2]),sum(s_trop_dry_g1[,3]),sum(s_trop_dry_g1[,4]),
           sum(s_warmt_moist_g1[,1]),sum(s_warmt_moist_g1[,2]),sum(s_warmt_moist_g1[,3]),sum(s_warmt_moist_g1[,4]),
           sum(s_warmt_dry_g1[,1]),sum(s_warmt_dry_g1[,2]),sum(s_warmt_dry_g1[,3]),sum(s_warmt_dry_g1[,4]),
           sum(s_coldt_moist_g1[,1]),sum(s_coldt_moist_g1[,2]),sum(s_coldt_moist_g1[,3]),sum(s_coldt_moist_g1[,4]),
           sum(s_coldt_dry_g1[,1]),sum(s_coldt_dry_g1[,2]),sum(s_coldt_dry_g1[,3]),sum(s_coldt_dry_g1[,4]),
           sum(s_boreal_moist_g1[,1]),sum(s_boreal_moist_g1[,2]),sum(s_boreal_moist_g1[,3]),sum(s_boreal_moist_g1[,4]),
           sum(s_boreal_dry_g1[,1]),sum(s_boreal_dry_g1[,2]),sum(s_boreal_dry_g1[,3]),sum(s_boreal_dry_g1[,4]))

#ipsl26:
st2<-(read_data_annual("soilc_agr.bin",outp=outpath_i1,syr=2099,eyr=2099,ryr=2006)+
        read_data_annual("litc_agr.bin",outp=outpath_i1,syr=2099,eyr=2099,ryr=2006))*
  (area_cropland2005_32*is_tropical_wet)/1e+15 #convert to Petagramm (Pg)
st1<-(read_data_annual("soilc_agr.bin",outp=outpath_i1,syr=2006,eyr=2006,ryr=2006)+
        read_data_annual("litc_agr.bin",outp=outpath_i1,syr=2006,eyr=2006,ryr=2006))*
  (area_cropland2005_32*is_tropical_wet)/1e+15
s_trop_wet_i1<-st2-st1
rm(st1,st2)

st2<-(read_data_annual("soilc_agr.bin",outp=outpath_i1,syr=2099,eyr=2099,ryr=2006)+
        read_data_annual("litc_agr.bin",outp=outpath_i1,syr=2099,eyr=2099,ryr=2006))*
  (area_cropland2005_32*is_tropical_moist)/1e+15 #convert to Petagramm (Pg)
st1<-(read_data_annual("soilc_agr.bin",outp=outpath_i1,syr=2006,eyr=2006,ryr=2006)+
        read_data_annual("litc_agr.bin",outp=outpath_i1,syr=2006,eyr=2006,ryr=2006))*
  (area_cropland2005_32*is_tropical_moist)/1e+15
s_trop_moist_i1<-st2-st1
rm(st1,st2)

st2<-(read_data_annual("soilc_agr.bin",outp=outpath_i1,syr=2099,eyr=2099,ryr=2006)+
        read_data_annual("litc_agr.bin",outp=outpath_i1,syr=2099,eyr=2099,ryr=2006))*
  (area_cropland2005_32*is_tropical_dry)/1e+15 #convert to Petagramm (Pg)
st1<-(read_data_annual("soilc_agr.bin",outp=outpath_i1,syr=2006,eyr=2006,ryr=2006)+
        read_data_annual("litc_agr.bin",outp=outpath_i1,syr=2006,eyr=2006,ryr=2006))*
  (area_cropland2005_32*is_tropical_dry)/1e+15
s_trop_dry_i1<-st2-st1
rm(st1,st2)

st2<-(read_data_annual("soilc_agr.bin",outp=outpath_i1,syr=2099,eyr=2099,ryr=2006)+
        read_data_annual("litc_agr.bin",outp=outpath_i1,syr=2099,eyr=2099,ryr=2006))*
  (area_cropland2005_32*is_warm_temperate_moist)/1e+15 #convert to Petagramm (Pg)
st1<-(read_data_annual("soilc_agr.bin",outp=outpath_i1,syr=2006,eyr=2006,ryr=2006)+
        read_data_annual("litc_agr.bin",outp=outpath_i1,syr=2006,eyr=2006,ryr=2006))*
  (area_cropland2005_32*is_warm_temperate_moist)/1e+15
s_warmt_moist_i1<-st2-st1
rm(st1,st2)

st2<-(read_data_annual("soilc_agr.bin",outp=outpath_i1,syr=2099,eyr=2099,ryr=2006)+
        read_data_annual("litc_agr.bin",outp=outpath_i1,syr=2099,eyr=2099,ryr=2006))*
  (area_cropland2005_32*is_warm_temperate_dry)/1e+15 #convert to Petagramm (Pg)
st1<-(read_data_annual("soilc_agr.bin",outp=outpath_i1,syr=2006,eyr=2006,ryr=2006)+
        read_data_annual("litc_agr.bin",outp=outpath_i1,syr=2006,eyr=2006,ryr=2006))*
  (area_cropland2005_32*is_warm_temperate_dry)/1e+15
s_warmt_dry_i1<-st2-st1
rm(st1,st2)

st2<-(read_data_annual("soilc_agr.bin",outp=outpath_i1,syr=2099,eyr=2099,ryr=2006)+
        read_data_annual("litc_agr.bin",outp=outpath_i1,syr=2099,eyr=2099,ryr=2006))*
  (area_cropland2005_32*is_cold_temperate_moist)/1e+15 #convert to Petagramm (Pg)
st1<-(read_data_annual("soilc_agr.bin",outp=outpath_i1,syr=2006,eyr=2006,ryr=2006)+
        read_data_annual("litc_agr.bin",outp=outpath_i1,syr=2006,eyr=2006,ryr=2006))*
  (area_cropland2005_32*is_cold_temperate_moist)/1e+15
s_coldt_moist_i1<-st2-st1
rm(st1,st2)

st2<-(read_data_annual("soilc_agr.bin",outp=outpath_i1,syr=2099,eyr=2099,ryr=2006)+
        read_data_annual("litc_agr.bin",outp=outpath_i1,syr=2099,eyr=2099,ryr=2006))*
  (area_cropland2005_32*is_cold_temperate_dry)/1e+15 #convert to Petagramm (Pg)
st1<-(read_data_annual("soilc_agr.bin",outp=outpath_i1,syr=2006,eyr=2006,ryr=2006)+
        read_data_annual("litc_agr.bin",outp=outpath_i1,syr=2006,eyr=2006,ryr=2006))*
  (area_cropland2005_32*is_cold_temperate_dry)/1e+15
s_coldt_dry_i1<-st2-st1
rm(st1,st2)

st2<-(read_data_annual("soilc_agr.bin",outp=outpath_i1,syr=2099,eyr=2099,ryr=2006)+
        read_data_annual("litc_agr.bin",outp=outpath_i1,syr=2099,eyr=2099,ryr=2006))*
  (area_cropland2005_32*is_boreal_moist)/1e+15 #convert to Petagramm (Pg)
st1<-(read_data_annual("soilc_agr.bin",outp=outpath_i1,syr=2006,eyr=2006,ryr=2006)+
        read_data_annual("litc_agr.bin",outp=outpath_i1,syr=2006,eyr=2006,ryr=2006))*
  (area_cropland2005_32*is_boreal_moist)/1e+15
s_boreal_moist_i1<-st2-st1
rm(st1,st2)

st2<-(read_data_annual("soilc_agr.bin",outp=outpath_i1,syr=2099,eyr=2099,ryr=2006)+
        read_data_annual("litc_agr.bin",outp=outpath_i1,syr=2099,eyr=2099,ryr=2006))*
  (area_cropland2005_32*is_boreal_dry)/1e+15 #convert to Petagramm (Pg)
st1<-(read_data_annual("soilc_agr.bin",outp=outpath_i1,syr=2006,eyr=2006,ryr=2006)+
        read_data_annual("litc_agr.bin",outp=outpath_i1,syr=2006,eyr=2006,ryr=2006))*
  (area_cropland2005_32*is_boreal_dry)/1e+15
s_boreal_dry_i1<-st2-st1
rm(st1,st2)

st2<-(read_data_annual("soilc_agr.bin",outp=outpath_i1,syr=2099,eyr=2099,ryr=2006)+
        read_data_annual("litc_agr.bin",outp=outpath_i1,syr=2099,eyr=2099,ryr=2006))*
  (area_cropland2005_32)/1e+15 #convert to Petagramm (Pg)
st1<-(read_data_annual("soilc_agr.bin",outp=outpath_i1,syr=2006,eyr=2006,ryr=2006)+
        read_data_annual("litc_agr.bin",outp=outpath_i1,syr=2006,eyr=2006,ryr=2006))*
  (area_cropland2005_32)/1e+15
s_allc_i1<-st2-st1
rm(st1,st2)

df_s_i1<-c(sum(s_trop_wet_i1[,1]),sum(s_trop_wet_i1[,2]),sum(s_trop_wet_i1[,3]),sum(s_trop_wet_i1[,4]),
           sum(s_trop_moist_i1[,1]),sum(s_trop_moist_i1[,2]),sum(s_trop_moist_i1[,3]),sum(s_trop_moist_i1[,4]),
           sum(s_trop_dry_i1[,1]),sum(s_trop_dry_i1[,2]),sum(s_trop_dry_i1[,3]),sum(s_trop_dry_i1[,4]),
           sum(s_warmt_moist_i1[,1]),sum(s_warmt_moist_i1[,2]),sum(s_warmt_moist_i1[,3]),sum(s_warmt_moist_i1[,4]),
           sum(s_warmt_dry_i1[,1]),sum(s_warmt_dry_i1[,2]),sum(s_warmt_dry_i1[,3]),sum(s_warmt_dry_i1[,4]),
           sum(s_coldt_moist_i1[,1]),sum(s_coldt_moist_i1[,2]),sum(s_coldt_moist_i1[,3]),sum(s_coldt_moist_i1[,4]),
           sum(s_coldt_dry_i1[,1]),sum(s_coldt_dry_i1[,2]),sum(s_coldt_dry_i1[,3]),sum(s_coldt_dry_i1[,4]),
           sum(s_boreal_moist_i1[,1]),sum(s_boreal_moist_i1[,2]),sum(s_boreal_moist_i1[,3]),sum(s_boreal_moist_i1[,4]),
           sum(s_boreal_dry_i1[,1]),sum(s_boreal_dry_i1[,2]),sum(s_boreal_dry_i1[,3]),sum(s_boreal_dry_i1[,4]))

#miroc26:
st2<-(read_data_annual("soilc_agr.bin",outp=outpath_m1,syr=2099,eyr=2099,ryr=2006)+
        read_data_annual("litc_agr.bin",outp=outpath_m1,syr=2099,eyr=2099,ryr=2006))*
  (area_cropland2005_32*is_tropical_wet)/1e+15 #convert to Petagramm (Pg)
st1<-(read_data_annual("soilc_agr.bin",outp=outpath_m1,syr=2006,eyr=2006,ryr=2006)+
        read_data_annual("litc_agr.bin",outp=outpath_m1,syr=2006,eyr=2006,ryr=2006))*
  (area_cropland2005_32*is_tropical_wet)/1e+15
s_trop_wet_m1<-st2-st1
rm(st1,st2)

st2<-(read_data_annual("soilc_agr.bin",outp=outpath_m1,syr=2099,eyr=2099,ryr=2006)+
        read_data_annual("litc_agr.bin",outp=outpath_m1,syr=2099,eyr=2099,ryr=2006))*
  (area_cropland2005_32*is_tropical_moist)/1e+15 #convert to Petagramm (Pg)
st1<-(read_data_annual("soilc_agr.bin",outp=outpath_m1,syr=2006,eyr=2006,ryr=2006)+
        read_data_annual("litc_agr.bin",outp=outpath_m1,syr=2006,eyr=2006,ryr=2006))*
  (area_cropland2005_32*is_tropical_moist)/1e+15
s_trop_moist_m1<-st2-st1
rm(st1,st2)

st2<-(read_data_annual("soilc_agr.bin",outp=outpath_m1,syr=2099,eyr=2099,ryr=2006)+
        read_data_annual("litc_agr.bin",outp=outpath_m1,syr=2099,eyr=2099,ryr=2006))*
  (area_cropland2005_32*is_tropical_dry)/1e+15 #convert to Petagramm (Pg)
st1<-(read_data_annual("soilc_agr.bin",outp=outpath_m1,syr=2006,eyr=2006,ryr=2006)+
        read_data_annual("litc_agr.bin",outp=outpath_m1,syr=2006,eyr=2006,ryr=2006))*
  (area_cropland2005_32*is_tropical_dry)/1e+15
s_trop_dry_m1<-st2-st1
rm(st1,st2)

st2<-(read_data_annual("soilc_agr.bin",outp=outpath_m1,syr=2099,eyr=2099,ryr=2006)+
        read_data_annual("litc_agr.bin",outp=outpath_m1,syr=2099,eyr=2099,ryr=2006))*
  (area_cropland2005_32*is_warm_temperate_moist)/1e+15 #convert to Petagramm (Pg)
st1<-(read_data_annual("soilc_agr.bin",outp=outpath_m1,syr=2006,eyr=2006,ryr=2006)+
        read_data_annual("litc_agr.bin",outp=outpath_m1,syr=2006,eyr=2006,ryr=2006))*
  (area_cropland2005_32*is_warm_temperate_moist)/1e+15
s_warmt_moist_m1<-st2-st1
rm(st1,st2)

st2<-(read_data_annual("soilc_agr.bin",outp=outpath_m1,syr=2099,eyr=2099,ryr=2006)+
        read_data_annual("litc_agr.bin",outp=outpath_m1,syr=2099,eyr=2099,ryr=2006))*
  (area_cropland2005_32*is_warm_temperate_dry)/1e+15 #convert to Petagramm (Pg)
st1<-(read_data_annual("soilc_agr.bin",outp=outpath_m1,syr=2006,eyr=2006,ryr=2006)+
        read_data_annual("litc_agr.bin",outp=outpath_m1,syr=2006,eyr=2006,ryr=2006))*
  (area_cropland2005_32*is_warm_temperate_dry)/1e+15
s_warmt_dry_m1<-st2-st1
rm(st1,st2)

st2<-(read_data_annual("soilc_agr.bin",outp=outpath_m1,syr=2099,eyr=2099,ryr=2006)+
        read_data_annual("litc_agr.bin",outp=outpath_m1,syr=2099,eyr=2099,ryr=2006))*
  (area_cropland2005_32*is_cold_temperate_moist)/1e+15 #convert to Petagramm (Pg)
st1<-(read_data_annual("soilc_agr.bin",outp=outpath_m1,syr=2006,eyr=2006,ryr=2006)+
        read_data_annual("litc_agr.bin",outp=outpath_m1,syr=2006,eyr=2006,ryr=2006))*
  (area_cropland2005_32*is_cold_temperate_moist)/1e+15
s_coldt_moist_m1<-st2-st1
rm(st1,st2)

st2<-(read_data_annual("soilc_agr.bin",outp=outpath_m1,syr=2099,eyr=2099,ryr=2006)+
        read_data_annual("litc_agr.bin",outp=outpath_m1,syr=2099,eyr=2099,ryr=2006))*
  (area_cropland2005_32*is_cold_temperate_dry)/1e+15 #convert to Petagramm (Pg)
st1<-(read_data_annual("soilc_agr.bin",outp=outpath_m1,syr=2006,eyr=2006,ryr=2006)+
        read_data_annual("litc_agr.bin",outp=outpath_m1,syr=2006,eyr=2006,ryr=2006))*
  (area_cropland2005_32*is_cold_temperate_dry)/1e+15
s_coldt_dry_m1<-st2-st1
rm(st1,st2)

st2<-(read_data_annual("soilc_agr.bin",outp=outpath_m1,syr=2099,eyr=2099,ryr=2006)+
        read_data_annual("litc_agr.bin",outp=outpath_m1,syr=2099,eyr=2099,ryr=2006))*
  (area_cropland2005_32*is_boreal_moist)/1e+15 #convert to Petagramm (Pg)
st1<-(read_data_annual("soilc_agr.bin",outp=outpath_m1,syr=2006,eyr=2006,ryr=2006)+
        read_data_annual("litc_agr.bin",outp=outpath_m1,syr=2006,eyr=2006,ryr=2006))*
  (area_cropland2005_32*is_boreal_moist)/1e+15
s_boreal_moist_m1<-st2-st1
rm(st1,st2)

st2<-(read_data_annual("soilc_agr.bin",outp=outpath_m1,syr=2099,eyr=2099,ryr=2006)+
        read_data_annual("litc_agr.bin",outp=outpath_m1,syr=2099,eyr=2099,ryr=2006))*
  (area_cropland2005_32*is_boreal_dry)/1e+15 #convert to Petagramm (Pg)
st1<-(read_data_annual("soilc_agr.bin",outp=outpath_m1,syr=2006,eyr=2006,ryr=2006)+
        read_data_annual("litc_agr.bin",outp=outpath_m1,syr=2006,eyr=2006,ryr=2006))*
  (area_cropland2005_32*is_boreal_dry)/1e+15
s_boreal_dry_m1<-st2-st1
rm(st1,st2)

st2<-(read_data_annual("soilc_agr.bin",outp=outpath_m1,syr=2099,eyr=2099,ryr=2006)+
        read_data_annual("litc_agr.bin",outp=outpath_m1,syr=2099,eyr=2099,ryr=2006))*
  (area_cropland2005_32)/1e+15 #convert to Petagramm (Pg)
st1<-(read_data_annual("soilc_agr.bin",outp=outpath_m1,syr=2006,eyr=2006,ryr=2006)+
        read_data_annual("litc_agr.bin",outp=outpath_m1,syr=2006,eyr=2006,ryr=2006))*
  (area_cropland2005_32)/1e+15
s_allc_m1<-st2-st1
rm(st1,st2)

df_s_m1<-c(sum(s_trop_wet_m1[,1]),sum(s_trop_wet_m1[,2]),sum(s_trop_wet_m1[,3]),sum(s_trop_wet_m1[,4]),
           sum(s_trop_moist_m1[,1]),sum(s_trop_moist_m1[,2]),sum(s_trop_moist_m1[,3]),sum(s_trop_moist_m1[,4]),
           sum(s_trop_dry_m1[,1]),sum(s_trop_dry_m1[,2]),sum(s_trop_dry_m1[,3]),sum(s_trop_dry_m1[,4]),
           sum(s_warmt_moist_m1[,1]),sum(s_warmt_moist_m1[,2]),sum(s_warmt_moist_m1[,3]),sum(s_warmt_moist_m1[,4]),
           sum(s_warmt_dry_m1[,1]),sum(s_warmt_dry_m1[,2]),sum(s_warmt_dry_m1[,3]),sum(s_warmt_dry_m1[,4]),
           sum(s_coldt_moist_m1[,1]),sum(s_coldt_moist_m1[,2]),sum(s_coldt_moist_m1[,3]),sum(s_coldt_moist_m1[,4]),
           sum(s_coldt_dry_m1[,1]),sum(s_coldt_dry_m1[,2]),sum(s_coldt_dry_m1[,3]),sum(s_coldt_dry_m1[,4]),
           sum(s_boreal_moist_m1[,1]),sum(s_boreal_moist_m1[,2]),sum(s_boreal_moist_m1[,3]),sum(s_boreal_moist_m1[,4]),
           sum(s_boreal_dry_m1[,1]),sum(s_boreal_dry_m1[,2]),sum(s_boreal_dry_m1[,3]),sum(s_boreal_dry_m1[,4]))

df_s_average_26<-(df_s_h1+df_s_g1+df_s_i1+df_s_m1)/4

df_names2<-c("Tropical wet","","","","Tropical moist","","","","Tropical dry","",
             "","","","Warm temperate moist","","","","Warm temperate dry","","",
             "","","Cold temperate moist","","","","","Cold temperate dry",
             "","","","Boreal moist","","","","Boreal dry")

pdf(paste0(plots_path,"barplot_total_SOC_change_regions_averages_rcp26_B.pdf"))
barplot(df_s_average_26,ylim=c(-13,3), col=myCol, yaxt="n",ylab="Cropland SOC change 2006-2099 (Pg C)",
        cex.lab=0.8,#main="RCP2.6",cex.main=0.8,
        # width=0.8,xlim=c(0,36),
        space=c(0,rep(0,3),0.5,(rep(0,3)),0.5,(rep(0,3)),0.5,
                (rep(0,3)),0.5,(rep(0,3)),0.5,(rep(0,3)),0.5,
                (rep(0,3)),0.5,(rep(0,3)),0.5,(rep(0,3))))
axis(2, at=seq(-12, 2, 2),las=2)
axis(side=1,at=c(4.25,8.75,13.25,17.75,22.25,26.75,31.25,35.75),labels=FALSE)
abline(0,0,lty=1)
abline(v=4.25,lty=2)
abline(v=8.75,lty=2)
abline(v=13.25,lty=2)
abline(v=17.75,lty=2)
abline(v=22.25,lty=2)
abline(v=26.75,lty=2)
abline(v=31.25,lty=2)
abline(v=35.75,lty=2)
text(x=1:36+1.5,y=par("usr")[3]-0.5,labels=df_names2,xpd=NA,srt=25,adj=0.9,cex=0.9)
text(x=1,y=3.5,labels="(B)",xpd=NA)
legend("bottomleft",legend=c("T_R", "NT_R", "T_NR","NT_NR"),pch=c(15,15,15,15),col=c("#80cdc1","#018571","#dfc27d","#a6611a"),bty="n",cex=0.85)
box()
dev.off()

#Total stocks RCP8.5:

#hadgem26:
st2<-(read_data_annual("soilc_agr.bin",outp=outpath_h2,syr=2099,eyr=2099,ryr=2006)+
        read_data_annual("litc_agr.bin",outp=outpath_h2,syr=2099,eyr=2099,ryr=2006))*
  (area_cropland2005_32*is_tropical_wet)/1e+15 #convert to Petagramm (Pg)
st1<-(read_data_annual("soilc_agr.bin",outp=outpath_h2,syr=2006,eyr=2006,ryr=2006)+
        read_data_annual("litc_agr.bin",outp=outpath_h2,syr=2006,eyr=2006,ryr=2006))*
  (area_cropland2005_32*is_tropical_wet)/1e+15
s_trop_wet_h2<-st2-st1
rm(st1,st2)

st2<-(read_data_annual("soilc_agr.bin",outp=outpath_h2,syr=2099,eyr=2099,ryr=2006)+
        read_data_annual("litc_agr.bin",outp=outpath_h2,syr=2099,eyr=2099,ryr=2006))*
  (area_cropland2005_32*is_tropical_moist)/1e+15 #convert to Petagramm (Pg)
st1<-(read_data_annual("soilc_agr.bin",outp=outpath_h2,syr=2006,eyr=2006,ryr=2006)+
        read_data_annual("litc_agr.bin",outp=outpath_h2,syr=2006,eyr=2006,ryr=2006))*
  (area_cropland2005_32*is_tropical_moist)/1e+15
s_trop_moist_h2<-st2-st1
rm(st1,st2)

st2<-(read_data_annual("soilc_agr.bin",outp=outpath_h2,syr=2099,eyr=2099,ryr=2006)+
        read_data_annual("litc_agr.bin",outp=outpath_h2,syr=2099,eyr=2099,ryr=2006))*
  (area_cropland2005_32*is_tropical_dry)/1e+15 #convert to Petagramm (Pg)
st1<-(read_data_annual("soilc_agr.bin",outp=outpath_h2,syr=2006,eyr=2006,ryr=2006)+
        read_data_annual("litc_agr.bin",outp=outpath_h2,syr=2006,eyr=2006,ryr=2006))*
  (area_cropland2005_32*is_tropical_dry)/1e+15
s_trop_dry_h2<-st2-st1
rm(st1,st2)

st2<-(read_data_annual("soilc_agr.bin",outp=outpath_h2,syr=2099,eyr=2099,ryr=2006)+
        read_data_annual("litc_agr.bin",outp=outpath_h2,syr=2099,eyr=2099,ryr=2006))*
  (area_cropland2005_32*is_warm_temperate_moist)/1e+15 #convert to Petagramm (Pg)
st1<-(read_data_annual("soilc_agr.bin",outp=outpath_h2,syr=2006,eyr=2006,ryr=2006)+
        read_data_annual("litc_agr.bin",outp=outpath_h2,syr=2006,eyr=2006,ryr=2006))*
  (area_cropland2005_32*is_warm_temperate_moist)/1e+15
s_warmt_moist_h2<-st2-st1
rm(st1,st2)

st2<-(read_data_annual("soilc_agr.bin",outp=outpath_h2,syr=2099,eyr=2099,ryr=2006)+
        read_data_annual("litc_agr.bin",outp=outpath_h2,syr=2099,eyr=2099,ryr=2006))*
  (area_cropland2005_32*is_warm_temperate_dry)/1e+15 #convert to Petagramm (Pg)
st1<-(read_data_annual("soilc_agr.bin",outp=outpath_h2,syr=2006,eyr=2006,ryr=2006)+
        read_data_annual("litc_agr.bin",outp=outpath_h2,syr=2006,eyr=2006,ryr=2006))*
  (area_cropland2005_32*is_warm_temperate_dry)/1e+15
s_warmt_dry_h2<-st2-st1
rm(st1,st2)

st2<-(read_data_annual("soilc_agr.bin",outp=outpath_h2,syr=2099,eyr=2099,ryr=2006)+
        read_data_annual("litc_agr.bin",outp=outpath_h2,syr=2099,eyr=2099,ryr=2006))*
  (area_cropland2005_32*is_cold_temperate_moist)/1e+15 #convert to Petagramm (Pg)
st1<-(read_data_annual("soilc_agr.bin",outp=outpath_h2,syr=2006,eyr=2006,ryr=2006)+
        read_data_annual("litc_agr.bin",outp=outpath_h2,syr=2006,eyr=2006,ryr=2006))*
  (area_cropland2005_32*is_cold_temperate_moist)/1e+15
s_coldt_moist_h2<-st2-st1
rm(st1,st2)

st2<-(read_data_annual("soilc_agr.bin",outp=outpath_h2,syr=2099,eyr=2099,ryr=2006)+
        read_data_annual("litc_agr.bin",outp=outpath_h2,syr=2099,eyr=2099,ryr=2006))*
  (area_cropland2005_32*is_cold_temperate_dry)/1e+15 #convert to Petagramm (Pg)
st1<-(read_data_annual("soilc_agr.bin",outp=outpath_h2,syr=2006,eyr=2006,ryr=2006)+
        read_data_annual("litc_agr.bin",outp=outpath_h2,syr=2006,eyr=2006,ryr=2006))*
  (area_cropland2005_32*is_cold_temperate_dry)/1e+15
s_coldt_dry_h2<-st2-st1
rm(st1,st2)

st2<-(read_data_annual("soilc_agr.bin",outp=outpath_h2,syr=2099,eyr=2099,ryr=2006)+
        read_data_annual("litc_agr.bin",outp=outpath_h2,syr=2099,eyr=2099,ryr=2006))*
  (area_cropland2005_32*is_boreal_moist)/1e+15 #convert to Petagramm (Pg)
st1<-(read_data_annual("soilc_agr.bin",outp=outpath_h2,syr=2006,eyr=2006,ryr=2006)+
        read_data_annual("litc_agr.bin",outp=outpath_h2,syr=2006,eyr=2006,ryr=2006))*
  (area_cropland2005_32*is_boreal_moist)/1e+15
s_boreal_moist_h2<-st2-st1
rm(st1,st2)

st2<-(read_data_annual("soilc_agr.bin",outp=outpath_h2,syr=2099,eyr=2099,ryr=2006)+
        read_data_annual("litc_agr.bin",outp=outpath_h2,syr=2099,eyr=2099,ryr=2006))*
  (area_cropland2005_32*is_boreal_dry)/1e+15 #convert to Petagramm (Pg)
st1<-(read_data_annual("soilc_agr.bin",outp=outpath_h2,syr=2006,eyr=2006,ryr=2006)+
        read_data_annual("litc_agr.bin",outp=outpath_h2,syr=2006,eyr=2006,ryr=2006))*
  (area_cropland2005_32*is_boreal_dry)/1e+15
s_boreal_dry_h2<-st2-st1
rm(st1,st2)

st2<-(read_data_annual("soilc_agr.bin",outp=outpath_h2,syr=2099,eyr=2099,ryr=2006)+
        read_data_annual("litc_agr.bin",outp=outpath_h2,syr=2099,eyr=2099,ryr=2006))*
  (area_cropland2005_32)/1e+15 #convert to Petagramm (Pg)
st1<-(read_data_annual("soilc_agr.bin",outp=outpath_h2,syr=2006,eyr=2006,ryr=2006)+
        read_data_annual("litc_agr.bin",outp=outpath_h2,syr=2006,eyr=2006,ryr=2006))*
  (area_cropland2005_32)/1e+15
s_allc_h2<-st2-st1
rm(st1,st2)

df_s_h2<-c(sum(s_trop_wet_h2[,1]),sum(s_trop_wet_h2[,2]),sum(s_trop_wet_h2[,3]),sum(s_trop_wet_h2[,4]),
           sum(s_trop_moist_h2[,1]),sum(s_trop_moist_h2[,2]),sum(s_trop_moist_h2[,3]),sum(s_trop_moist_h2[,4]),
           sum(s_trop_dry_h2[,1]),sum(s_trop_dry_h2[,2]),sum(s_trop_dry_h2[,3]),sum(s_trop_dry_h2[,4]),
           sum(s_warmt_moist_h2[,1]),sum(s_warmt_moist_h2[,2]),sum(s_warmt_moist_h2[,3]),sum(s_warmt_moist_h2[,4]),
           sum(s_warmt_dry_h2[,1]),sum(s_warmt_dry_h2[,2]),sum(s_warmt_dry_h2[,3]),sum(s_warmt_dry_h2[,4]),
           sum(s_coldt_moist_h2[,1]),sum(s_coldt_moist_h2[,2]),sum(s_coldt_moist_h2[,3]),sum(s_coldt_moist_h2[,4]),
           sum(s_coldt_dry_h2[,1]),sum(s_coldt_dry_h2[,2]),sum(s_coldt_dry_h2[,3]),sum(s_coldt_dry_h2[,4]),
           sum(s_boreal_moist_h2[,1]),sum(s_boreal_moist_h2[,2]),sum(s_boreal_moist_h2[,3]),sum(s_boreal_moist_h2[,4]),
           sum(s_boreal_dry_h2[,1]),sum(s_boreal_dry_h2[,2]),sum(s_boreal_dry_h2[,3]),sum(s_boreal_dry_h2[,4]))

#gfdl85:
st2<-(read_data_annual("soilc_agr.bin",outp=outpath_g2,syr=2099,eyr=2099,ryr=2006)+
        read_data_annual("litc_agr.bin",outp=outpath_g2,syr=2099,eyr=2099,ryr=2006))*
  (area_cropland2005_32*is_tropical_wet)/1e+15 #convert to Petagramm (Pg)
st1<-(read_data_annual("soilc_agr.bin",outp=outpath_g2,syr=2006,eyr=2006,ryr=2006)+
        read_data_annual("litc_agr.bin",outp=outpath_g2,syr=2006,eyr=2006,ryr=2006))*
  (area_cropland2005_32*is_tropical_wet)/1e+15
s_trop_wet_g2<-st2-st1
rm(st1,st2)

st2<-(read_data_annual("soilc_agr.bin",outp=outpath_g2,syr=2099,eyr=2099,ryr=2006)+
        read_data_annual("litc_agr.bin",outp=outpath_g2,syr=2099,eyr=2099,ryr=2006))*
  (area_cropland2005_32*is_tropical_moist)/1e+15 #convert to Petagramm (Pg)
st1<-(read_data_annual("soilc_agr.bin",outp=outpath_g2,syr=2006,eyr=2006,ryr=2006)+
        read_data_annual("litc_agr.bin",outp=outpath_g2,syr=2006,eyr=2006,ryr=2006))*
  (area_cropland2005_32*is_tropical_moist)/1e+15
s_trop_moist_g2<-st2-st1
rm(st1,st2)

st2<-(read_data_annual("soilc_agr.bin",outp=outpath_g2,syr=2099,eyr=2099,ryr=2006)+
        read_data_annual("litc_agr.bin",outp=outpath_g2,syr=2099,eyr=2099,ryr=2006))*
  (area_cropland2005_32*is_tropical_dry)/1e+15 #convert to Petagramm (Pg)
st1<-(read_data_annual("soilc_agr.bin",outp=outpath_g2,syr=2006,eyr=2006,ryr=2006)+
        read_data_annual("litc_agr.bin",outp=outpath_g2,syr=2006,eyr=2006,ryr=2006))*
  (area_cropland2005_32*is_tropical_dry)/1e+15
s_trop_dry_g2<-st2-st1
rm(st1,st2)

st2<-(read_data_annual("soilc_agr.bin",outp=outpath_g2,syr=2099,eyr=2099,ryr=2006)+
        read_data_annual("litc_agr.bin",outp=outpath_g2,syr=2099,eyr=2099,ryr=2006))*
  (area_cropland2005_32*is_warm_temperate_moist)/1e+15 #convert to Petagramm (Pg)
st1<-(read_data_annual("soilc_agr.bin",outp=outpath_g2,syr=2006,eyr=2006,ryr=2006)+
        read_data_annual("litc_agr.bin",outp=outpath_g2,syr=2006,eyr=2006,ryr=2006))*
  (area_cropland2005_32*is_warm_temperate_moist)/1e+15
s_warmt_moist_g2<-st2-st1
rm(st1,st2)

st2<-(read_data_annual("soilc_agr.bin",outp=outpath_g2,syr=2099,eyr=2099,ryr=2006)+
        read_data_annual("litc_agr.bin",outp=outpath_g2,syr=2099,eyr=2099,ryr=2006))*
  (area_cropland2005_32*is_warm_temperate_dry)/1e+15 #convert to Petagramm (Pg)
st1<-(read_data_annual("soilc_agr.bin",outp=outpath_g2,syr=2006,eyr=2006,ryr=2006)+
        read_data_annual("litc_agr.bin",outp=outpath_g2,syr=2006,eyr=2006,ryr=2006))*
  (area_cropland2005_32*is_warm_temperate_dry)/1e+15
s_warmt_dry_g2<-st2-st1
rm(st1,st2)

st2<-(read_data_annual("soilc_agr.bin",outp=outpath_g2,syr=2099,eyr=2099,ryr=2006)+
        read_data_annual("litc_agr.bin",outp=outpath_g2,syr=2099,eyr=2099,ryr=2006))*
  (area_cropland2005_32*is_cold_temperate_moist)/1e+15 #convert to Petagramm (Pg)
st1<-(read_data_annual("soilc_agr.bin",outp=outpath_g2,syr=2006,eyr=2006,ryr=2006)+
        read_data_annual("litc_agr.bin",outp=outpath_g2,syr=2006,eyr=2006,ryr=2006))*
  (area_cropland2005_32*is_cold_temperate_moist)/1e+15
s_coldt_moist_g2<-st2-st1
rm(st1,st2)

st2<-(read_data_annual("soilc_agr.bin",outp=outpath_g2,syr=2099,eyr=2099,ryr=2006)+
        read_data_annual("litc_agr.bin",outp=outpath_g2,syr=2099,eyr=2099,ryr=2006))*
  (area_cropland2005_32*is_cold_temperate_dry)/1e+15 #convert to Petagramm (Pg)
st1<-(read_data_annual("soilc_agr.bin",outp=outpath_g2,syr=2006,eyr=2006,ryr=2006)+
        read_data_annual("litc_agr.bin",outp=outpath_g2,syr=2006,eyr=2006,ryr=2006))*
  (area_cropland2005_32*is_cold_temperate_dry)/1e+15
s_coldt_dry_g2<-st2-st1
rm(st1,st2)

st2<-(read_data_annual("soilc_agr.bin",outp=outpath_g2,syr=2099,eyr=2099,ryr=2006)+
        read_data_annual("litc_agr.bin",outp=outpath_g2,syr=2099,eyr=2099,ryr=2006))*
  (area_cropland2005_32*is_boreal_moist)/1e+15 #convert to Petagramm (Pg)
st1<-(read_data_annual("soilc_agr.bin",outp=outpath_g2,syr=2006,eyr=2006,ryr=2006)+
        read_data_annual("litc_agr.bin",outp=outpath_g2,syr=2006,eyr=2006,ryr=2006))*
  (area_cropland2005_32*is_boreal_moist)/1e+15
s_boreal_moist_g2<-st2-st1
rm(st1,st2)

st2<-(read_data_annual("soilc_agr.bin",outp=outpath_g2,syr=2099,eyr=2099,ryr=2006)+
        read_data_annual("litc_agr.bin",outp=outpath_g2,syr=2099,eyr=2099,ryr=2006))*
  (area_cropland2005_32*is_boreal_dry)/1e+15 #convert to Petagramm (Pg)
st1<-(read_data_annual("soilc_agr.bin",outp=outpath_g2,syr=2006,eyr=2006,ryr=2006)+
        read_data_annual("litc_agr.bin",outp=outpath_g2,syr=2006,eyr=2006,ryr=2006))*
  (area_cropland2005_32*is_boreal_dry)/1e+15
s_boreal_dry_g2<-st2-st1
rm(st1,st2)

st2<-(read_data_annual("soilc_agr.bin",outp=outpath_g2,syr=2099,eyr=2099,ryr=2006)+
        read_data_annual("litc_agr.bin",outp=outpath_g2,syr=2099,eyr=2099,ryr=2006))*
  (area_cropland2005_32)/1e+15 #convert to Petagramm (Pg)
st1<-(read_data_annual("soilc_agr.bin",outp=outpath_g2,syr=2006,eyr=2006,ryr=2006)+
        read_data_annual("litc_agr.bin",outp=outpath_g2,syr=2006,eyr=2006,ryr=2006))*
  (area_cropland2005_32)/1e+15
s_allc_g2<-st2-st1
rm(st1,st2)

df_s_g2<-c(sum(s_trop_wet_g2[,1]),sum(s_trop_wet_g2[,2]),sum(s_trop_wet_g2[,3]),sum(s_trop_wet_g2[,4]),
           sum(s_trop_moist_g2[,1]),sum(s_trop_moist_g2[,2]),sum(s_trop_moist_g2[,3]),sum(s_trop_moist_g2[,4]),
           sum(s_trop_dry_g2[,1]),sum(s_trop_dry_g2[,2]),sum(s_trop_dry_g2[,3]),sum(s_trop_dry_g2[,4]),
           sum(s_warmt_moist_g2[,1]),sum(s_warmt_moist_g2[,2]),sum(s_warmt_moist_g2[,3]),sum(s_warmt_moist_g2[,4]),
           sum(s_warmt_dry_g2[,1]),sum(s_warmt_dry_g2[,2]),sum(s_warmt_dry_g2[,3]),sum(s_warmt_dry_g2[,4]),
           sum(s_coldt_moist_g2[,1]),sum(s_coldt_moist_g2[,2]),sum(s_coldt_moist_g2[,3]),sum(s_coldt_moist_g2[,4]),
           sum(s_coldt_dry_g2[,1]),sum(s_coldt_dry_g2[,2]),sum(s_coldt_dry_g2[,3]),sum(s_coldt_dry_g2[,4]),
           sum(s_boreal_moist_g2[,1]),sum(s_boreal_moist_g2[,2]),sum(s_boreal_moist_g2[,3]),sum(s_boreal_moist_g2[,4]),
           sum(s_boreal_dry_g2[,1]),sum(s_boreal_dry_g2[,2]),sum(s_boreal_dry_g2[,3]),sum(s_boreal_dry_g2[,4]))

#ipsl85:
st2<-(read_data_annual("soilc_agr.bin",outp=outpath_i2,syr=2099,eyr=2099,ryr=2006)+
        read_data_annual("litc_agr.bin",outp=outpath_i2,syr=2099,eyr=2099,ryr=2006))*
  (area_cropland2005_32*is_tropical_wet)/1e+15 #convert to Petagramm (Pg)
st1<-(read_data_annual("soilc_agr.bin",outp=outpath_i2,syr=2006,eyr=2006,ryr=2006)+
        read_data_annual("litc_agr.bin",outp=outpath_i2,syr=2006,eyr=2006,ryr=2006))*
  (area_cropland2005_32*is_tropical_wet)/1e+15
s_trop_wet_i2<-st2-st1
rm(st1,st2)

st2<-(read_data_annual("soilc_agr.bin",outp=outpath_i2,syr=2099,eyr=2099,ryr=2006)+
        read_data_annual("litc_agr.bin",outp=outpath_i2,syr=2099,eyr=2099,ryr=2006))*
  (area_cropland2005_32*is_tropical_moist)/1e+15 #convert to Petagramm (Pg)
st1<-(read_data_annual("soilc_agr.bin",outp=outpath_i2,syr=2006,eyr=2006,ryr=2006)+
        read_data_annual("litc_agr.bin",outp=outpath_i2,syr=2006,eyr=2006,ryr=2006))*
  (area_cropland2005_32*is_tropical_moist)/1e+15
s_trop_moist_i2<-st2-st1
rm(st1,st2)

st2<-(read_data_annual("soilc_agr.bin",outp=outpath_i2,syr=2099,eyr=2099,ryr=2006)+
        read_data_annual("litc_agr.bin",outp=outpath_i2,syr=2099,eyr=2099,ryr=2006))*
  (area_cropland2005_32*is_tropical_dry)/1e+15 #convert to Petagramm (Pg)
st1<-(read_data_annual("soilc_agr.bin",outp=outpath_i2,syr=2006,eyr=2006,ryr=2006)+
        read_data_annual("litc_agr.bin",outp=outpath_i2,syr=2006,eyr=2006,ryr=2006))*
  (area_cropland2005_32*is_tropical_dry)/1e+15
s_trop_dry_i2<-st2-st1
rm(st1,st2)

st2<-(read_data_annual("soilc_agr.bin",outp=outpath_i2,syr=2099,eyr=2099,ryr=2006)+
        read_data_annual("litc_agr.bin",outp=outpath_i2,syr=2099,eyr=2099,ryr=2006))*
  (area_cropland2005_32*is_warm_temperate_moist)/1e+15 #convert to Petagramm (Pg)
st1<-(read_data_annual("soilc_agr.bin",outp=outpath_i2,syr=2006,eyr=2006,ryr=2006)+
        read_data_annual("litc_agr.bin",outp=outpath_i2,syr=2006,eyr=2006,ryr=2006))*
  (area_cropland2005_32*is_warm_temperate_moist)/1e+15
s_warmt_moist_i2<-st2-st1
rm(st1,st2)

st2<-(read_data_annual("soilc_agr.bin",outp=outpath_i2,syr=2099,eyr=2099,ryr=2006)+
        read_data_annual("litc_agr.bin",outp=outpath_i2,syr=2099,eyr=2099,ryr=2006))*
  (area_cropland2005_32*is_warm_temperate_dry)/1e+15 #convert to Petagramm (Pg)
st1<-(read_data_annual("soilc_agr.bin",outp=outpath_i2,syr=2006,eyr=2006,ryr=2006)+
        read_data_annual("litc_agr.bin",outp=outpath_i2,syr=2006,eyr=2006,ryr=2006))*
  (area_cropland2005_32*is_warm_temperate_dry)/1e+15
s_warmt_dry_i2<-st2-st1
rm(st1,st2)

st2<-(read_data_annual("soilc_agr.bin",outp=outpath_i2,syr=2099,eyr=2099,ryr=2006)+
        read_data_annual("litc_agr.bin",outp=outpath_i2,syr=2099,eyr=2099,ryr=2006))*
  (area_cropland2005_32*is_cold_temperate_moist)/1e+15 #convert to Petagramm (Pg)
st1<-(read_data_annual("soilc_agr.bin",outp=outpath_i2,syr=2006,eyr=2006,ryr=2006)+
        read_data_annual("litc_agr.bin",outp=outpath_i2,syr=2006,eyr=2006,ryr=2006))*
  (area_cropland2005_32*is_cold_temperate_moist)/1e+15
s_coldt_moist_i2<-st2-st1
rm(st1,st2)

st2<-(read_data_annual("soilc_agr.bin",outp=outpath_i2,syr=2099,eyr=2099,ryr=2006)+
        read_data_annual("litc_agr.bin",outp=outpath_i2,syr=2099,eyr=2099,ryr=2006))*
  (area_cropland2005_32*is_cold_temperate_dry)/1e+15 #convert to Petagramm (Pg)
st1<-(read_data_annual("soilc_agr.bin",outp=outpath_i2,syr=2006,eyr=2006,ryr=2006)+
        read_data_annual("litc_agr.bin",outp=outpath_i2,syr=2006,eyr=2006,ryr=2006))*
  (area_cropland2005_32*is_cold_temperate_dry)/1e+15
s_coldt_dry_i2<-st2-st1
rm(st1,st2)

st2<-(read_data_annual("soilc_agr.bin",outp=outpath_i2,syr=2099,eyr=2099,ryr=2006)+
        read_data_annual("litc_agr.bin",outp=outpath_i2,syr=2099,eyr=2099,ryr=2006))*
  (area_cropland2005_32*is_boreal_moist)/1e+15 #convert to Petagramm (Pg)
st1<-(read_data_annual("soilc_agr.bin",outp=outpath_i2,syr=2006,eyr=2006,ryr=2006)+
        read_data_annual("litc_agr.bin",outp=outpath_i2,syr=2006,eyr=2006,ryr=2006))*
  (area_cropland2005_32*is_boreal_moist)/1e+15
s_boreal_moist_i2<-st2-st1
rm(st1,st2)

st2<-(read_data_annual("soilc_agr.bin",outp=outpath_i2,syr=2099,eyr=2099,ryr=2006)+
        read_data_annual("litc_agr.bin",outp=outpath_i2,syr=2099,eyr=2099,ryr=2006))*
  (area_cropland2005_32*is_boreal_dry)/1e+15 #convert to Petagramm (Pg)
st1<-(read_data_annual("soilc_agr.bin",outp=outpath_i2,syr=2006,eyr=2006,ryr=2006)+
        read_data_annual("litc_agr.bin",outp=outpath_i2,syr=2006,eyr=2006,ryr=2006))*
  (area_cropland2005_32*is_boreal_dry)/1e+15
s_boreal_dry_i2<-st2-st1
rm(st1,st2)

st2<-(read_data_annual("soilc_agr.bin",outp=outpath_i2,syr=2099,eyr=2099,ryr=2006)+
        read_data_annual("litc_agr.bin",outp=outpath_i2,syr=2099,eyr=2099,ryr=2006))*
  (area_cropland2005_32)/1e+15 #convert to Petagramm (Pg)
st1<-(read_data_annual("soilc_agr.bin",outp=outpath_i2,syr=2006,eyr=2006,ryr=2006)+
        read_data_annual("litc_agr.bin",outp=outpath_i2,syr=2006,eyr=2006,ryr=2006))*
  (area_cropland2005_32)/1e+15
s_allc_i2<-st2-st1
rm(st1,st2)

df_s_i2<-c(sum(s_trop_wet_i2[,1]),sum(s_trop_wet_i2[,2]),sum(s_trop_wet_i2[,3]),sum(s_trop_wet_i2[,4]),
           sum(s_trop_moist_i2[,1]),sum(s_trop_moist_i2[,2]),sum(s_trop_moist_i2[,3]),sum(s_trop_moist_i2[,4]),
           sum(s_trop_dry_i2[,1]),sum(s_trop_dry_i2[,2]),sum(s_trop_dry_i2[,3]),sum(s_trop_dry_i2[,4]),
           sum(s_warmt_moist_i2[,1]),sum(s_warmt_moist_i2[,2]),sum(s_warmt_moist_i2[,3]),sum(s_warmt_moist_i2[,4]),
           sum(s_warmt_dry_i2[,1]),sum(s_warmt_dry_i2[,2]),sum(s_warmt_dry_i2[,3]),sum(s_warmt_dry_i2[,4]),
           sum(s_coldt_moist_i2[,1]),sum(s_coldt_moist_i2[,2]),sum(s_coldt_moist_i2[,3]),sum(s_coldt_moist_i2[,4]),
           sum(s_coldt_dry_i2[,1]),sum(s_coldt_dry_i2[,2]),sum(s_coldt_dry_i2[,3]),sum(s_coldt_dry_i2[,4]),
           sum(s_boreal_moist_i2[,1]),sum(s_boreal_moist_i2[,2]),sum(s_boreal_moist_i2[,3]),sum(s_boreal_moist_i2[,4]),
           sum(s_boreal_dry_i2[,1]),sum(s_boreal_dry_i2[,2]),sum(s_boreal_dry_i2[,3]),sum(s_boreal_dry_i2[,4]))

#miroc85:
st2<-(read_data_annual("soilc_agr.bin",outp=outpath_m2,syr=2099,eyr=2099,ryr=2006)+
        read_data_annual("litc_agr.bin",outp=outpath_m2,syr=2099,eyr=2099,ryr=2006))*
  (area_cropland2005_32*is_tropical_wet)/1e+15 #convert to Petagramm (Pg)
st1<-(read_data_annual("soilc_agr.bin",outp=outpath_m2,syr=2006,eyr=2006,ryr=2006)+
        read_data_annual("litc_agr.bin",outp=outpath_m2,syr=2006,eyr=2006,ryr=2006))*
  (area_cropland2005_32*is_tropical_wet)/1e+15
s_trop_wet_m2<-st2-st1
rm(st1,st2)

st2<-(read_data_annual("soilc_agr.bin",outp=outpath_m2,syr=2099,eyr=2099,ryr=2006)+
        read_data_annual("litc_agr.bin",outp=outpath_m2,syr=2099,eyr=2099,ryr=2006))*
  (area_cropland2005_32*is_tropical_moist)/1e+15 #convert to Petagramm (Pg)
st1<-(read_data_annual("soilc_agr.bin",outp=outpath_m2,syr=2006,eyr=2006,ryr=2006)+
        read_data_annual("litc_agr.bin",outp=outpath_m2,syr=2006,eyr=2006,ryr=2006))*
  (area_cropland2005_32*is_tropical_moist)/1e+15
s_trop_moist_m2<-st2-st1
rm(st1,st2)

st2<-(read_data_annual("soilc_agr.bin",outp=outpath_m2,syr=2099,eyr=2099,ryr=2006)+
        read_data_annual("litc_agr.bin",outp=outpath_m2,syr=2099,eyr=2099,ryr=2006))*
  (area_cropland2005_32*is_tropical_dry)/1e+15 #convert to Petagramm (Pg)
st1<-(read_data_annual("soilc_agr.bin",outp=outpath_m2,syr=2006,eyr=2006,ryr=2006)+
        read_data_annual("litc_agr.bin",outp=outpath_m2,syr=2006,eyr=2006,ryr=2006))*
  (area_cropland2005_32*is_tropical_dry)/1e+15
s_trop_dry_m2<-st2-st1
rm(st1,st2)

st2<-(read_data_annual("soilc_agr.bin",outp=outpath_m2,syr=2099,eyr=2099,ryr=2006)+
        read_data_annual("litc_agr.bin",outp=outpath_m2,syr=2099,eyr=2099,ryr=2006))*
  (area_cropland2005_32*is_warm_temperate_moist)/1e+15 #convert to Petagramm (Pg)
st1<-(read_data_annual("soilc_agr.bin",outp=outpath_m2,syr=2006,eyr=2006,ryr=2006)+
        read_data_annual("litc_agr.bin",outp=outpath_m2,syr=2006,eyr=2006,ryr=2006))*
  (area_cropland2005_32*is_warm_temperate_moist)/1e+15
s_warmt_moist_m2<-st2-st1
rm(st1,st2)

st2<-(read_data_annual("soilc_agr.bin",outp=outpath_m2,syr=2099,eyr=2099,ryr=2006)+
        read_data_annual("litc_agr.bin",outp=outpath_m2,syr=2099,eyr=2099,ryr=2006))*
  (area_cropland2005_32*is_warm_temperate_dry)/1e+15 #convert to Petagramm (Pg)
st1<-(read_data_annual("soilc_agr.bin",outp=outpath_m2,syr=2006,eyr=2006,ryr=2006)+
        read_data_annual("litc_agr.bin",outp=outpath_m2,syr=2006,eyr=2006,ryr=2006))*
  (area_cropland2005_32*is_warm_temperate_dry)/1e+15
s_warmt_dry_m2<-st2-st1
rm(st1,st2)

st2<-(read_data_annual("soilc_agr.bin",outp=outpath_m2,syr=2099,eyr=2099,ryr=2006)+
        read_data_annual("litc_agr.bin",outp=outpath_m2,syr=2099,eyr=2099,ryr=2006))*
  (area_cropland2005_32*is_cold_temperate_moist)/1e+15 #convert to Petagramm (Pg)
st1<-(read_data_annual("soilc_agr.bin",outp=outpath_m2,syr=2006,eyr=2006,ryr=2006)+
        read_data_annual("litc_agr.bin",outp=outpath_m2,syr=2006,eyr=2006,ryr=2006))*
  (area_cropland2005_32*is_cold_temperate_moist)/1e+15
s_coldt_moist_m2<-st2-st1
rm(st1,st2)

st2<-(read_data_annual("soilc_agr.bin",outp=outpath_m2,syr=2099,eyr=2099,ryr=2006)+
        read_data_annual("litc_agr.bin",outp=outpath_m2,syr=2099,eyr=2099,ryr=2006))*
  (area_cropland2005_32*is_cold_temperate_dry)/1e+15 #convert to Petagramm (Pg)
st1<-(read_data_annual("soilc_agr.bin",outp=outpath_m2,syr=2006,eyr=2006,ryr=2006)+
        read_data_annual("litc_agr.bin",outp=outpath_m2,syr=2006,eyr=2006,ryr=2006))*
  (area_cropland2005_32*is_cold_temperate_dry)/1e+15
s_coldt_dry_m2<-st2-st1
rm(st1,st2)

st2<-(read_data_annual("soilc_agr.bin",outp=outpath_m2,syr=2099,eyr=2099,ryr=2006)+
        read_data_annual("litc_agr.bin",outp=outpath_m2,syr=2099,eyr=2099,ryr=2006))*
  (area_cropland2005_32*is_boreal_moist)/1e+15 #convert to Petagramm (Pg)
st1<-(read_data_annual("soilc_agr.bin",outp=outpath_m2,syr=2006,eyr=2006,ryr=2006)+
        read_data_annual("litc_agr.bin",outp=outpath_m2,syr=2006,eyr=2006,ryr=2006))*
  (area_cropland2005_32*is_boreal_moist)/1e+15
s_boreal_moist_m2<-st2-st1
rm(st1,st2)

st2<-(read_data_annual("soilc_agr.bin",outp=outpath_m2,syr=2099,eyr=2099,ryr=2006)+
        read_data_annual("litc_agr.bin",outp=outpath_m2,syr=2099,eyr=2099,ryr=2006))*
  (area_cropland2005_32*is_boreal_dry)/1e+15 #convert to Petagramm (Pg)
st1<-(read_data_annual("soilc_agr.bin",outp=outpath_m2,syr=2006,eyr=2006,ryr=2006)+
        read_data_annual("litc_agr.bin",outp=outpath_m2,syr=2006,eyr=2006,ryr=2006))*
  (area_cropland2005_32*is_boreal_dry)/1e+15
s_boreal_dry_m2<-st2-st1
rm(st1,st2)

st2<-(read_data_annual("soilc_agr.bin",outp=outpath_m2,syr=2099,eyr=2099,ryr=2006)+
        read_data_annual("litc_agr.bin",outp=outpath_m2,syr=2099,eyr=2099,ryr=2006))*
  (area_cropland2005_32)/1e+15 #convert to Petagramm (Pg)
st1<-(read_data_annual("soilc_agr.bin",outp=outpath_m2,syr=2006,eyr=2006,ryr=2006)+
        read_data_annual("litc_agr.bin",outp=outpath_m2,syr=2006,eyr=2006,ryr=2006))*
  (area_cropland2005_32)/1e+15
s_allc_m2<-st2-st1
rm(st1,st2)

df_s_m2<-c(sum(s_trop_wet_m2[,1]),sum(s_trop_wet_m2[,2]),sum(s_trop_wet_m2[,3]),sum(s_trop_wet_m2[,4]),
           sum(s_trop_moist_m2[,1]),sum(s_trop_moist_m2[,2]),sum(s_trop_moist_m2[,3]),sum(s_trop_moist_m2[,4]),
           sum(s_trop_dry_m2[,1]),sum(s_trop_dry_m2[,2]),sum(s_trop_dry_m2[,3]),sum(s_trop_dry_m2[,4]),
           sum(s_warmt_moist_m2[,1]),sum(s_warmt_moist_m2[,2]),sum(s_warmt_moist_m2[,3]),sum(s_warmt_moist_m2[,4]),
           sum(s_warmt_dry_m2[,1]),sum(s_warmt_dry_m2[,2]),sum(s_warmt_dry_m2[,3]),sum(s_warmt_dry_m2[,4]),
           sum(s_coldt_moist_m2[,1]),sum(s_coldt_moist_m2[,2]),sum(s_coldt_moist_m2[,3]),sum(s_coldt_moist_m2[,4]),
           sum(s_coldt_dry_m2[,1]),sum(s_coldt_dry_m2[,2]),sum(s_coldt_dry_m2[,3]),sum(s_coldt_dry_m2[,4]),
           sum(s_boreal_moist_m2[,1]),sum(s_boreal_moist_m2[,2]),sum(s_boreal_moist_m2[,3]),sum(s_boreal_moist_m2[,4]),
           sum(s_boreal_dry_m2[,1]),sum(s_boreal_dry_m2[,2]),sum(s_boreal_dry_m2[,3]),sum(s_boreal_dry_m2[,4]))

df_s_average_85<-(df_s_h2+df_s_g2+df_s_i2+df_s_m2)/4

df_names2<-c("Tropical wet","","","","Tropical moist","","","","Tropical dry","",
             "","","","Warm temperate moist","","","","Warm temperate dry","","",
             "","","Cold temperate moist","","","","","Cold temperate dry",
             "","","","Boreal moist","","","","Boreal dry")

pdf(paste0(plots_path,"barplot_total_SOC_change_regions_averages_rcp85_D.pdf"))
barplot(df_s_average_85,ylim=c(-13,3), col=myCol, yaxt="n",ylab="Cropland SOC change 2006-2099 (Pg C)",
        cex.lab=0.8,#main="RCP8.5",cex.main=0.8,
        # width=0.8,xlim=c(0,36),
        space=c(0,rep(0,3),0.5,(rep(0,3)),0.5,(rep(0,3)),0.5,
                (rep(0,3)),0.5,(rep(0,3)),0.5,(rep(0,3)),0.5,
                (rep(0,3)),0.5,(rep(0,3)),0.5,(rep(0,3))))
axis(2, at=seq(-12, 2, 2),las=2)
axis(side=1,at=c(4.25,8.75,13.25,17.75,22.25,26.75,31.25,35.75),labels=FALSE)
abline(0,0,lty=1)
abline(v=4.25,lty=2)
abline(v=8.75,lty=2)
abline(v=13.25,lty=2)
abline(v=17.75,lty=2)
abline(v=22.25,lty=2)
abline(v=26.75,lty=2)
abline(v=31.25,lty=2)
abline(v=35.75,lty=2)
text(x=1:36+1.5,y=par("usr")[3]-0.5,labels=df_names2,xpd=NA,srt=25,adj=0.9,cex=0.9)
text(x=1,y=3.5,labels="(D)",xpd=NA)
legend("bottomleft",legend=c("T_R", "NT_R", "T_NR","NT_NR"),pch=c(15,15,15,15),col=c("#80cdc1","#018571","#dfc27d","#a6611a"),bty="n",cex=0.85)
box()
dev.off()
