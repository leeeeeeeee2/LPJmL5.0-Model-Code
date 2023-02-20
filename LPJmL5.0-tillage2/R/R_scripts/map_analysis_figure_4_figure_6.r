###############################################################################################
#                                                                                             #
#  Script by Tobias Herzfeld (04.2020) for:                                                   #
# "SOC sequestration potentials for agricultural management practices under climate change"   #
#                                                                                             #
#  Reads in the output data created by LPJmL5.0-tillage2 model runs and creates the SOC       #
#  difference maps figure 4 and 6.                                                            #
#                                                                                             #
###############################################################################################

rm(list=ls(all=T));gc()

library(oce)
library(Hmisc)
library(weights)
library(plotrix)

#define paths:
output_path<-"path/.../"                #path of model output data for rcps
plots_path<-"path/.../"                 #path for output plots created
sfpath<-"path/.../R/shapefiles/"        #path for shhapefiles
landuse_path<-"path/.../"               #path for landuse input data
grid_path<-"path/.../R/grid/"           #path for lpj grid

#define scens, NCELL and years:
scens <- c("T_NR","T_R","NT_NR","NT_R")
NCELL <- 67420
refyear<-2006
syear <- 2099
eyear <- 2099
nyears <- length(syear:eyear)

dirnames <- array(dim=length(scens),dimnames=list(scens))

dirnames["T_NR"] <- "OUTPUT_T_NR"
dirnames["T_R"] <- "OUTPUT_T_R"
dirnames["NT_NR"] <- "OUTPUT_NT_NR"
dirnames["NT_R"] <- "OUTPUT_NT_R"

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
#functions:

read_data_annual <- function(filename="",syear,eyear)
{
  nyears=length(syear:eyear)
  data.out <- array(data=0,dim=c(NCELL,length(scens)),dimnames=list(NULL,scens))
  for(scen in scens)
  {
    f <- file(paste0(outpath,dirnames[scen],"/",filename),"rb")
    seek(f,where=NCELL*(syear-refyear)*4,origin="start")
    for(iy in 1:nyears) data.out[,scen] <- data.out[,scen] + readBin(f,double(),n=NCELL,size=4)
    close(f)
  }
  data.out/nyears
}

read_data_annual_total <- function(filename="",syear,eyear,area=cellarea,scaling=1e+15)
{
  nyears=length(syear:eyear)
  data.out <- array(data=0,dim=c(NCELL,length(scens)),dimnames=list(NULL,scens))
  for(scen in scens)
  {
    f <- file(paste0(outpath,dirnames[scen],"/",filename),"rb")
    seek(f,where=NCELL*(syear-refyear)*4,origin="start")
    for(iy in 1:nyears) data.out[,scen] <- data.out[,scen] + readBin(f,double(),n=NCELL,size=4)
    data.out<-data.out*area/scaling
    close(f)
  }
  data.out/nyears
}

read_data_landuse <- function(filename="",syear=2000, refyear=1700, lu=64)
{
  data.out <- array(data=0,dim=c(NCELL,lu)) #
  f <- file(filename,"rb")
  seek(f,where=43+NCELL*(syear-refyear)*2*lu,origin="start")
  data.out <- matrix(readBin(f,integer(),n=NCELL*64,size=2), ncol=lu,byrow=T)
  colnames(data.out) <- bbnames_cft
  close(f)
  data.out/1000
}

##read grid + calcutate cellarea:

zz <- file(paste0(grid_path,"grid.bin"),"rb")
seek(zz,where=43,origin = "start")
x <- readBin(zz, integer(), n=2*NCELL, size=2) / 100
lon <- x[c(1:NCELL)*2-1]
lat <- x[c(1:NCELL)*2]
cellarea <- (111e3*0.5)*(111e3*0.5)*cos(lat/180*pi)
ilon <- as.integer((lon+180)/0.5 + 1.01)
ilat <- as.integer((lat+90)/0.5 + 1.01)
mat.ilonilat <- matrix(c(ilon,ilat),nrow=NCELL,ncol=2,byrow=F)
close(zz);rm(zz,x);gc()

##read landuse data

data_landuse <- read_data_landuse(paste0(landuse_path,"cft1700_2005_irrigation_systems_64bands.bin"))
frac_cropland <- rowSums(data_landuse[,-c(rep(13:16,4)+rep(16*(0:3),each=4))])
is_cropland <- ifelse(frac_cropland>0,1,NA)
area_cropland <- frac_cropland*cellarea # in m2

#############################################################################################################

library(maptools)
library(RColorBrewer)
library(fields)

world_wgs<-readShapePoly(paste0(sfpath,"world_borders/world_country_admin_boundary_shapefile_with_fips_codes.shp"),proj4string=CRS("+proj=longlat +datum=WGS84"))
load(paste0(sfpath,"world_borders/mask_land.RData"))

plot_maps <- function(data.plot,data.area,mapbreaks=c(-1,-0.4,-0.2,-0.1,-0.05,0.025,0.025,0.05,0.1,0.2,0.4,1e30),histbreaks=seq(-1,1,0.1),col=brewer.pal(11,"Spectral"),filename="test.png",main="")
{
  map.plot <- array(dim=c(720,360))
  
  bitmap(filename,type="png16m",height=14+4,width=30,res=200)
  layout(matrix(c(1,2,5,3,4,5),3,2),heights=c(0.5,0.5,4/14))
  par(mar=c(0,0,2.1,0))
  ylim <- c(-50,90)
  panels <- character(0)
  panels[c("T_NR","NT_NR","T_R","NT_R")] <- c("A","B","C","D")
  for(scen in names(panels))
  {
    map.plot[mat.ilonilat] <- ifelse(data.area>0,1,NA)*(data.plot[,scen])
    #map.plot[mat.ilonilat] <- vec_cropland(=is_cropland)
    plot(mask_land,border=NA,col="white",ylim=ylim,main=paste0(panels[scen],": ",scen))
    image(x=seq(-180,180,0.5),y=seq(-90,90,0.5),map.plot,ylim=ylim,xlab="",ylab="",xaxt="n",yaxt="n",breaks=mapbreaks,col=col,add=T)
    plot(mask_land,border=NA,col="white",ylim=ylim,add=T)
    plot(world_wgs,ylim=ylim,add=T)
    plotInset(xleft=-180,ybottom=-60,xright=-90,ytop=20,expr={wtd.hist(x=pmin(pmax(data.plot[,scen],min(histbreaks)),max(histbreaks)),weight=area_cropland*1e-10,breaks=histbreaks,xlab="",ylab="",col=grey(0.5),cex.axis=0.8,main="");mtext("Cropland in Mha",2,line=2,cex=0.6)})
  }
  
  plot.new()
  par(fig=pmax(par("fig"),0)) # this is needed since fig is set incorrectly by plot.new()
  breaks.legend <- seq(-1,1,len=length(mapbreaks))
  labs.legend <- c(paste0("<",mapbreaks[2]),mapbreaks[2:(length(mapbreaks)-1)],paste0(">",mapbreaks[length(mapbreaks)-1]))
  image.plot(zlim=c(-1,1),breaks=breaks.legend,col=col,legend.only=TRUE,horizontal=TRUE,legend.shrink=0.8,legend.width=3,legend.mar=12.1,axis.args=list(at=breaks.legend,labels=labs.legend))
  text(x=0.5,y=0.1,labels=main,cex=0.8,font=2)
  
  dev.off()
}

##read data and create maps for soc_agr,
##difference between 2099-2006 for each rcp26 and rcp85:

col.allc_change <- brewer.pal(11,"Spectral")
histbreaks.allc <- seq(0,20,1)
mapbreaks.allc_mapchange <- c(-20,-10,-5,-1,-0.5,-0.1,0.1,0.5,1,5,10,1e30)
histbreaks.allc_mapchange <- seq(-20,20,1)
col.allc <- brewer.pal(11,"Spectral")

#rcp2.6:
outpath<-paste0(output_path,"rcp26_hadgem/")
data.allc1_hadgem26 <- read_data_annual("soilc_agr.bin",syear=2006,eyear=2006)+read_data_annual("litc_agr.bin",syear=2006,eyear=2006)
data.allc2_hadgem26 <- read_data_annual("soilc_agr.bin",syear=2099,eyear=2099)+read_data_annual("litc_agr.bin",syear=2099,eyear=2099)
outpath<-paste0(output_path,"rcp26_gfdl/")
data.allc1_gfdl26 <- read_data_annual("soilc_agr.bin",syear=2006,eyear=2006)+read_data_annual("litc_agr.bin",syear=2006,eyear=2006)
data.allc2_gfdl26 <- read_data_annual("soilc_agr.bin",syear=2099,eyear=2099)+read_data_annual("litc_agr.bin",syear=2099,eyear=2099)
outpath<-paste0(output_path,"rcp26_ipsl/")
data.allc1_ipsl26 <- read_data_annual("soilc_agr.bin",syear=2006,eyear=2006)+read_data_annual("litc_agr.bin",syear=2006,eyear=2006)
data.allc2_ipsl26 <- read_data_annual("soilc_agr.bin",syear=2099,eyear=2099)+read_data_annual("litc_agr.bin",syear=2099,eyear=2099)
outpath<-paste0(output_path,"rcp26_miroc/")
data.allc1_miroc26 <- read_data_annual("soilc_agr.bin",syear=2006,eyear=2006)+read_data_annual("litc_agr.bin",syear=2006,eyear=2006)
data.allc2_miroc26 <- read_data_annual("soilc_agr.bin",syear=2099,eyear=2099)+read_data_annual("litc_agr.bin",syear=2099,eyear=2099)

data.allc_agr_hadgem26<-(data.allc2_hadgem26 - data.allc1_hadgem26)/1000
data.allc_agr_gfdl26<-(data.allc2_gfdl26 - data.allc1_gfdl26)/1000
data.allc_agr_ipsl26<-(data.allc2_ipsl26 - data.allc1_ipsl26)/1000
data.allc_agr_miroc26<-(data.allc2_miroc26 - data.allc1_miroc26)/1000

plot_maps(data.allc_agr_hadgem26,area_cropland,mapbreaks=mapbreaks.allc_mapchange,histbreaks=histbreaks.allc_mapchange,
          col=col.allc,filename=paste0(plots_path,"maps_allc_AGR_ensemble_RCP26_hadgem.png"),
          main="Simulated cropland SOC density change (kg/"~m^2~") between year 2006 and 2099 for RCP2.6 (HadGEM2_ES)")

plot_maps(data.allc_agr_gfdl26,area_cropland,mapbreaks=mapbreaks.allc_mapchange,histbreaks=histbreaks.allc_mapchange,
          col=col.allc,filename=paste0(plots_path,"maps_allc_AGR_ensemble_RCP26_gfdl.png"),
          main="Simulated cropland SOC density change (kg/"~m^2~") between year 2006 and 2099 for RCP2.6 (GFDL-ESM2M)")

plot_maps(data.allc_agr_ipsl26,area_cropland,mapbreaks=mapbreaks.allc_mapchange,histbreaks=histbreaks.allc_mapchange,
          col=col.allc,filename=paste0(plots_path,"maps_allc_AGR_ensemble_RCP26_ipsl.png"),
          main="Simulated cropland SOC density change (kg/"~m^2~") between year 2006 and 2099 for RCP2.6 (IPSL-CM5A-LR)")

plot_maps(data.allc_agr_miroc26,area_cropland,mapbreaks=mapbreaks.allc_mapchange,histbreaks=histbreaks.allc_mapchange,
          col=col.allc,filename=paste0(plots_path,"maps_allc_AGR_ensemble_RCP26_miroc.png"),
          main="Simulated cropland SOC density change (kg/"~m^2~") between year 2006 and 2099 for RCP2.6 (MIROC5)")


data.allc1_26<-(data.allc1_hadgem26+data.allc1_gfdl26+data.allc1_ipsl26+data.allc1_miroc26)/4
data.allc2_26<-(data.allc2_hadgem26+data.allc2_gfdl26+data.allc2_ipsl26+data.allc2_miroc26)/4
data.allc_agr26 <- data.allc2_26 - data.allc1_26

#rcp8.5:
outpath<-paste0(output_path,"rcp85_hadgem/")
data.allc1_hadgem85 <- read_data_annual("soilc_agr.bin",syear=2006,eyear=2006)+read_data_annual("litc_agr.bin",syear=2006,eyear=2006)
data.allc2_hadgem85 <- read_data_annual("soilc_agr.bin",syear=2099,eyear=2099)+read_data_annual("litc_agr.bin",syear=2099,eyear=2099)
outpath<-paste0(output_path,"rcp85_gfdl/")
data.allc1_gfdl85 <- read_data_annual("soilc_agr.bin",syear=2006,eyear=2006)+read_data_annual("litc_agr.bin",syear=2006,eyear=2006)
data.allc2_gfdl85 <- read_data_annual("soilc_agr.bin",syear=2099,eyear=2099)+read_data_annual("litc_agr.bin",syear=2099,eyear=2099)
outpath<-paste0(output_path,"rcp85_ipsl/")
data.allc1_ipsl85 <- read_data_annual("soilc_agr.bin",syear=2006,eyear=2006)+read_data_annual("litc_agr.bin",syear=2006,eyear=2006)
data.allc2_ipsl85 <- read_data_annual("soilc_agr.bin",syear=2099,eyear=2099)+read_data_annual("litc_agr.bin",syear=2099,eyear=2099)
outpath<-paste0(output_path,"rcp85_miroc/")
data.allc1_miroc85 <- read_data_annual("soilc_agr.bin",syear=2006,eyear=2006)+read_data_annual("litc_agr.bin",syear=2006,eyear=2006)
data.allc2_miroc85 <- read_data_annual("soilc_agr.bin",syear=2099,eyear=2099)+read_data_annual("litc_agr.bin",syear=2099,eyear=2099)

data.allc_agr_hadgem85<-(data.allc2_hadgem85 - data.allc1_hadgem85)/1000
data.allc_agr_gfdl85<-(data.allc2_gfdl85 - data.allc1_gfdl85)/1000
data.allc_agr_ipsl85<-(data.allc2_ipsl85 - data.allc1_ipsl85)/1000
data.allc_agr_miroc85<-(data.allc2_miroc85 - data.allc1_miroc85)/1000

plot_maps(data.allc_agr_hadgem85,area_cropland,mapbreaks=mapbreaks.allc_mapchange,histbreaks=histbreaks.allc_mapchange,
          col=col.allc,filename=paste0(plots_path,"maps_allc_AGR_ensemble_RCP85_hadgem.png"),
          main="Simulated cropland SOC density change (kg/"~m^2~") between year 2006 and 2099 for RCP8.5 (HadGEM2_ES)")

plot_maps(data.allc_agr_gfdl85,area_cropland,mapbreaks=mapbreaks.allc_mapchange,histbreaks=histbreaks.allc_mapchange,
          col=col.allc,filename=paste0(plots_path,"maps_allc_AGR_ensemble_RCP85_gfdl.png"),
          main="Simulated cropland SOC density change (kg/"~m^2~") between year 2006 and 2099 for RCP8.5 (GFDL-ESM2M)")

plot_maps(data.allc_agr_ipsl85,area_cropland,mapbreaks=mapbreaks.allc_mapchange,histbreaks=histbreaks.allc_mapchange,
          col=col.allc,filename=paste0(plots_path,"maps_allc_AGR_ensemble_RCP85_ipsl.png"),
          main="Simulated cropland SOC density change (kg/"~m^2~") between year 2006 and 2099 for RCP8.5 (IPSL-CM5A-LR)")

plot_maps(data.allc_agr_miroc85,area_cropland,mapbreaks=mapbreaks.allc_mapchange,histbreaks=histbreaks.allc_mapchange,
          col=col.allc,filename=paste0(plots_path,"maps_allc_AGR_ensemble_RCP85_miroc.png"),
          main="Simulated cropland SOC density change (kg/"~m^2~") between year 2006 and 2099 for RCP8.5 (MIROC5)")


data.allc1_85<-(data.allc1_hadgem85+data.allc1_gfdl85+data.allc1_ipsl85+data.allc1_miroc85)/4
data.allc2_85<-(data.allc2_hadgem85+data.allc2_gfdl85+data.allc2_ipsl85+data.allc2_miroc85)/4
data.allc_agr85 <- data.allc2_85 - data.allc1_85

#create change maps for same area affected:

data.affect<-ifelse(data.allc_agr26<0 & data.allc_agr85<0 , 1, 
                  ifelse(data.allc_agr26>0 & data.allc_agr85>0,2, 
                         ifelse(data.allc_agr26==0 & data.allc_agr85==0 , 0, 3)))

data.affect_hadgem<-ifelse(data.allc_agr_hadgem26<0 & data.allc_agr_hadgem85<0 , 1, 
                         ifelse(data.allc_agr_hadgem26>0 & data.allc_agr_hadgem85>0,2, 
                                ifelse(data.allc_agr_hadgem26==0 & data.allc_agr_hadgem85==0 , 0, 3)))

mapbreaks.allc_affect <- c(-0.5,0.5,1.5,2.5,1e30)
histbreaks.allc_affect <- seq(-0.5,3.5,1)
col.allc_change <- c("#ffffbf","#fdae61","#2c7bb6","#d7191c")

plot_maps(data.affect,area_cropland,mapbreaks=mapbreaks.allc_affect,histbreaks=histbreaks.allc_affect,
          col=col.allc_change,filename=paste0(plots_path,"maps_change_areas_allc_AGR_RCP26_85.png"),
          main="Maps of equal areas affected in SOC stocks change in RCP2.6 and RCP8.5 between 2006-2099\n(0=no change, 1=both decrease, 2=both increase, 3=opposite direction) )")

plot_maps(data.affect_hadgem,area_cropland,mapbreaks=mapbreaks.allc_affect,histbreaks=histbreaks.allc_affect,
          col=col.allc_change,filename=paste0(plots_path,"maps_change_areas_allc_AGR_RCP26_85_hadgem.png"),
          main="Maps of equal areas affected in SOC stocks change in RCP2.6 and RCP8.5 between 2006-2099 for HadGEM_ES\n(0=no change, 1=both decrease, 2=both increase, 3=opposite direction) )")

