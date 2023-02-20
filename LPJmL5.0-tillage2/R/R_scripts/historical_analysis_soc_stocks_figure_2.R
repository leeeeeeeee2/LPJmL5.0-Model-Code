###############################################################################################
#                                                                                             #
#  Script by Tobias Herzfeld (04.2020) for:                                                   #
# "SOC sequestration potentials for agricultural management practices under climate change"   #
#                                                                                             #
#  Reads in the data created by LPJmL5.0-tillage2 model runs and creates the plot cropland    #
#  outputs in figure 2. Calculates cropland SOC density for PNV and SOC losses from land-use  #
#  change in soc_cultivate.                                                                   #
#                                                                                             #
###############################################################################################

rm(list=ls(all=T));gc()

library(oce)
library(Hmisc)
library(weights)
library(plotrix)

##define paths:
output_path<-"path/.../"          #path of model output data
grid_path<-"path/.../R/grid/"     #path for lpj grid
landuse_path<-"path/.../"         #path for landuse input data
plots_path<-"path/.../"           #path for output plots created

#paths for simulations:
outp_pnv<-paste0(output_path,"OUTPUT_1700_PNV_1840_RESTART_7000/")
outp_lu<-paste0(output_path,"OUTPUT_1700_TILLAGE_1840_CRU403/")
outp_clu<-paste0(output_path,"OUTPUT_1700_TILLAGE_1840_CRU403_CONST_LU/")

##define NCELL and years:
NCELL <- 67420
refyear <- 1700
syear <- 1700
eyear <- 2018

nyears <- length(syear:eyear)
years<-eyear-syear+1

bnames_cft <- c("wheat","rice","maize","millet","pulses","sugarbeet","cassava","sunflower","soybeans","groundnuts","rapeseed","sugarcane","others","grasses","begrasses","betrees")
bbnames_cft<- c(paste0(bnames_cft,"_rf"),paste0(bnames_cft,"_surf"),paste0(bnames_cft,"_sprink"), paste0(bnames_cft,"_drip"))
bnames_cft <- c(paste0(bnames_cft,"_rf"),paste0(bnames_cft,"_ir"))
nbands_cft <- length(bnames_cft)

bnames_soil <- c("1","2","3","4","5")
nbands_soil <- length(bnames_soil)

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

#read landuse for single for array(cell,cropland area)
read_data_landuse_single_year <- function(filename="",startyear=2005, refyear=1700, NCELL=67420, lu=32)
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

#read landuse input and create matrix(cell,cropland area per year) for LUH2 (32 bands) from 1700-2018
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
  dummy<-cbind(dummy,dummy2,dummy2,dummy2)               #because LU in year 2015=...=2018
  data.out<-dummy[,2:(NYEARS+1)]                         #remove first column
  colnames(data.out)<-c(startyear:(startyear+NYEARS-1))  #name colums
  close(f);rm(f,dummy,dummy2)
  data.out <- data.out*area/1000                         #scaling 1000 for LU
}

###read transient LU:
landuse32bands<-paste0(landuse_path,"lu_madrat_850-2015_32bands.clm")
landuse_cropland_area_32<-read_landuse_cropland_area_to_years_32(landuse32bands)

data_landuse2005_32 <- read_data_landuse_single_year(landuse32bands,startyear=2005,refyear=850,lu=32)
frac_cropland2005_32 <- rowSums(data_landuse2005_32[,-c(rep(14:16,2)+rep(16*(0:1),each=3))])
is_cropland2005_32 <- ifelse(frac_cropland2005_32>0,1,NA)
area_cropland2005_32 <- frac_cropland2005_32*cellarea # in m2

cropcells<-colSums(landuse_cropland_area_32 != 0)
attributes(cropcells)<-NULL
pnvcells<-rep(67420,319)-cropcells

##read output fuctions:
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

read_data_soillayer <- function(filename="",NCELL=67420,NYEARS=1,syear=1700, refyear=2018,size=4,area=cellarea,scaling=1e+15)
{
  data.out <- array(data=0,dim=c(NCELL,nbands_soil,NYEARS))
  f <- file(filename,"rb")
  seek(f,where=NCELL*(syear-refyear)*nbands_soil*size*NYEARS,origin="start")
  data.out<-matrix(readBin(f,double(),n=NCELL*nbands_soil*NYEARS,size), ncol=nbands_soil,byrow=F)
  close(f);rm(f);gc()
  data.out<-data.out*area/scaling
}

read_data_soillayer_singleyear <- function(filename="",NCELL=67420,yr=1700,size=4,area=cellarea,scaling=1e+15)
{
  data.out <- array(data=0,dim=c(NCELL,nbands_soil))
  f <- file(filename,"rb")
  seek(f,where=NCELL*(yr-1700)*nbands_soil*size,origin="start")
  data.out<-matrix(readBin(f,double(),n=NCELL*nbands_soil,size), ncol=nbands_soil,byrow=F)
  close(f);rm(f);gc()
  data.out<-data.out*area/scaling
  colSums(data.out)
}

read_output <- function(filename="", NCELL=67420,syear=1700, refyear=1700, NYEARS=years,size=4)
{
  data.out <- array(data=0,dim=c(NCELL,NYEARS))
  f <- file(filename,"rb")
  seek(f,where=NCELL*(syear-refyear)*size,origin="start")
  data.out<-matrix(readBin(f,double(),n=NCELL*NYEARS,size), ncol=NYEARS,byrow=F)
  colnames(data.out) <- c(syear:(syear+NYEARS-1))
  close(f)
  data.out
}

#calculate pnv area:
sum_cropland2005<-sum(area_cropland2005_32)
sum_cropland2005<-rep(sum_cropland2005,319)
sum_cropland_lu<-colSums(landuse_cropland_area_32)
diff_cropland_2005_lu<-sum_cropland2005-sum_cropland_lu
area_all<-matrix(rep(cellarea,319),ncol=319,byrow=F)
area_pnv<-area_all-landuse_cropland_area_32
area_pnv[area_pnv<=0]<-1e+6

#read cropland data:
npp_agr_wlu<-read_output_global_years(paste0(outp_lu,"/anpp_agr.bin"),area=landuse_cropland_area_32)
soilc_agr_wlu<-read_output_global_years(paste0(outp_lu,"/soilc_agr.bin"),area=landuse_cropland_area_32)
litc_agr_wlu<-read_output_global_years(paste0(outp_lu,"/litc_agr.bin"),area=landuse_cropland_area_32)
allc_agr_wlu<-soilc_agr_wlu+litc_agr_wlu
harv_agr_wlu<-read_output_global_years(paste0(outp_lu,"/flux_harvest.bin"),area=landuse_cropland_area_32)

#calclate ratio harvest to cropland npp
harvest_data<-read_output(paste0(outp_lu,"/flux_harvest.bin"))
harvest_years<-colSums(harvest_data)
npp_data<-read_output(paste0(outp_lu,"/anpp_agr.bin"))
npp_years<-colSums(npp_data)
ratio_harvest_npp<-harvest_years/npp_years

#read const_lu data:
npp_agr_clu<-read_output_global_years(paste0(outp_clu,"/anpp_agr.bin"),area=area_cropland2005_32)
soilc_agr_clu<-read_output_global_years(paste0(outp_clu,"/soilc_agr.bin"),area=area_cropland2005_32)
litc_agr_clu<-read_output_global_years(paste0(outp_clu,"/litc_agr.bin"),area=area_cropland2005_32)
allc_agr_clu<-soilc_agr_clu+litc_agr_clu
npp_clu<-read_output_global_years(paste0(outp_clu,"/anpp.bin"),area=cellarea)
soilc_clu<-read_output_global_years(paste0(outp_clu,"/soilc.bin"),area=cellarea)
litc_clu<-read_output_global_years(paste0(outp_clu,"/litc.bin"),area=cellarea)
allc_clu<-soilc_clu+litc_clu

#read density:
d_soilc<-read_output(paste0(outp_lu,"soilc.bin")) #density soil carbon g/m2
d_litc<-read_output(paste0(outp_lu,"litc.bin"))
d_soilc_agr<-read_output(paste0(outp_lu,"soilc_agr.bin"))
d_litc_agr<-read_output(paste0(outp_lu,"litc_agr.bin"))
d_npp_cells<-read_output(paste0(outp_lu,"anpp.bin")) #density npp g/m2
d_npp_agr_cells<-read_output(paste0(outp_lu,"anpp_agr.bin")) #density npp g/m2

#calculate d_pnv:
soc_soilc<-d_soilc*area_all
soc_litc<-d_litc*area_all
npp_total_cells<-d_npp_cells*area_all
soc_total_cells<-soc_soilc+soc_litc
soc_soilc_agr<-d_soilc_agr*landuse_cropland_area_32
soc_litc_agr<-d_litc_agr*landuse_cropland_area_32
npp_agr_total_cells<-d_npp_agr_cells*landuse_cropland_area_32
soc_agr_total_cells<-soc_soilc_agr+soc_litc_agr
soc_agr_total_global<-colSums(soc_agr_total_cells)/1e+15
npp_agr_total_global<-colSums(npp_agr_total_cells)/1e+15
d_soc_pnv_cells<-((soc_total_cells)-(soc_agr_total_cells))/area_pnv
d_soc_pnv_global<-colSums(d_soc_pnv_cells)/pnvcells
d_npp_pnv_cells<-((npp_total_cells)-(npp_agr_total_cells))/area_pnv
d_npp_pnv_global<-colSums(d_npp_pnv_cells)/pnvcells

soc_pnv_cells<-d_soc_pnv_cells*area_pnv
soc_pnv_global<-colSums(soc_pnv_cells)/1e+15

##calculate soc cultivate: SOC loss form LUC
soc_cultivate<-(d_soc_pnv_global*(diff_cropland_2005_lu)/1e+15)+soc_agr_total_global
npp_cultivate<-(d_npp_pnv_global*(diff_cropland_2005_lu)/1e+15)+npp_agr_total_global
ratio_harvest_npp<-ratio_harvest_npp*100  #convert ratio to %

##plot in pdf:
pdf(file=paste0(plots_path,"plot_historical_cropland_NPP_SOC_figure2.pdf"),width=10,height=4)
par(mfrow=c(1,3))
plot(cbind(c(1700:2018),npp_agr_wlu),type="l",xlim=c(1700,2020),ylim=c(0,5),col="black",xlab="Year",ylab="")
lines(cbind(c(1700:2018),npp_agr_clu),col="turquoise3")
lines(cbind(c(1700:2018),harv_agr_wlu),col="green3")
mtext(text="PgC "~a^-1~"",side = 2, line = 2.2,cex=0.7)
legend("topleft",c("h_dLU cropland NPP","h_dLU harvested C","h_cLU cropland NPP"),col=c("black","green3","turquoise3"),
       lty=1,cex=0.9,bty="n")
mtext("(A)",3,cex=0.7)

plot(cbind(c(1700:2018),ratio_harvest_npp),type="l",ylim=c(10,40),xlab="Year",ylab="Harvested C to cropland NPP (%) in h_dLU")
mtext("(B)",3,cex=0.7)
plot(cbind(c(1700:2018),allc_agr_wlu),type="l",xlim=c(1700,2020),ylim=c(0,500),col="black",xlab = "Year",ylab = "PgC")
lines(cbind(c(1700:2018),soc_cultivate),col="red")
lines(cbind(c(1700:2018),allc_agr_clu),col="turquoise3")
legend("topleft",c("h_dLU cropland SOC","h_dLU_area05 SOC","h_cLU cropland SOC"),
       col=c("black","red","turquoise3"),lty=1,cex=0.9,bty="n")
mtext("(C)",3,cex=0.7)
dev.off()

# calculate soil layer SOC analysis (for table 2):
sc_agr_l_1700<-read_data_soillayer_singleyear(paste0(outp_lu,"soilc_agr_layer.bin"),yr=1700,area=landuse_cropland_area_32[,1])
sc_agr_l_2000<-read_data_soillayer_singleyear(paste0(outp_lu,"soilc_agr_layer.bin"),yr=2000,area=landuse_cropland_area_32[,301])
sc_agr_l_2001<-read_data_soillayer_singleyear(paste0(outp_lu,"soilc_agr_layer.bin"),yr=2001,area=landuse_cropland_area_32[,302])
sc_agr_l_2002<-read_data_soillayer_singleyear(paste0(outp_lu,"soilc_agr_layer.bin"),yr=2002,area=landuse_cropland_area_32[,303])
sc_agr_l_2003<-read_data_soillayer_singleyear(paste0(outp_lu,"soilc_agr_layer.bin"),yr=2003,area=landuse_cropland_area_32[,304])
sc_agr_l_2004<-read_data_soillayer_singleyear(paste0(outp_lu,"soilc_agr_layer.bin"),yr=2004,area=landuse_cropland_area_32[,305])
sc_agr_l_2005<-read_data_soillayer_singleyear(paste0(outp_lu,"soilc_agr_layer.bin"),yr=2005,area=landuse_cropland_area_32[,306])
sc_agr_l_2006<-read_data_soillayer_singleyear(paste0(outp_lu,"soilc_agr_layer.bin"),yr=2006,area=landuse_cropland_area_32[,307])
sc_agr_l_2007<-read_data_soillayer_singleyear(paste0(outp_lu,"soilc_agr_layer.bin"),yr=2007,area=landuse_cropland_area_32[,308])
sc_agr_l_2008<-read_data_soillayer_singleyear(paste0(outp_lu,"soilc_agr_layer.bin"),yr=2008,area=landuse_cropland_area_32[,309])
sc_agr_l_2009<-read_data_soillayer_singleyear(paste0(outp_lu,"soilc_agr_layer.bin"),yr=2009,area=landuse_cropland_area_32[,310])
sc_agr_l_2015<-read_data_soillayer_singleyear(paste0(outp_lu,"soilc_agr_layer.bin"),yr=2015,area=landuse_cropland_area_32[,316])
sc_agr_l_2018<-read_data_soillayer_singleyear(paste0(outp_lu,"soilc_agr_layer.bin"),yr=2018,area=landuse_cropland_area_32[,319])

litc_agr_mean_00_09<-mean(litc_agr_wlu[301:310])

sc_agr_l_00_09<-(sc_agr_l_2000+sc_agr_l_2001+sc_agr_l_2002+sc_agr_l_2003+sc_agr_l_2004+
                   sc_agr_l_2005+sc_agr_l_2006+sc_agr_l_2007+sc_agr_l_2008+sc_agr_l_2009)/10
sc_agr_20<-sc_agr_l_00_09[1]+litc_agr_mean_00_09
sc_agr_30<-sc_agr_l_00_09[1]+(1/3*sc_agr_l_00_09[2])+litc_agr_mean_00_09
sc_agr_50<-sc_agr_l_00_09[1]+sc_agr_l_00_09[2]+litc_agr_mean_00_09
sc_agr_100<-sc_agr_l_00_09[1]+sc_agr_l_00_09[2]+sc_agr_l_00_09[3]+litc_agr_mean_00_09
sc_agr_200<-sc_agr_l_00_09[1]+sc_agr_l_00_09[2]+sc_agr_l_00_09[3]+sc_agr_l_00_09[4]+litc_agr_mean_00_09
sc_agr_300<-sc_agr_l_00_09[1]+sc_agr_l_00_09[2]+sc_agr_l_00_09[3]+sc_agr_l_00_09[4]+sc_agr_l_00_09[5]+litc_agr_mean_00_09

sc_l_2000<-read_data_soillayer_singleyear(paste0(outp_lu,"soilc_layer.bin"),yr=2000,area=cellarea)
sc_l_2001<-read_data_soillayer_singleyear(paste0(outp_lu,"soilc_layer.bin"),yr=2001,area=cellarea)
sc_l_2002<-read_data_soillayer_singleyear(paste0(outp_lu,"soilc_layer.bin"),yr=2002,area=cellarea)
sc_l_2003<-read_data_soillayer_singleyear(paste0(outp_lu,"soilc_layer.bin"),yr=2003,area=cellarea)
sc_l_2004<-read_data_soillayer_singleyear(paste0(outp_lu,"soilc_layer.bin"),yr=2004,area=cellarea)
sc_l_2005<-read_data_soillayer_singleyear(paste0(outp_lu,"soilc_layer.bin"),yr=2005,area=cellarea)
sc_l_2006<-read_data_soillayer_singleyear(paste0(outp_lu,"soilc_layer.bin"),yr=2006,area=cellarea)
sc_l_2007<-read_data_soillayer_singleyear(paste0(outp_lu,"soilc_layer.bin"),yr=2007,area=cellarea)
sc_l_2008<-read_data_soillayer_singleyear(paste0(outp_lu,"soilc_layer.bin"),yr=2008,area=cellarea)
sc_l_2009<-read_data_soillayer_singleyear(paste0(outp_lu,"soilc_layer.bin"),yr=2009,area=cellarea)

litc_mean_00_09<-mean(litc_wlu[301:310])

sc_l_00_09<-(sc_l_2000+sc_l_2001+sc_l_2002+sc_l_2003+sc_l_2004+
               sc_l_2005+sc_l_2006+sc_l_2007+sc_l_2008+sc_l_2009)/10
sc_20<-sc_l_00_09[1]+litc_mean_00_09
sc_30<-sc_l_00_09[1]+(1/3*sc_l_00_09[2])+litc_mean_00_09
sc_100<-sc_l_00_09[1]+sc_l_00_09[2]+sc_l_00_09[3]+litc_mean_00_09
sc_200<-sc_l_00_09[1]+sc_l_00_09[2]+sc_l_00_09[3]+sc_l_00_09[4]+litc_mean_00_09
sc_300<-sc_l_00_09[1]+sc_l_00_09[2]+sc_l_00_09[3]+sc_l_00_09[4]+sc_l_00_09[5]+litc_mean_00_09

scp_l_2000<-read_data_soillayer_singleyear(paste0(outp_pnv,"soilc_layer.bin"),yr=2000,area=cellarea)
scp_l_2001<-read_data_soillayer_singleyear(paste0(outp_pnv,"soilc_layer.bin"),yr=2001,area=cellarea)
scp_l_2002<-read_data_soillayer_singleyear(paste0(outp_pnv,"soilc_layer.bin"),yr=2002,area=cellarea)
scp_l_2003<-read_data_soillayer_singleyear(paste0(outp_pnv,"soilc_layer.bin"),yr=2003,area=cellarea)
scp_l_2004<-read_data_soillayer_singleyear(paste0(outp_pnv,"soilc_layer.bin"),yr=2004,area=cellarea)
scp_l_2005<-read_data_soillayer_singleyear(paste0(outp_pnv,"soilc_layer.bin"),yr=2005,area=cellarea)
scp_l_2006<-read_data_soillayer_singleyear(paste0(outp_pnv,"soilc_layer.bin"),yr=2006,area=cellarea)
scp_l_2007<-read_data_soillayer_singleyear(paste0(outp_pnv,"soilc_layer.bin"),yr=2007,area=cellarea)
scp_l_2008<-read_data_soillayer_singleyear(paste0(outp_pnv,"soilc_layer.bin"),yr=2008,area=cellarea)
scp_l_2009<-read_data_soillayer_singleyear(paste0(outp_pnv,"soilc_layer.bin"),yr=2009,area=cellarea)

litc_pnv<-read_output_global_years(paste0(outp_pnv,"/litc.bin"),area=cellarea)
litc_mean_pnv_00_09<-mean(litc_pnv[301:310])

scp_l_00_09<-(scp_l_2000+scp_l_2001+scp_l_2002+scp_l_2003+scp_l_2004+
                scp_l_2005+scp_l_2006+scp_l_2007+scp_l_2008+scp_l_2009)/10
scp_20<-scp_l_00_09[1]+litc_mean_pnv_00_09
scp_30<-scp_l_00_09[1]+(1/3*sc_l_00_09[2])+litc_mean_pnv_00_09
scp_100<-scp_l_00_09[1]+scp_l_00_09[2]+scp_l_00_09[3]+litc_mean_pnv_00_09
scp_200<-scp_l_00_09[1]+scp_l_00_09[2]+scp_l_00_09[3]+scp_l_00_09[4]+litc_mean_pnv_00_09
scp_300<-scp_l_00_09[1]+scp_l_00_09[2]+scp_l_00_09[3]+scp_l_00_09[4]+scp_l_00_09[5]+litc_mean_pnv_00_09
