###############################################################################################
#                                                                                             #
#  Script by Tobias Herzfeld (04.2020) for:                                                   #
# "SOC sequestration potentials for agricultural management practices under climate change"   #
#                                                                                             #
#  reads in the csv files created by the script analysis_convert_data2csv_rcps.R from         #
#  LPJmL5.0-tillage2 model runs and creates plots using ggplot for figure 3                   #
#                                                                                             #
###############################################################################################

rm(list=ls(all=T));gc()

library(ggplot2)
library(plotrix)
library(dplyr)

##define paths:
output_path<-"path/.../"    #path of csv output data
plots_path<-"path/.../"     #path for output plots created

##read data:
sim1<-"OUTPUT_HIST_HADGEM_7000"
sim2<-"OUTPUT_HIST_GFDL_7000"
sim3<-"OUTPUT_HIST_IPSL_7000"
sim4<-"OUTPUT_HIST_MIROC_7000"
df1<-read.csv(paste0(output_path,sim1,"_LU32.csv"))
df2<-read.csv(paste0(output_path,sim2,"_LU32.csv"))
df3<-read.csv(paste0(output_path,sim3,"_LU32.csv"))
df4<-read.csv(paste0(output_path,sim4,"_LU32.csv"))
HIST_RCP<-(df1+df2+df3+df4)/4
rm(df1,df2,df3,df4)

sim1<-"OUTPUT_TRc05_rcp26_hadgem"
sim2<-"OUTPUT_TRc05_rcp26_gfdl"
sim3<-"OUTPUT_TRc05_rcp26_ipsl"
sim4<-"OUTPUT_TRc05_rcp26_miroc"
sim5<-"OUTPUT_TRc05_rcp85_hadgem"
sim6<-"OUTPUT_TRc05_rcp85_gfdl"
sim7<-"OUTPUT_TRc05_rcp85_ipsl"
sim8<-"OUTPUT_TRc05_rcp85_miroc"
df1<-read.csv(paste0(output_path,sim1,".csv"))
df2<-read.csv(paste0(output_path,sim2,".csv"))
df3<-read.csv(paste0(output_path,sim3,".csv"))
df4<-read.csv(paste0(output_path,sim4,".csv"))
CLU_rcp26<-(df1+df2+df3+df4)/4

rcp26_CLU_min_npp_agr<-min(df1$npp_agr[94],df2$npp_agr[94],df3$npp_agr[94],df4$npp_agr[94])
rcp26_CLU_max_npp_agr<-max(df1$npp_agr[94],df2$npp_agr[94],df3$npp_agr[94],df4$npp_agr[94])
rcp26_CLU_min_allc_agr<-min(df1$allc_agr[94],df2$allc_agr[94],df3$allc_agr[94],df4$allc_agr[94])
rcp26_CLU_max_allc_agr<-max(df1$allc_agr[94],df2$allc_agr[94],df3$allc_agr[94],df4$allc_agr[94])
rcp26_CLU_min_mrt_agr<-100/min(df1$mrt_agr[94],df2$mrt_agr[94],df3$mrt_agr[94],df4$mrt_agr[94])
rcp26_CLU_max_mrt_agr<-100/max(df1$mrt_agr[94],df2$mrt_agr[94],df3$mrt_agr[94],df4$mrt_agr[94])
rcp26_CLU_min_litf_agr<-min(df1$litf_agr[94],df2$litf_agr[94],df3$litf_agr[94],df4$litf_agr[94])
rcp26_CLU_max_litf_agr<-max(df1$litf_agr[94],df2$litf_agr[94],df3$litf_agr[94],df4$litf_agr[94])
rm(df1,df2,df3,df4)

df1<-read.csv(paste0(output_path,sim5,".csv"))
df2<-read.csv(paste0(output_path,sim6,".csv"))
df3<-read.csv(paste0(output_path,sim7,".csv"))
df4<-read.csv(paste0(output_path,sim8,".csv"))
CLU_rcp85<-(df1+df2+df3+df4)/4

rcp85_CLU_min_npp_agr<-min(df1$npp_agr[94],df2$npp_agr[94],df3$npp_agr[94],df4$npp_agr[94])
rcp85_CLU_max_npp_agr<-max(df1$npp_agr[94],df2$npp_agr[94],df3$npp_agr[94],df4$npp_agr[94])
rcp85_CLU_min_allc_agr<-min(df1$allc_agr[94],df2$allc_agr[94],df3$allc_agr[94],df4$allc_agr[94])
rcp85_CLU_max_allc_agr<-max(df1$allc_agr[94],df2$allc_agr[94],df3$allc_agr[94],df4$allc_agr[94])
rcp85_CLU_min_mrt_agr<-100/min(df1$mrt_agr[94],df2$mrt_agr[94],df3$mrt_agr[94],df4$mrt_agr[94])
rcp85_CLU_max_mrt_agr<-100/max(df1$mrt_agr[94],df2$mrt_agr[94],df3$mrt_agr[94],df4$mrt_agr[94])
rcp85_CLU_min_litf_agr<-min(df1$litf_agr[94],df2$litf_agr[94],df3$litf_agr[94],df4$litf_agr[94])
rcp85_CLU_max_litf_agr<-max(df1$litf_agr[94],df2$litf_agr[94],df3$litf_agr[94],df4$litf_agr[94])
rm(df1,df2,df3,df4)

sim1<-"rcp26_hadgem"
sim2<-"rcp26_gfdl"
sim3<-"rcp26_ipsl"
sim4<-"rcp26_miroc"
sim5<-"rcp85_hadgem"
sim6<-"rcp85_gfdl"
sim7<-"rcp85_ipsl"
sim8<-"rcp85_miroc"

TR1 <- read.csv(paste0(output_path,"OUTPUT_T_R_df_",sim1,".csv"))
TR2 <- read.csv(paste0(output_path,"OUTPUT_T_R_df_",sim2,".csv"))
TR3 <- read.csv(paste0(output_path,"OUTPUT_T_R_df_",sim3,".csv"))
TR4 <- read.csv(paste0(output_path,"OUTPUT_T_R_df_",sim4,".csv"))
TR_rcp26<-(TR1+TR2+TR3+TR4)/4

rcp26_TR_min_npp_agr<-min(TR1$npp_agr[94],TR2$npp_agr[94],TR3$npp_agr[94],TR4$npp_agr[94])
rcp26_TR_max_npp_agr<-max(TR1$npp_agr[94],TR2$npp_agr[94],TR3$npp_agr[94],TR4$npp_agr[94])
rcp26_TR_min_allc_agr<-min(TR1$allc_agr[94],TR2$allc_agr[94],TR3$allc_agr[94],TR4$allc_agr[94])
rcp26_TR_max_allc_agr<-max(TR1$allc_agr[94],TR2$allc_agr[94],TR3$allc_agr[94],TR4$allc_agr[94])
rcp26_TR_min_mrt_agr<-100/min(TR1$mrt_agr[94],TR2$mrt_agr[94],TR3$mrt_agr[94],TR4$mrt_agr[94])
rcp26_TR_max_mrt_agr<-100/max(TR1$mrt_agr[94],TR2$mrt_agr[94],TR3$mrt_agr[94],TR4$mrt_agr[94])
rcp26_TR_min_litf_agr<-min(TR1$litf_agr[94],TR2$litf_agr[94],TR3$litf_agr[94],TR4$litf_agr[94])
rcp26_TR_max_litf_agr<-max(TR1$litf_agr[94],TR2$litf_agr[94],TR3$litf_agr[94],TR4$litf_agr[94])
rm(TR1,TR2,TR3,TR4)

TR1 <- read.csv(paste0(output_path,"OUTPUT_T_R_df_",sim5,".csv"))
TR2 <- read.csv(paste0(output_path,"OUTPUT_T_R_df_",sim6,".csv"))
TR3 <- read.csv(paste0(output_path,"OUTPUT_T_R_df_",sim7,".csv"))
TR4 <- read.csv(paste0(output_path,"OUTPUT_T_R_df_",sim8,".csv"))
TR_rcp85<-(TR1+TR2+TR3+TR4)/4

rcp85_TR_min_npp_agr<-min(TR1$npp_agr[94],TR2$npp_agr[94],TR3$npp_agr[94],TR4$npp_agr[94])
rcp85_TR_max_npp_agr<-max(TR1$npp_agr[94],TR2$npp_agr[94],TR3$npp_agr[94],TR4$npp_agr[94])
rcp85_TR_min_allc_agr<-min(TR1$allc_agr[94],TR2$allc_agr[94],TR3$allc_agr[94],TR4$allc_agr[94])
rcp85_TR_max_allc_agr<-max(TR1$allc_agr[94],TR2$allc_agr[94],TR3$allc_agr[94],TR4$allc_agr[94])
rcp85_TR_min_mrt_agr<-100/min(TR1$mrt_agr[94],TR2$mrt_agr[94],TR3$mrt_agr[94],TR4$mrt_agr[94])
rcp85_TR_max_mrt_agr<-100/max(TR1$mrt_agr[94],TR2$mrt_agr[94],TR3$mrt_agr[94],TR4$mrt_agr[94])
rcp85_TR_min_litf_agr<-min(TR1$litf_agr[94],TR2$litf_agr[94],TR3$litf_agr[94],TR4$litf_agr[94])
rcp85_TR_max_litf_agr<-max(TR1$litf_agr[94],TR2$litf_agr[94],TR3$litf_agr[94],TR4$litf_agr[94])
rm(TR1,TR2,TR3,TR4)

NTR1 <- read.csv(paste0(output_path,"OUTPUT_NT_R_df_",sim1,".csv"))
NTR2 <- read.csv(paste0(output_path,"OUTPUT_NT_R_df_",sim2,".csv"))
NTR3 <- read.csv(paste0(output_path,"OUTPUT_NT_R_df_",sim3,".csv"))
NTR4 <- read.csv(paste0(output_path,"OUTPUT_NT_R_df_",sim4,".csv"))
NTR_rcp26<-(NTR1+NTR2+NTR3+NTR4)/4

rcp26_NTR_min_npp_agr<-min(NTR1$npp_agr[94],NTR2$npp_agr[94],NTR3$npp_agr[94],NTR4$npp_agr[94])
rcp26_NTR_max_npp_agr<-max(NTR1$npp_agr[94],NTR2$npp_agr[94],NTR3$npp_agr[94],NTR4$npp_agr[94])
rcp26_NTR_min_allc_agr<-min(NTR1$allc_agr[94],NTR2$allc_agr[94],NTR3$allc_agr[94],NTR4$allc_agr[94])
rcp26_NTR_max_allc_agr<-max(NTR1$allc_agr[94],NTR2$allc_agr[94],NTR3$allc_agr[94],NTR4$allc_agr[94])
rcp26_NTR_min_mrt_agr<-100/min(NTR1$mrt_agr[94],NTR2$mrt_agr[94],NTR3$mrt_agr[94],NTR4$mrt_agr[94])
rcp26_NTR_max_mrt_agr<-100/max(NTR1$mrt_agr[94],NTR2$mrt_agr[94],NTR3$mrt_agr[94],NTR4$mrt_agr[94])
rcp26_NTR_min_litf_agr<-min(NTR1$litf_agr[94],NTR2$litf_agr[94],NTR3$litf_agr[94],NTR4$litf_agr[94])
rcp26_NTR_max_litf_agr<-max(NTR1$litf_agr[94],NTR2$litf_agr[94],NTR3$litf_agr[94],NTR4$litf_agr[94])

NTR1 <- read.csv(paste0(output_path,"OUTPUT_NT_R_df_",sim5,".csv"))
NTR2 <- read.csv(paste0(output_path,"OUTPUT_NT_R_df_",sim6,".csv"))
NTR3 <- read.csv(paste0(output_path,"OUTPUT_NT_R_df_",sim7,".csv"))
NTR4 <- read.csv(paste0(output_path,"OUTPUT_NT_R_df_",sim8,".csv"))
NTR_rcp85<-(NTR1+NTR2+NTR3+NTR4)/4

rcp85_NTR_min_npp_agr<-min(NTR1$npp_agr[94],NTR2$npp_agr[94],NTR3$npp_agr[94],NTR4$npp_agr[94])
rcp85_NTR_max_npp_agr<-max(NTR1$npp_agr[94],NTR2$npp_agr[94],NTR3$npp_agr[94],NTR4$npp_agr[94])
rcp85_NTR_min_allc_agr<-min(NTR1$allc_agr[94],NTR2$allc_agr[94],NTR3$allc_agr[94],NTR4$allc_agr[94])
rcp85_NTR_max_allc_agr<-max(NTR1$allc_agr[94],NTR2$allc_agr[94],NTR3$allc_agr[94],NTR4$allc_agr[94])
rcp85_NTR_min_mrt_agr<-100/min(NTR1$mrt_agr[94],NTR2$mrt_agr[94],NTR3$mrt_agr[94],NTR4$mrt_agr[94])
rcp85_NTR_max_mrt_agr<-100/max(NTR1$mrt_agr[94],NTR2$mrt_agr[94],NTR3$mrt_agr[94],NTR4$mrt_agr[94])
rcp85_NTR_min_litf_agr<-min(NTR1$litf_agr[94],NTR2$litf_agr[94],NTR3$litf_agr[94],NTR4$litf_agr[94])
rcp85_NTR_max_litf_agr<-max(NTR1$litf_agr[94],NTR2$litf_agr[94],NTR3$litf_agr[94],NTR4$litf_agr[94])
rm(NTR1,NTR2,NTR3,NTR4)

TNR1 <- read.csv(paste0(output_path,"OUTPUT_T_NR_df_",sim1,".csv"))
TNR2 <- read.csv(paste0(output_path,"OUTPUT_T_NR_df_",sim2,".csv"))
TNR3 <- read.csv(paste0(output_path,"OUTPUT_T_NR_df_",sim3,".csv"))
TNR4 <- read.csv(paste0(output_path,"OUTPUT_T_NR_df_",sim4,".csv"))
TNR_rcp26<-(TNR1+TNR2+TNR3+TNR4)/4

rcp26_TNR_min_npp_agr<-min(TNR1$npp_agr[94],TNR2$npp_agr[94],TNR3$npp_agr[94],TNR4$npp_agr[94])
rcp26_TNR_max_npp_agr<-max(TNR1$npp_agr[94],TNR2$npp_agr[94],TNR3$npp_agr[94],TNR4$npp_agr[94])
rcp26_TNR_min_allc_agr<-min(TNR1$allc_agr[94],TNR2$allc_agr[94],TNR3$allc_agr[94],TNR4$allc_agr[94])
rcp26_TNR_max_allc_agr<-max(TNR1$allc_agr[94],TNR2$allc_agr[94],TNR3$allc_agr[94],TNR4$allc_agr[94])
rcp26_TNR_min_mrt_agr<-100/min(TNR1$mrt_agr[94],TNR2$mrt_agr[94],TNR3$mrt_agr[94],TNR4$mrt_agr[94])
rcp26_TNR_max_mrt_agr<-100/max(TNR1$mrt_agr[94],TNR2$mrt_agr[94],TNR3$mrt_agr[94],TNR4$mrt_agr[94])
rcp26_TNR_min_litf_agr<-min(TNR1$litf_agr[94],TNR2$litf_agr[94],TNR3$litf_agr[94],TNR4$litf_agr[94])
rcp26_TNR_max_litf_agr<-max(TNR1$litf_agr[94],TNR2$litf_agr[94],TNR3$litf_agr[94],TNR4$litf_agr[94])

TNR1 <- read.csv(paste0(output_path,"OUTPUT_T_NR_df_",sim5,".csv"))
TNR2 <- read.csv(paste0(output_path,"OUTPUT_T_NR_df_",sim6,".csv"))
TNR3 <- read.csv(paste0(output_path,"OUTPUT_T_NR_df_",sim7,".csv"))
TNR4 <- read.csv(paste0(output_path,"OUTPUT_T_NR_df_",sim8,".csv"))
TNR_rcp85<-(TNR1+TNR2+TNR3+TNR4)/4

rcp85_TNR_min_npp_agr<-min(TNR1$npp_agr[94],TNR2$npp_agr[94],TNR3$npp_agr[94],TNR4$npp_agr[94])
rcp85_TNR_max_npp_agr<-max(TNR1$npp_agr[94],TNR2$npp_agr[94],TNR3$npp_agr[94],TNR4$npp_agr[94])
rcp85_TNR_min_allc_agr<-min(TNR1$allc_agr[94],TNR2$allc_agr[94],TNR3$allc_agr[94],TNR4$allc_agr[94])
rcp85_TNR_max_allc_agr<-max(TNR1$allc_agr[94],TNR2$allc_agr[94],TNR3$allc_agr[94],TNR4$allc_agr[94])
rcp85_TNR_min_mrt_agr<-100/min(TNR1$mrt_agr[94],TNR2$mrt_agr[94],TNR3$mrt_agr[94],TNR4$mrt_agr[94])
rcp85_TNR_max_mrt_agr<-100/max(TNR1$mrt_agr[94],TNR2$mrt_agr[94],TNR3$mrt_agr[94],TNR4$mrt_agr[94])
rcp85_TNR_min_litf_agr<-min(TNR1$litf_agr[94],TNR2$litf_agr[94],TNR3$litf_agr[94],TNR4$litf_agr[94])
rcp85_TNR_max_litf_agr<-max(TNR1$litf_agr[94],TNR2$litf_agr[94],TNR3$litf_agr[94],TNR4$litf_agr[94])

rm(TNR1,TNR2,TNR3,TNR4)

NTNR1 <- read.csv(paste0(output_path,"OUTPUT_NT_NR_df_",sim1,".csv"))
NTNR2 <- read.csv(paste0(output_path,"OUTPUT_NT_NR_df_",sim2,".csv"))
NTNR3 <- read.csv(paste0(output_path,"OUTPUT_NT_NR_df_",sim3,".csv"))
NTNR4 <- read.csv(paste0(output_path,"OUTPUT_NT_NR_df_",sim4,".csv"))
NTNR_rcp26<-(NTNR1+NTNR2+NTNR3+NTNR4)/4

rcp26_NTNR_min_npp_agr<-min(NTNR1$npp_agr[94],NTNR2$npp_agr[94],NTNR3$npp_agr[94],NTNR4$npp_agr[94])
rcp26_NTNR_max_npp_agr<-max(NTNR1$npp_agr[94],NTNR2$npp_agr[94],NTNR3$npp_agr[94],NTNR4$npp_agr[94])
rcp26_NTNR_min_allc_agr<-min(NTNR1$allc_agr[94],NTNR2$allc_agr[94],NTNR3$allc_agr[94],NTNR4$allc_agr[94])
rcp26_NTNR_max_allc_agr<-max(NTNR1$allc_agr[94],NTNR2$allc_agr[94],NTNR3$allc_agr[94],NTNR4$allc_agr[94])
rcp26_NTNR_min_mrt_agr<-100/min(NTNR1$mrt_agr[94],NTNR2$mrt_agr[94],NTNR3$mrt_agr[94],NTNR4$mrt_agr[94])
rcp26_NTNR_max_mrt_agr<-100/max(NTNR1$mrt_agr[94],NTNR2$mrt_agr[94],NTNR3$mrt_agr[94],NTNR4$mrt_agr[94])
rcp26_NTNR_min_litf_agr<-min(NTNR1$litf_agr[94],NTNR2$litf_agr[94],NTNR3$litf_agr[94],NTNR4$litf_agr[94])
rcp26_NTNR_max_litf_agr<-max(NTNR1$litf_agr[94],NTNR2$litf_agr[94],NTNR3$litf_agr[94],NTNR4$litf_agr[94])

NTNR1 <- read.csv(paste0(output_path,"OUTPUT_NT_NR_df_",sim5,".csv"))
NTNR2 <- read.csv(paste0(output_path,"OUTPUT_NT_NR_df_",sim6,".csv"))
NTNR3 <- read.csv(paste0(output_path,"OUTPUT_NT_NR_df_",sim7,".csv"))
NTNR4 <- read.csv(paste0(output_path,"OUTPUT_NT_NR_df_",sim8,".csv"))
NTNR_rcp85<-(NTNR1+NTNR2+NTNR3+NTNR4)/4

rcp85_NTNR_min_npp_agr<-min(NTNR1$npp_agr[94],NTNR2$npp_agr[94],NTNR3$npp_agr[94],NTNR4$npp_agr[94])
rcp85_NTNR_max_npp_agr<-max(NTNR1$npp_agr[94],NTNR2$npp_agr[94],NTNR3$npp_agr[94],NTNR4$npp_agr[94])
rcp85_NTNR_min_allc_agr<-min(NTNR1$allc_agr[94],NTNR2$allc_agr[94],NTNR3$allc_agr[94],NTNR4$allc_agr[94])
rcp85_NTNR_max_allc_agr<-max(NTNR1$allc_agr[94],NTNR2$allc_agr[94],NTNR3$allc_agr[94],NTNR4$allc_agr[94])
rcp85_NTNR_min_mrt_agr<-100/min(NTNR1$mrt_agr[94],NTNR2$mrt_agr[94],NTNR3$mrt_agr[94],NTNR4$mrt_agr[94])
rcp85_NTNR_max_mrt_agr<-100/max(NTNR1$mrt_agr[94],NTNR2$mrt_agr[94],NTNR3$mrt_agr[94],NTNR4$mrt_agr[94])
rcp85_NTNR_min_litf_agr<-min(NTNR1$litf_agr[94],NTNR2$litf_agr[94],NTNR3$litf_agr[94],NTNR4$litf_agr[94])
rcp85_NTNR_max_litf_agr<-max(NTNR1$litf_agr[94],NTNR2$litf_agr[94],NTNR3$litf_agr[94],NTNR4$litf_agr[94])
rm(NTNR1,NTNR2,NTNR3,NTNR4)

syear<-2006
eyear<-2099
scaling<-1

#history included
total_TR26<-rbind(HIST_RCP,TR_rcp26)
TR26_100<-total_TR26[301:400,]

total_NTR26<-rbind(HIST_RCP,NTR_rcp26)
NTR26_100<-total_NTR26[301:400,]

total_TNR26<-rbind(HIST_RCP,TNR_rcp26)
TNR26_100<-total_TNR26[301:400,]

total_NTNR26<-rbind(HIST_RCP,NTNR_rcp26)
NTNR26_100<-total_NTNR26[301:400,]

total_TR85<-rbind(HIST_RCP,TR_rcp85)
TR85_100<-total_TR85[301:400,]

total_NTR85<-rbind(HIST_RCP,NTR_rcp85)
NTR85_100<-total_NTR85[301:400,]

total_TNR85<-rbind(HIST_RCP,TNR_rcp85)
TNR85_100<-total_TNR85[301:400,]

total_NTNR85<-rbind(HIST_RCP,NTNR_rcp85)
NTNR85_100<-total_NTNR85[301:400,]

total_CLU26<-rbind(HIST_RCP,CLU_rcp26)
CLU26_100<-total_CLU26[301:400,]
total_CLU85<-rbind(HIST_RCP,CLU_rcp85)
CLU85_100<-total_CLU85[301:400,]

#convert mrt_agr
TR26_100$mrt_agr<-100/TR26_100$mrt_agr
NTR26_100$mrt_agr<-100/NTR26_100$mrt_agr
TNR26_100$mrt_agr<-100/TNR26_100$mrt_agr
NTNR26_100$mrt_agr<-100/NTNR26_100$mrt_agr
TR85_100$mrt_agr<-100/TR85_100$mrt_agr
NTR85_100$mrt_agr<-100/NTR85_100$mrt_agr
TNR85_100$mrt_agr<-100/TNR85_100$mrt_agr
NTNR85_100$mrt_agr<-100/NTNR85_100$mrt_agr
CLU26_100$mrt_agr<-100/CLU26_100$mrt_agr
CLU85_100$mrt_agr<-100/CLU85_100$mrt_agr

#differences:
diff_max<-TR85_100$allc_agr[100]-NTNR85_100$allc_agr[100]
diff_min<-TR26_100$allc_agr[100]-TNR85_100$allc_agr[100]

#relative change in 2099
rel.change_TR_NTR_26<-((TR26_100$allc_agr[100]-NTR26_100$allc_agr[100])/TR26_100$allc_agr[100])*100
rel.change_TR_NTR_85<-((TR85_100$allc_agr[100]-NTR85_100$allc_agr[100])/TR85_100$allc_agr[100])*100
rel.change_TNR_NTNR_26<-((TNR26_100$allc_agr[100]-NTNR26_100$allc_agr[100])/TNR26_100$allc_agr[100])*100
rel.change_TNR_NTNR_85<-((TNR85_100$allc_agr[100]-NTNR85_100$allc_agr[100])/TNR85_100$allc_agr[100])*100

rel.change_TR_TNR_26<-((TR26_100$allc_agr[100]-TNR26_100$allc_agr[100])/TR26_100$allc_agr[100])*100
rel.change_TR_TNR_85<-((TR85_100$allc_agr[100]-TNR85_100$allc_agr[100])/TR85_100$allc_agr[100])*100
rel.change_NTR_NTNR_26<-((NTR26_100$allc_agr[100]-NTNR26_100$allc_agr[100])/NTR26_100$allc_agr[100])*100
rel.change_NTR_NTNR_85<-((NTR85_100$allc_agr[100]-NTNR85_100$allc_agr[100])/NTR85_100$allc_agr[100])*100

#rel change between rcps with same management:
rel.change_TR_85_26<-((TR85_100$allc_agr[100]-TR26_100$allc_agr[100])/TR26_100$allc_agr[100])*100
rel.change_NTR_85_26<-((NTR85_100$allc_agr[100]-NTR26_100$allc_agr[100])/NTR26_100$allc_agr[100])*100
rel.change_TNR_85_26<-((TNR85_100$allc_agr[100]-TNR26_100$allc_agr[100])/TNR26_100$allc_agr[100])*100
rel.change_NTNR_85_26<-((NTNR85_100$allc_agr[100]-NTNR26_100$allc_agr[100])/NTNR26_100$allc_agr[100])*100

#absolute and relative changes 2099-2005
abs.change_TR26<-TR26_100$allc_agr[100]-TR26_100$allc_agr[6]
abs.change_NTR26<-NTR26_100$allc_agr[100]-NTR26_100$allc_agr[6]
abs.change_TNR26<-TNR26_100$allc_agr[100]-TNR26_100$allc_agr[6]
abs.change_NTNR26<-NTNR26_100$allc_agr[100]-NTNR26_100$allc_agr[6]
rel.change_TR26_2005_2099<-((TR26_100$allc_agr[100]-TR26_100$allc_agr[6])/TR26_100$allc_agr[6])*100
rel.change_NTR26_2005_2099<-((NTR26_100$allc_agr[100]-NTR26_100$allc_agr[6])/NTR26_100$allc_agr[6])*100
rel.change_TNR26_2005_2099<-((TNR26_100$allc_agr[100]-TNR26_100$allc_agr[6])/TNR26_100$allc_agr[6])*100
rel.change_NTNR26_2005_2099<-((NTNR26_100$allc_agr[100]-NTNR26_100$allc_agr[6])/NTNR26_100$allc_agr[6])*100

abs.change_TR85<-TR85_100$allc_agr[100]-TR85_100$allc_agr[6]
abs.change_NTR85<-NTR85_100$allc_agr[100]-NTR85_100$allc_agr[6]
abs.change_TNR85<-TNR85_100$allc_agr[100]-TNR85_100$allc_agr[6]
abs.change_NTNR85<-NTNR85_100$allc_agr[100]-NTNR85_100$allc_agr[6]
rel.change_TR85_2005_2099<-((TR85_100$allc_agr[100]-TR85_100$allc_agr[6])/TR85_100$allc_agr[6])*100
rel.change_NTR85_2005_2099<-((NTR85_100$allc_agr[100]-NTR85_100$allc_agr[6])/NTR85_100$allc_agr[6])*100
rel.change_TNR85_2005_2099<-((TNR85_100$allc_agr[100]-TNR85_100$allc_agr[6])/TNR85_100$allc_agr[6])*100
rel.change_NTNR85_2005_2099<-((NTNR85_100$allc_agr[100]-NTNR85_100$allc_agr[6])/NTNR85_100$allc_agr[6])*100

#absolute and relative changes 2099-2005
abs.change_TR26_2006_2099<-TR26_100$allc_agr[100]-TR26_100$allc_agr[7]
abs.change_NTR26_2006_2099<-NTR26_100$allc_agr[100]-NTR26_100$allc_agr[7]
abs.change_TNR26_2006_2099<-TNR26_100$allc_agr[100]-TNR26_100$allc_agr[7]
abs.change_NTNR26_2006_2099<-NTNR26_100$allc_agr[100]-NTNR26_100$allc_agr[7]
rel.change_TR26_2006_2099<-((TR26_100$allc_agr[100]-TR26_100$allc_agr[7])/TR26_100$allc_agr[7])*100
rel.change_NTR26_2006_2099<-((NTR26_100$allc_agr[100]-NTR26_100$allc_agr[7])/NTR26_100$allc_agr[7])*100
rel.change_TNR26_2006_2099<-((TNR26_100$allc_agr[100]-TNR26_100$allc_agr[7])/TNR26_100$allc_agr[7])*100
rel.change_NTNR26_2006_2099<-((NTNR26_100$allc_agr[100]-NTNR26_100$allc_agr[7])/NTNR26_100$allc_agr[7])*100

abs.change_TR85_2006_2099<-TR85_100$allc_agr[100]-TR85_100$allc_agr[7]
abs.change_NTR85_2006_2099<-NTR85_100$allc_agr[100]-NTR85_100$allc_agr[7]
abs.change_TNR85_2006_2099<-TNR85_100$allc_agr[100]-TNR85_100$allc_agr[7]
abs.change_NTNR85_2006_2099<-NTNR85_100$allc_agr[100]-NTNR85_100$allc_agr[7]
rel.change_TR85_2006_2099<-((TR85_100$allc_agr[100]-TR85_100$allc_agr[7])/TR85_100$allc_agr[7])*100
rel.change_NTR85_2006_2099<-((NTR85_100$allc_agr[100]-NTR85_100$allc_agr[7])/NTR85_100$allc_agr[7])*100
rel.change_TNR85_2006_2099<-((TNR85_100$allc_agr[100]-TNR85_100$allc_agr[7])/TNR85_100$allc_agr[7])*100
rel.change_NTNR85_2006_2099<-((NTNR85_100$allc_agr[100]-NTNR85_100$allc_agr[7])/NTNR85_100$allc_agr[7])*100

abs.change_CLU26<-CLU26_100$allc_agr[100]-CLU26_100$allc_agr[6]
abs.change_CLU85<-CLU85_100$allc_agr[100]-CLU85_100$allc_agr[6]
rel.change_CLU26<-((CLU26_100$allc_agr[100]-CLU26_100$allc_agr[6])/CLU26_100$allc_agr[6])*100
rel.change_CLU85<-((CLU85_100$allc_agr[100]-CLU85_100$allc_agr[6])/CLU85_100$allc_agr[6])*100

abs.change_CLU26_2006_2099<-CLU26_100$allc_agr[100]-CLU26_100$allc_agr[7]
abs.change_CLU85_2006_2099<-CLU85_100$allc_agr[100]-CLU85_100$allc_agr[7]
rel.change_CLU26_2006_2099<-((CLU26_100$allc_agr[100]-CLU26_100$allc_agr[7])/CLU26_100$allc_agr[7])*100
rel.change_CLU85_2006_2099<-((CLU85_100$allc_agr[100]-CLU85_100$allc_agr[7])/CLU85_100$allc_agr[7])*100

#plot data in pdf
library(plotrix)
require(gridExtra)
library(grid)

pdf(paste(plots_path,"ggplot_global_values_RCPs_CONST_LU_AVERAGES.pdf",sep=""), width=8, height=8)
LegendTitle<-"Scenario"

a <-ggplot(TR26_100,aes(x=X,y=npp_agr)) +
  ggtitle("(A)")+
  geom_line(aes(y=TR26_100$npp_agr,color="TR_26",linetype="TR_26"))+
  geom_errorbar(aes(x=2100, ymin=rcp26_TR_min_npp_agr, ymax=rcp26_TR_max_npp_agr), size=0.2, width=2, color="green4",linetype="dotdash")+
  geom_line(aes(y=NTR26_100$npp_agr,color="NTR_26",linetype="NTR_26"))+
  geom_errorbar(aes(x=2102, ymin=rcp26_NTR_min_npp_agr, ymax=rcp26_NTR_max_npp_agr), size=0.2, width=2, color="green",linetype="dotdash")+
  geom_line(aes(y=TNR26_100$npp_agr,color="TNR_26",linetype="TNR_26"))+
  geom_errorbar(aes(x=2100, ymin=rcp26_TNR_min_npp_agr, ymax=rcp26_TNR_max_npp_agr), size=0.2, width=2, color="red4",linetype="dotdash")+
  geom_line(aes(y=NTNR26_100$npp_agr,color="NTNR_26",linetype="NTNR_26"))+
  geom_errorbar(aes(x=2102, ymin=rcp26_NTNR_min_npp_agr, ymax=rcp26_NTNR_max_npp_agr), size=0.2, width=2, color="red",linetype="dotdash")+
  geom_line(aes(y=TR85_100$npp_agr,color="TR_85",linetype="TR_85"))+
  geom_errorbar(aes(x=2100, ymin=rcp85_TR_min_npp_agr, ymax=rcp85_TR_max_npp_agr), size=0.2, width=2, color="green4")+
  geom_line(aes(y=NTR85_100$npp_agr,color="NTR_85",linetype="NTR_85"))+
  geom_errorbar(aes(x=2102, ymin=rcp85_NTR_min_npp_agr, ymax=rcp85_NTR_max_npp_agr), size=0.2, width=2, color="green")+
  geom_line(aes(y=TNR85_100$npp_agr,color="TNR_85",linetype="TNR_85"))+
  geom_errorbar(aes(x=2100, ymin=rcp85_TNR_min_npp_agr, ymax=rcp85_TNR_max_npp_agr), size=0.2, width=2, color="red4")+
  geom_line(aes(y=NTNR85_100$npp_agr,color="NTNR_85",linetype="NTNR_85"))+
  geom_errorbar(aes(x=2102, ymin=rcp85_NTNR_min_npp_agr, ymax=rcp85_NTNR_max_npp_agr), size=0.2, width=2, color="red")+
  geom_line(aes(y=CLU26_100$npp_agr,color="CLU_26",linetype="CLU_26"))+
  geom_errorbar(aes(x=2104, ymin=rcp26_CLU_min_npp_agr, ymax=rcp26_CLU_max_npp_agr), size=0.2, width=2, color="blue",linetype="dotdash")+
  geom_line(aes(y=CLU85_100$npp_agr,color="CLU_85",linetype="CLU_85"))+
  geom_errorbar(aes(x=2104, ymin=rcp85_CLU_min_npp_agr, ymax=rcp85_CLU_max_npp_agr), size=0.2, width=2, color="blue",linetype="solid")+
  scale_color_manual(name=LegendTitle,values=c("TR_26"="green4","NTR_26"="green","TNR_26"="red4","NTNR_26"="red",
                                               "TR_85"="green4","NTR_85"="green","TNR_85"="red4","NTNR_85"="red",
                                               "CLU_26"="blue","CLU_85"="blue"))+
  scale_linetype_manual(name=LegendTitle,values=c("TR_26"="dotdash","NTR_26"="dotdash","TNR_26"="dotdash","NTNR_26"="dotdash",
                                                  "TR_85"="solid","NTR_85"="solid","TNR_85"="solid","NTNR_85"="solid",
                                                  "CLU_26"="dotted","CLU_85"="solid"))+
  ylab("Cropland NPP (PgC "~a^-1~")")+
  xlab("Year")+
  theme_bw()+
  theme(legend.position="none",
        axis.title.y=element_text(size=9),
        axis.title.x=element_text(size=9),
        plot.title=element_text(size=9))


b <-ggplot(TR26_100,aes(x=X,y=allc_agr)) +
  ggtitle("(B)")+
  geom_line(aes(y=TR26_100$allc_agr,color="TR_26",linetype="TR_26"))+
  geom_errorbar(aes(x=2100, ymin=rcp26_TR_min_allc_agr, ymax=rcp26_TR_max_allc_agr), size=0.2, width=2, color="green4",linetype="dotdash")+
  
  geom_line(aes(y=NTR26_100$allc_agr,color="NTR_26",linetype="NTR_26"))+
  geom_errorbar(aes(x=2102, ymin=rcp26_NTR_min_allc_agr, ymax=rcp26_NTR_max_allc_agr), size=0.2, width=2, color="green",linetype="dotdash")+
  
  geom_line(aes(y=TNR26_100$allc_agr,color="TNR_26",linetype="TNR_26"))+
  geom_errorbar(aes(x=2100, ymin=rcp26_TNR_min_allc_agr, ymax=rcp26_TNR_max_allc_agr), size=0.2, width=2, color="red4",linetype="dotdash")+
  
  geom_line(aes(y=NTNR26_100$allc_agr,color="NTNR_26",linetype="NTNR_26"))+
  geom_errorbar(aes(x=2102, ymin=rcp26_NTNR_min_allc_agr, ymax=rcp26_NTNR_max_allc_agr), size=0.2, width=2, color="red",linetype="dotdash")+
  
  geom_line(aes(y=TR85_100$allc_agr,color="TR_85",linetype="TR_85"))+
  geom_errorbar(aes(x=2101, ymin=rcp85_TR_min_allc_agr, ymax=rcp85_TR_max_allc_agr), size=0.2, width=2, color="green4")+
  
  geom_line(aes(y=NTR85_100$allc_agr,color="NTR_85",linetype="NTR_85"))+
  geom_errorbar(aes(x=2103, ymin=rcp85_NTR_min_allc_agr, ymax=rcp85_NTR_max_allc_agr), size=0.2, width=2, color="green")+
  
  geom_line(aes(y=TNR85_100$allc_agr,color="TNR_85",linetype="TNR_85"))+
  geom_errorbar(aes(x=2101, ymin=rcp85_TNR_min_allc_agr, ymax=rcp85_TNR_max_allc_agr), size=0.2, width=2, color="red4")+
  
  geom_line(aes(y=NTNR85_100$allc_agr,color="NTNR_85",linetype="NTNR_85"))+
  geom_errorbar(aes(x=2103, ymin=rcp85_NTNR_min_allc_agr, ymax=rcp85_NTNR_max_allc_agr), size=0.2, width=2, color="red")+
  
  geom_line(aes(y=CLU26_100$allc_agr,color="CLU_26",linetype="CLU_26"))+
  geom_errorbar(aes(x=2100, ymin=rcp26_CLU_min_allc_agr, ymax=rcp26_CLU_max_allc_agr), size=0.2, width=2, color="blue",linetype="dotdash")+
  geom_line(aes(y=CLU85_100$allc_agr,color="CLU_85",linetype="CLU_85"))+
  geom_errorbar(aes(x=2102, ymin=rcp85_CLU_min_allc_agr, ymax=rcp85_CLU_max_allc_agr), size=0.2, width=2, color="blue",linetype="solid")+
  scale_color_manual(name=LegendTitle,values=c("TR_26"="green4","NTR_26"="green","TNR_26"="red4","NTNR_26"="red",
                                               "TR_85"="green4","NTR_85"="green","TNR_85"="red4","NTNR_85"="red",
                                               "CLU_26"="blue","CLU_85"="blue"))+
  scale_linetype_manual(name=LegendTitle,values=c("TR_26"="dotdash","NTR_26"="dotdash","TNR_26"="dotdash","NTNR_26"="dotdash",
                                                  "TR_85"="solid","NTR_85"="solid","TNR_85"="solid","NTNR_85"="solid",
                                                  "CLU_26"="dotted","CLU_85"="solid"))+ 
  
  ylab("Cropland SOC (PgC)")+
  xlab("Year")+
  theme_bw()+
  theme(legend.position="none",
        axis.title.y=element_text(size=9),
        axis.title.x=element_text(size=9),
        plot.title=element_text(size=9))

c <-ggplot(TR26_100,aes(x=X,y=mrt_agr)) +
  ggtitle("(C)")+
  #geom_text(x=2000, y=max(NTR85_100$mrt_agr), label="(C)",size=4)+
  geom_line(aes(y=TR26_100$mrt_agr,colour="TR_26",linetype="TR_26"))+
  geom_errorbar(aes(x=2100, ymin=rcp26_TR_min_mrt_agr, ymax=rcp26_TR_max_mrt_agr), size=0.2, width=2, color="green4",linetype="dotdash")+
  
  geom_line(aes(y=NTR26_100$mrt_agr,colour="NTR_26",linetype="NTR_26"))+
  geom_errorbar(aes(x=2102, ymin=rcp26_NTR_min_mrt_agr, ymax=rcp26_NTR_max_mrt_agr), size=0.2, width=2, color="green",linetype="dotdash")+
  
  geom_line(aes(y=TNR26_100$mrt_agr,colour="TNR_26",linetype="TNR_26"))+
  geom_errorbar(aes(x=2100, ymin=rcp26_TNR_min_mrt_agr, ymax=rcp26_TNR_max_mrt_agr), size=0.2, width=2, color="red4",linetype="dotdash")+
  
  geom_line(aes(y=NTNR26_100$mrt_agr,colour="NTNR_26",linetype="NTNR_26"))+
  geom_errorbar(aes(x=2102, ymin=rcp26_NTNR_min_mrt_agr, ymax=rcp26_NTNR_max_mrt_agr), size=0.2, width=2, color="red",linetype="dotdash")+
  
  geom_line(aes(y=TR85_100$mrt_agr,colour="TR_85",linetype="TR_85"))+
  geom_errorbar(aes(x=2100, ymin=rcp85_TR_min_mrt_agr, ymax=rcp85_TR_max_mrt_agr), size=0.2, width=2, color="green4")+
  
  geom_line(aes(y=NTR85_100$mrt_agr,colour="NTR_85",linetype="NTR_85"))+
  geom_errorbar(aes(x=2102, ymin=rcp85_NTR_min_mrt_agr, ymax=rcp85_NTR_max_mrt_agr), size=0.2, width=2, color="green")+
  
  geom_line(aes(y=TNR85_100$mrt_agr,colour="TNR_85",linetype="TNR_85"))+
  geom_errorbar(aes(x=2100, ymin=rcp85_TNR_min_mrt_agr, ymax=rcp85_TNR_max_mrt_agr), size=0.2, width=2, color="red4")+
  
  geom_line(aes(y=NTNR85_100$mrt_agr,colour="NTNR_85",linetype="NTNR_85"))+
  geom_errorbar(aes(x=2102, ymin=rcp85_NTNR_min_mrt_agr, ymax=rcp85_NTNR_max_mrt_agr), size=0.2, width=2, color="red")+
  
  geom_line(aes(y=CLU26_100$mrt_agr,colour="CLU_26",linetype="CLU_26"))+
  geom_errorbar(aes(x=2100, ymin=rcp26_CLU_min_mrt_agr, ymax=rcp26_CLU_max_mrt_agr), size=0.2, width=2, color="blue",linetype="dotdash")+
  
  geom_line(aes(y=CLU85_100$mrt_agr,colour="CLU_85",linetype="CLU_85"))+
  geom_errorbar(aes(x=2100, ymin=rcp85_CLU_min_mrt_agr, ymax=rcp85_CLU_max_mrt_agr), size=0.2, width=2, color="blue",linetype="solid")+
  
  scale_color_manual(name=LegendTitle,values=c("TR_26"="green4","NTR_26"="green","TNR_26"="red4","NTNR_26"="red",
                                               "TR_85"="green4","NTR_85"="green","TNR_85"="red4","NTNR_85"="red",
                                               "CLU_26"="blue","CLU_85"="blue"))+
  scale_linetype_manual(name=LegendTitle,values=c("TR_26"="dotdash","NTR_26"="dotdash","TNR_26"="dotdash","NTNR_26"="dotdash",
                                                  "TR_85"="solid","NTR_85"="solid","TNR_85"="solid","NTNR_85"="solid",
                                                  "CLU_26"="dotted","CLU_85"="solid"))+ 
  ylab("Cropland turnover rate (%)")+
  xlab("Year")+
  theme_bw()+
  theme(legend.position="none",
        axis.title.y=element_text(size=9),
        axis.title.x=element_text(size=9),
        plot.title=element_text(size=9))


d <-ggplot(TR26_100,aes(x=X,y=litf_agr))+
  ggtitle("(D)")+
  geom_line(aes(y=TR26_100$litf_agr,colour="TR_26",linetype="TR_26"))+
  geom_errorbar(aes(x=2100, ymin=rcp26_TR_min_litf_agr, ymax=rcp26_TR_max_litf_agr), size=0.2, width=2, color="green4",linetype="dotdash")+
  
  geom_line(aes(y=NTR26_100$litf_agr,colour="NTR_26",linetype="NTR_26"))+
  geom_errorbar(aes(x=2102, ymin=rcp26_NTR_min_litf_agr, ymax=rcp26_NTR_max_litf_agr), size=0.2, width=2, color="green",linetype="dotdash")+
  
  geom_line(aes(y=TNR26_100$litf_agr,colour="TNR_26",linetype="TNR_26"))+
  geom_errorbar(aes(x=2100, ymin=rcp26_TNR_min_litf_agr, ymax=rcp26_TNR_max_litf_agr), size=0.2, width=2, color="red4",linetype="dotdash")+
  
  geom_line(aes(y=NTNR26_100$litf_agr,colour="NTNR_26",linetype="NTNR_26"))+
  geom_errorbar(aes(x=2102, ymin=rcp26_NTNR_min_litf_agr, ymax=rcp26_NTNR_max_litf_agr), size=0.2, width=2, color="red",linetype="dotdash")+
  
  geom_line(aes(y=TR85_100$litf_agr,colour="TR_85",linetype="TR_85"))+
  geom_errorbar(aes(x=2100, ymin=rcp85_TR_min_litf_agr, ymax=rcp85_TR_max_litf_agr), size=0.2, width=2, color="green4")+
  
  geom_line(aes(y=NTR85_100$litf_agr,colour="NTR_85",linetype="NTR_85"))+
  geom_errorbar(aes(x=2102, ymin=rcp85_NTR_min_litf_agr, ymax=rcp85_NTR_max_litf_agr), size=0.2, width=2, color="green")+
  
  geom_line(aes(y=TNR85_100$litf_agr,colour="TNR_85",linetype="TNR_85"))+
  geom_errorbar(aes(x=2100, ymin=rcp85_TNR_min_litf_agr, ymax=rcp85_TNR_max_litf_agr), size=0.2, width=2, color="red4")+
  
  geom_line(aes(y=NTNR85_100$litf_agr,colour="NTNR_85",linetype="NTNR_85"))+
  geom_errorbar(aes(x=2102, ymin=rcp85_NTNR_min_litf_agr, ymax=rcp85_NTNR_max_litf_agr), size=0.2, width=2, color="red")+
  
  geom_line(aes(y=CLU26_100$litf_agr,colour="CLU_26",linetype="CLU_26"))+
  geom_errorbar(aes(x=2100, ymin=rcp26_CLU_min_litf_agr, ymax=rcp26_CLU_max_litf_agr), size=0.2, width=2, color="blue",linetype="dotdash")+
  
  geom_line(aes(y=CLU85_100$litf_agr,colour="CLU_85",linetype="CLU_85"))+
  geom_errorbar(aes(x=2100, ymin=rcp85_CLU_min_litf_agr, ymax=rcp85_CLU_max_litf_agr), size=0.2, width=2, color="blue",linetype="solid")+
  
  scale_color_manual(name=LegendTitle,values=c("TR_26"="green4","NTR_26"="green","TNR_26"="red4","NTNR_26"="red",
                                               "TR_85"="green4","NTR_85"="green","TNR_85"="red4","NTNR_85"="red",
                                               "CLU_26"="blue","CLU_85"="blue"))+
  scale_linetype_manual(name=LegendTitle,values=c("TR_26"="dotdash","NTR_26"="dotdash","TNR_26"="dotdash","NTNR_26"="dotdash",
                                                  "TR_85"="solid","NTR_85"="solid","TNR_85"="solid","NTNR_85"="solid",
                                                  "CLU_26"="dotted","CLU_85"="solid"))+
  ylab("Cropland litterfall (PgC "~a^-1~")")+
  xlab("Year")+
  theme_bw()+
  theme(legend.position="none",
        axis.title.y=element_text(size=9),
        axis.title.x=element_text(size=9),
        plot.title=element_text(size=9))
grid.arrange(a,b,c,d)
dev.off()

#plot for legend in pdf:
pdf(paste(plots_path,"ggplot_global_values_RCPs_CONST_LU_LEGEND.pdf",sep=""), width=8, height=6)
LegendTitle<-"Scenario"
ggplot(TR26_100,aes(x=X,y=litf_agr)) +
  geom_text(x=2000, y=max(NTR85_100$litf_agr), label="(D)",size=4)+
  geom_line(aes(y=TR26_100$litf_agr,colour="T_R_26",linetype="T_R_26"))+
  geom_errorbar(aes(x=2100, ymin=rcp26_TR_min_litf_agr, ymax=rcp26_TR_max_litf_agr), size=0.2, width=2, color="green4",linetype="dotdash")+
  
  geom_line(aes(y=NTR26_100$litf_agr,colour="NT_R_26",linetype="NT_R_26"))+
  geom_errorbar(aes(x=2102, ymin=rcp26_NTR_min_litf_agr, ymax=rcp26_NTR_max_litf_agr), size=0.2, width=2, color="green",linetype="dotdash")+
  
  geom_line(aes(y=TNR26_100$litf_agr,colour="T_NR_26",linetype="T_NR_26"))+
  geom_errorbar(aes(x=2100, ymin=rcp26_TNR_min_litf_agr, ymax=rcp26_TNR_max_litf_agr), size=0.2, width=2, color="red4",linetype="dotdash")+
  
  geom_line(aes(y=NTNR26_100$litf_agr,colour="NT_NR_26",linetype="NT_NR_26"))+
  geom_errorbar(aes(x=2102, ymin=rcp26_NTNR_min_litf_agr, ymax=rcp26_NTNR_max_litf_agr), size=0.2, width=2, color="red",linetype="dotdash")+
  
  geom_line(aes(y=TR85_100$litf_agr,colour="T_R_85",linetype="T_R_85"))+
  geom_errorbar(aes(x=2100, ymin=rcp85_TR_min_litf_agr, ymax=rcp85_TR_max_litf_agr), size=0.2, width=2, color="green4")+
  
  geom_line(aes(y=NTR85_100$litf_agr,colour="NT_R_85",linetype="NT_R_85"))+
  geom_errorbar(aes(x=2102, ymin=rcp85_NTR_min_litf_agr, ymax=rcp85_NTR_max_litf_agr), size=0.2, width=2, color="green")+
  
  geom_line(aes(y=TNR85_100$litf_agr,colour="T_NR_85",linetype="T_NR_85"))+
  geom_errorbar(aes(x=2100, ymin=rcp85_TNR_min_litf_agr, ymax=rcp85_TNR_max_litf_agr), size=0.2, width=2, color="red4")+
  
  geom_line(aes(y=NTNR85_100$litf_agr,colour="NT_NR_85",linetype="NT_NR_85"))+
  geom_errorbar(aes(x=2102, ymin=rcp85_NTNR_min_litf_agr, ymax=rcp85_NTNR_max_litf_agr), size=0.2, width=2, color="red")+
  
  geom_line(aes(y=CLU26_100$litf_agr,colour="TRc05_26",linetype="TRc05_26"))+
  geom_errorbar(aes(x=2100, ymin=rcp26_CLU_min_litf_agr, ymax=rcp26_CLU_max_litf_agr), size=0.2, width=2, color="blue",linetype="dotdash")+
  
  geom_line(aes(y=CLU85_100$litf_agr,colour="TRc05_85",linetype="TRc05_85"))+
  geom_errorbar(aes(x=2100, ymin=rcp85_CLU_min_litf_agr, ymax=rcp85_CLU_max_litf_agr), size=0.2, width=2, color="blue",linetype="solid")+
  
  scale_color_manual(name=LegendTitle,values=c("T_R_26"="green4","NT_R_26"="green","T_NR_26"="red4","NT_NR_26"="red",
                                               "T_R_85"="green4","NT_R_85"="green","T_NR_85"="red4","NT_NR_85"="red",
                                               "TRc05_26"="blue","TRc05_85"="blue"))+
  scale_linetype_manual(name=LegendTitle,values=c("T_R_26"="dotdash","NT_R_26"="dotdash","T_NR_26"="dotdash","NT_NR_26"="dotdash",
                                                  "T_R_85"="solid","NT_R_85"="solid","T_NR_85"="solid","NT_NR_85"="solid",
                                                  "TRc05_26"="dotted","TRc05_85"="solid"))+
  ylab("litfall agriculture (PgC)")+
  xlab("")+
  theme_bw()+
  theme(legend.title=element_text(color="white",size=5), 
        legend.background = element_rect(fill = "white"),
        legend.key.size = unit(.5, "cm"), 
        legend.key.width = unit(.5, "cm"),
        legend.position="bottom")
dev.off()