#!/bin/bash


gcms=("UKESM1-0-LL" "MRI-ESM2-0" "MPI-ESM1-2-HR" "IPSL-CM6A-LR" "GFDL-ESM4")
gc=("GCM1" "GCM2" "GCM3" "GCM4" "GCM5")
rcps=("ssp126" "ssp370" "ssp585")
Rs=("SSP126" "SSP370" "SSP585")
till=("TILL" "NOTILL")
fert=("NLIM0" "NLIM1" "NLIM2" "NLIM3" "NLIM4")
lsu=("LSU1" "LSU2" "LSU3" "LSU4" "LSU5")

#outdir="/p/projects/landuse/users/cmueller/ggcmi_phase3_nchecks_fbed5c8b"
outdir="../output"

# define -DOBS for GSWP3-W5E5
# define -DCRU for CRU4
# define -DHIST for GCM specific historic runs (use these for spinup runs)
# define -DFUT for GCM-specific future runs, these need to define the RCP

#mkdir -p $outdir/CRU4/historical/
#lpjsubmit -class standby -group macmit -o ${outdir}/CRU4/lpjml_log.out -e ${outdir}/CRU4/lpjml_log.err -blocking 16 128 -DCRU -DOUTDIR=$outdir/CRU4/historical lpjml_magpie.js
#lpjsubmit -class standby -group macmit -o ${outdir}/CRU4/historical/lpjml_log.out -e ${outdir}/CRU4/historical/lpjml_log.err -blocking 16 128 -DCRU -DOUTDIR=$outdir/CRU4/historical -DFROM_RESTART lpjml_magpie.js


for ((g=2;g<3;g+=1)) # only GCM3
do
  #mkdir -p $outdir/${gcms[g]}/historical/
  #lpjsubmit -class standby -group macmit -o ${outdir}/${gcms[g]}/lpjml_log.out -e ${outdir}/${gcms[g]}/lpjml_log.err -blocking 16 128 -D${gc[g]} -DHIST -DOUTDIR=$outdir/${gcms[g]}/historical lpjml_climasteppe.js
  #lpjsubmit -class standby -group macmit -o ${outdir}/${gcms[g]}/historical/lpjml_log.out -e ${outdir}/${gcms[g]}/historical/lpjml_log.err -blocking 16 160 -D${gc[g]} -DHIST -DOUTDIR=$outdir/${gcms[g]}/historical -DFROM_RESTART lpjml_magpie.js
  echo "historical"
  #lpjcheck -D${gc[g]} -DHIST -DOUTDIR=$outdir/${gcms[g]}/historical -DFROM_RESTART lpjml_climasteppe.js
  for ((l=0;l<5;l+=1))
  do
    mkdir -p ${outdir}/${gcms[g]}/historical_${lsu[l]}
    #lpjcheck -D${gc[g]} -D${lsu[l]} -DHIST -DOUTDIR=$outdir/${gcms[g]}/historical_${lsu[l]} -DFROM_RESTART lpjml_climasteppe.js
    lpjsubmit -class standby -wtime 10:00:00 -group macmit -o ${outdir}/${gcms[g]}/historical/lpjml_log.out -e ${outdir}/${gcms[g]}/historical/lpjml_log.err -blocking 16 128 -D${gc[g]} -D${lsu[l]} -DHIST -DOUTDIR=$outdir/${gcms[g]}/historical_${lsu[l]} -DFROM_RESTART lpjml_climasteppe.js
  done
  #for ((t=0;t<2;t+=1))
  #do
  # for notill loops over all nlim and lsus
  echo "loop over notill scenarios"
  t=1
    for ((k=0;k<5;k+=1))
    do
      for ((l=0;l<5;l+=1))
      do
        echo $t $k $l ${till[$t]} ${fert[k]} ${lsu[l]}
        mkdir -p ${outdir}/${gcms[g]}/${till[t]}_${fert[k]}_${lsu[l]}
        #lpjcheck -D${gc[g]} -DHIST -DOUTDIR=$outdir/${gcms[g]}/${till[t]}_${fert[k]}_${lsu[l]} -DFROM_RESTART -D${till[t]} -D${fert[k]} -D${lsu[l]} lpjml_climasteppe.js
        lpjsubmit -class standby -wtime 10:00:00 -group macmit -o ${outdir}/${gcms[g]}/${till[t]}_${fert[k]}_${lsu[l]}/lpjml_log.out -e ${outdir}/${gcms[g]}/${till[t]}_${fert[k]}_${lsu[l]}/lpjml_log.err -blocking 16 128 -D${gc[g]} -DHIST -DOUTDIR=$outdir/${gcms[g]}/${till[t]}_${fert[k]}_${lsu[l]} -DFROM_RESTART -D${till[t]} -D${fert[k]} -D${lsu[l]} lpjml_climasteppe.js
      done
    done
  #done
  # for tillage loop over all nlims for lsu1
  echo "loop over till scenarios"
  t=0
    for ((k=0;k<5;k+=1))
    do
      #for ((l=0;l<5;l+=1))
      #do
      l=0
        echo $t $k $l ${till[$t]} ${fert[k]} ${lsu[l]}
        mkdir -p ${outdir}/${gcms[g]}/${till[t]}_${fert[k]}_${lsu[l]}
        #lpjcheck -D${gc[g]} -DHIST -DOUTDIR=$outdir/${gcms[g]}/${till[t]}_${fert[k]}_${lsu[l]} -DFROM_RESTART -D${till[t]} -D${fert[k]} -D${lsu[l]} lpjml_climasteppe.js
        lpjsubmit -class standby -wtime 10:00:00 -group macmit -o ${outdir}/${gcms[g]}/${till[t]}_${fert[k]}_${lsu[l]}/lpjml_log.out -e ${outdir}/${gcms[g]}/${till[t]}_${fert[k]}_${lsu[l]}/lpjml_log.err -blocking 16 128 -D${gc[g]} -DHIST -DOUTDIR=$outdir/${gcms[g]}/${till[t]}_${fert[k]}_${lsu[l]} -DFROM_RESTART -D${till[t]} -D${fert[k]} -D${lsu[l]} lpjml_climasteppe.js
      #done
    done
##  for ((r=0;r<3;r+=1))
##  do
##    mkdir -p $outdir/${gcms[g]}/${rcps[r]}/
##    lpjcheck -D${gc[g]} -DHIST -D${Rs[r]} -DOUTDIR=$outdir/${gcms[g]}/${rcps[r]} lpjml_climasteppe.js
##  done
done

