#!/bin/bash

gcms=("UKESM1-0-LL" "MRI-ESM2-0" "MPI-ESM1-2-HR" "IPSL-CM6A-LR" "GFDL-ESM4")
gc=("GCM1" "GCM2" "GCM3" "GCM4" "GCM5")
rcps=("ssp126" "ssp370" "ssp585" "ssp119" "ssp245" "ssp460")
Rs=("SSP126" "SSP370" "SSP585" "SSP119" "SSP245" "SSP460")

outdir="/p/projects/landuse/users/cmueller/ggcmi_phase3_nchecks_9ca735cb"

# define -DOBS for GSWP3-W5E5
# define -DCRU for CRU4
# define -DHIST for GCM specific historic runs (use these for spinup runs)
# define -DFUT for GCM-specific future runs, these need to define the RCP

mkdir -p $outdir/GSWP3-W5E5/historical/
mkdir -p $outdir/GSWP3-W5E5/historical_pnv/
mkdir -p $outdir/GSWP3-W5E5/historical_gsadapt/
#lpjsubmit -class standby -group macmit -o ${outdir}/GSWP3-W5E5/lpjml_log.out -e ${outdir}/GSWP3-W5E5/lpjml_log.err -blocking 16 256 -DOBS -DOUTDIR=$outdir/GSWP3-W5E5/historical lpjml_magpie.js
lpjsubmit -group macmit -o ${outdir}/GSWP3-W5E5/historical/lpjml_log.out -e ${outdir}/GSWP3-W5E5/historical/lpjml_log.err -blocking 16 512 -DOBS -DOUTDIR=$outdir/GSWP3-W5E5/historical -DFROM_RESTART lpjml_magpie.js
#lpjsubmit -class standby -group macmit -o ${outdir}/GSWP3-W5E5/historical_gsadapt/lpjml_log.out -e ${outdir}/GSWP3-W5E5/historical_gsadapt/lpjml_log.err -blocking 16 128 -DOBS -DOUTDIR=$outdir/GSWP3-W5E5/historical_gsadapt -DGSADAPT -DFROM_RESTART lpjml_magpie.js
#lpjsubmit -class standby -group macmit -wtime 6:00:00 -o ${outdir}/GSWP3-W5E5/historical_pnv/lpjml_log_pnv.out -e ${outdir}/GSWP3-W5E5/historical_pnv/lpjml_log_pnv.err -blocking 16 128 -DOBS -DPNV -DOUTDIR=$outdir/GSWP3-W5E5/historical_pnv -DFROM_RESTART lpjml_magpie.js

#mkdir -p $outdir/CRU4/historical/
#lpjsubmit -class standby -group macmit -o ${outdir}/CRU4/lpjml_log.out -e ${outdir}/CRU4/lpjml_log.err -blocking 16 128 -DCRU -DOUTDIR=$outdir/CRU4/historical lpjml_magpie.js
#lpjsubmit -class standby -group macmit -o ${outdir}/CRU4/historical/lpjml_log.out -e ${outdir}/CRU4/historical/lpjml_log.err -blocking 16 128 -DCRU -DOUTDIR=$outdir/CRU4/historical -DFROM_RESTART lpjml_magpie.js


#for ((g=2;g<5;g+=1))
#do
  #mkdir -p $outdir/${gcms[g]}/historical_gsadapt/
  #lpjsubmit -class standby -group macmit -o ${outdir}/${gcms[g]}/lpjml_log.out -e ${outdir}/${gcms[g]}/lpjml_log.err -blocking 16 128 -D${gc[g]} -DHIST -DOUTDIR=$outdir/${gcms[g]}/historical lpjml_magpie.js
  #lpjsubmit -class standby -group macmit -o ${outdir}/${gcms[g]}/historical_gsadapt/lpjml_log.out -e ${outdir}/${gcms[g]}/historical_gsadapt/lpjml_log.err -blocking 16 160 -D${gc[g]} -DHIST -DOUTDIR=$outdir/${gcms[g]}/historical_gsadapt -DFROM_RESTART -DGSADAPT lpjml_magpie.js
  #for ((r=0;r<3;r+=1))
  #for ((r=5;r<6;r+=1))
  #do
    #mkdir -p $outdir/${gcms[g]}/${rcps[r]}_gsadapt/
    #lpjsubmit -group macmit -o ${outdir}/${gcms[g]}/${rcps[r]}_gsadapt/lpjml_log.out -e ${outdir}/${gcms[g]}/${rcps[r]}_gsadapt/lpjml_log.err -blocking 16 128 -D${gc[g]} -DFUT -D${Rs[r]} -DOUTDIR=$outdir/${gcms[g]}/${rcps[r]}_gsadapt -DFROM_RESTART -DGSADAPT lpjml_magpie2.js
  #done
#done

