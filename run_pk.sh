#!/bin/bash

### set the directroy name under whith Gadget output (e.g., run00001/) is contained
DATADIR=/home/hikage/work/recon/

### set realization number
nrun=4000

### program
runcode=run_pk

### compile 

make pk_recon
mpic++ -O2 -openmp -openmp_report0 ${runcode}.cpp -o ${runcode}

### output directory name
outdir=output/

### mkdir ${outdir} if ${outdir} doesn't exist
if [ ! -d ${outdir} ]; then
mkdir $outdir
fi

./$runcode $DATADIR $outdir $nrun
