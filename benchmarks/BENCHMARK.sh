#!/bin/bash -x

# Modify the below paths to R to test, one on each line
ALLR=(
  ~/Projects/R/R-Array-Hash/bin/R,R-Array-Hash
  ~/Projects/R/optimized/devel/bin/R,R-devel
)
DATASETS=(
  ~/Projects/Wikipedia/SKEW.1mil
  ~/Projects/Wikipedia/DISTINCT.500thou
  ~/Projects/Wikipedia/DISTINCT.1mil
)
RUNS=`seq 1 3`

HASHSIZES=`Rscript -e 'cat(2^(10:15),"\n")'`

OPTS='--vanilla --no-save'
RESULTS='results.csv'

rm -rf $RESULTS error.log out.log

for run in $RUNS; do
for hashsize in $HASHSIZES; do
for dataset in ${DATASETS[*]}; do
for progname in ${ALLR[*]}; do
    progR=`echo $progname | awk -F, '{print $1}'`
    shortname=`echo $progname | awk -F, '{print $2}'`
    R_GC_GROWFRAC=0.75 R_GC_GROWINCFRAC=0.8 \
    R_GC_NGROWINCFRAC=0.8 R_GC_VGROWINCFRAC=0.8 \
    $progR $OPTS --args $dataset $hashsize $run $RESULTS < RUN.R \
      1>>out.log 2>>error.log
    echo '-----------------------' >> error.log
done
done
done
done

rm -rf runs
mkdir runs

for run in $RUNS; do
for progname in ${ALLR[*]}; do
    progR=`echo $progname | awk -F, '{print $1}'`
    shortname=`echo $progname | awk -F, '{print $2}'`
    $progR --slave < jsonlite/jsonlite.R > runs/jsonlite.Rout.${shortname}.${run} 2>&1
    $progR --slave < MASS-Ex.R > runs/MASS-Ex.Rout.${shortname}.${run} 2>&1
    $progR --slave < R-benchmark-25.R > runs/R-benchmark-25.Rout.${shortname}.${run} 2>&1
    $progR --slave < simplesum.R > runs/simplesum.Rout.${shortname}.${run} 2>&1
done
done
