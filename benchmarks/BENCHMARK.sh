#!/bin/bash -x

# Modify the below paths to R to test, one on each line
ALLR=(
  ~/Projects/R-Array-Hash/bin/R
  ~/Projects/devel/bin/R
  Revo64
)
DATASETS=(
  SKEW.1mil
  DISTINCT.500thou
  DISTINCT.1mil
)
RUNS=`seq 1 10`

HASHSIZES=`Revoscript -e 'cat(2^(10:19),"\n")'`


OPTS='--vanilla --no-save'
RESULTS='results.csv'

rm -f $RESULTS error.log out.log

for run in $RUNS; do
for hashsize in $HASHSIZES; do
for dataset in ${DATASETS[*]}; do
for progname in ${ALLR[*]}; do
    R_GC_GROWFRAC=0.75 R_GC_GROWINCFRAC=0.8 \
    R_GC_NGROWINCFRAC=0.8 R_GC_VGROWINCFRAC=0.8 \
    $progname $OPTS --args $dataset $hashsize $run $RESULTS < RUN.R \
      1>>out.log 2>>error.log
    echo '-----------------------' >> error.log
done
done
done
done
