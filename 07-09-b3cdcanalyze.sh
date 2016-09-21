#!/bin/bash
set -o errexit
set -o nounset
SRCDIR=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )

cdctraj="$HOME/Code/photosynth/cdc-fmo/cdctraj.sh"
grodir="/home/jhaberstroh/data/2016-06-fmo500ns/2016-06-fmo500ns/output"
top="/mnt/pup1/jhaberstroh/data/2015-09-FMO_conf/4BCL_pp.top"

date=`date -I`


## 07-09 version
## for i in {1..1000}; do
##     TEMP=$(mktemp)
##     TEMPFILE=${TEMP}-${i}.txt
##     ATOMS=99548 TRJLEN=1 TOP=$top $cdctraj $grodir/ALL_md500-SEP-${i}.gro > $TEMPFILE
##     for l in {1..7}; do
##         head -n $l $TEMPFILE | tail -n 1 >> $SRCDIR/output/cdc_c${l}_${date}.txt
##     done
##     rm $TEMP
##     rm $TEMPFILE
## done

TEMP=$(mktemp)
TEMPFILE=${TEMP}-names.txt
ATOMS=99548 TRJLEN=1 PYARGS="-residuenames" TOP=$top $cdctraj $grodir/job525_md500-SEP-1.gro > $TEMPFILE
for l in {1..7}; do
    head -n 1 $TEMPFILE > $SRCDIR/output/cdc_c${l}_${date}.txt
done
rm $TEMP
rm $TEMPFILE

for i in {1..1700}; do
    TEMP=$(mktemp)
    TEMPFILE=${TEMP}-${i}.txt
    ATOMS=99548 TRJLEN=1 TOP=$top $cdctraj $grodir/job525_md500-SEP-${i}.gro > $TEMPFILE
    for l in {1..7}; do
        head -n $l $TEMPFILE | tail -n 1 >> $SRCDIR/output/cdc_c${l}_${date}.txt
    done
    rm $TEMP
    rm $TEMPFILE
done
