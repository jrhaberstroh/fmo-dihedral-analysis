#!/bin/bash
set -o errexit
set -o nounset
SRCDIR=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )

cdctraj="$HOME/Code/photosynth/cdc-fmo/cdctraj.sh"
grodir="/home/jhaberstroh/data/2016-06-fmo500ns/2016-06-fmo500ns/output"
top="/mnt/pup1/jhaberstroh/data/2015-09-FMO_conf/4BCL_pp.top"

date=`date -I`

for i in {1..1}; do
    TESTFILE=$SRCDIR/output/_TEST_file.txt
    PYARGS="-residuenames" ATOMS=99548 TRJLEN=1 TOP=$top $cdctraj $grodir/ALL_md500-SEP-${i}.gro > $TESTFILE
done

