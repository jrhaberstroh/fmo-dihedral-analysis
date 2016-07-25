#!/bin/bash
set -o nounset
set -o errexit
SRCDIR=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )
cd $SRCDIR
PC_OUT=$SRCDIR/output/gap_pc.txt


if [ -e $PC_OUT ]; then
    rm $PC_OUT
fi
for f in {1..7}; do
    file=$SRCDIR/output/cdc_c${f}_2016-07-14.txt
    python 07-14-b4gapRMSF.py -f $file -c $f -o $SUBGROUP_BASE/2016-07-15/c${f}_gapRMSF >> $PC_OUT
done

python 07-14-b4gap_pctable.py -f $PC_OUT
