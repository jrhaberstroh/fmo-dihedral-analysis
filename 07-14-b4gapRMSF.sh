#!/bin/bash
set -o nounset
set -o errexit
SRCDIR=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )

for f in {1..7}; do
    file=$SRCDIR/output/cdc_c${f}_2016-07-14.txt
    python 07-14-b4gapRMSF.py -f $file -c $f -o $SUBGROUP_BASE/2016-07-15/c${f}_gapRMSF
done
