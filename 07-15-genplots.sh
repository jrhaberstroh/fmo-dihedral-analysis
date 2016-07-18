#!/bin/bash
set -o nounset
set -o errexit
SRCDIR=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )

$SRCDIR/07-14-a4plotrmsf.py
$SRCDIR/07-14-ab4_crossRMSF.sh
$SRCDIR/07-14-b4gapRMSF.sh
$SRCDIR/07-14-b5gap_exclusion.py
