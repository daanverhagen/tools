#!/bin/bash

set -x # echo on

NANOPOLISH=$1
RAWDATA=$2
BASECALLS=$3
REFGENOME=$4
MAPPED_READS=$5
NANOPOLISH_OUTPUT=$6
TOOLS=$7
READ_STATS=$8
SITE_STATS=$9

python $TOOLS/calc_site_stats.py -c 2.5 -i $NANOPOLISH_OUTPUT > $SITE_STATS