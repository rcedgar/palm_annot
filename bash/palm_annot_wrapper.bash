#!/bin/bash -e

if [ x$PALM_ANNOT_DIR == x ] ; then
	echo PALM_ANNOT_DIR not set
	exit 1
fi

if [ x$1 == x ] ; then
	echo Missing arg
	exit 1
fi

infa=$1

if [ ! -s $infa ] ; then
	echo Not found infa=$infa
	exit 1
fi

threads=`grep -c processor /proc/cpuinfo`

palmscan2 \
  -search_pssms $infa \
  -model $PALM_ANNOT_DIR/pssms/palm.ppm \
  -fev palmscan2_pssms.fev

palm_hmm_motif_search.py \
  --input $infa \
  --palmcore palmcore.fa \
  --fev hmm_motif_search.fev \
  --threads $threads

palm_hmm_search.py \
  --input $infa \
  --fev hmm_search.fev \
  --threads $threads

palm_diamond_motif_search.py \
  --input $infa \
  --fev diamond_motif_search.fev \
  --threads $threads
