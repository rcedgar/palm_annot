#!/bin/bash -e

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

