#!/bin/bash -e

mkdir -p ../hmmdbs
cd ../hmmdbs

cat ../rdrp_plus_hmms/*.hmm \
  > rdrp_plus

cat ../rdrp_minus_hmms/*.hmm \
  > rdrp_minus

cat ../motif_hmms/*.hmm \
  > motifs

hmmpress rdrp_plus
hmmpress rdrp_minus
hmmpress motifs
