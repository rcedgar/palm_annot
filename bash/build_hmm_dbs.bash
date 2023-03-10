#!/bin/bash -e

rm -rf ../hmmdbs
mkdir -p ../hmmdbs
cd ../hmmdbs

cat ../rdrp_plus_hmms/*.hmm \
  > rdrp_plus

cat ../rdrp_minus_hmms/*.hmm \
  > rdrp_minus

cat ../motif_hmms/*.hmm \
  > motifs

cat \
  rdrp_plus \
  rdrp_minus \
  motifs \
  > palm

hmmpress rdrp_plus
hmmpress rdrp_minus
hmmpress motifs
hmmpress palm
