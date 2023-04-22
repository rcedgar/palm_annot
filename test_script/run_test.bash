#!/bin/bash -e

mkdir -p ../test_output
cd ../test_output

palm_annot.py \
  --input ../test_data/test_sequences.fna \
  --seqtype nt \
  --fev hits.fev \
  --rdrp test_sequences.palmcore.faa
