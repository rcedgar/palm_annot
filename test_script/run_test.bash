#!/bin/bash -e

mkdir -p ../test_output
cd ../test_output

palm_nuc_search.py \
  --input ../test_data/test_sequences.fna \
  --output nuc_search_output.fna

palmscan2 \
  -fasta_xlat nuc_search_output.fna \
  -trunclabels \
  -sep _frame= \
  -fastaout nuc_search.sixframe.faa

palm_annot.py \
  --input nuc_search.sixframe.faa \
  --fev hits.fev \
  --fasta rdrp.faa
