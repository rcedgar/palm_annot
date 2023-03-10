## palm_annot
Scripts, HMMs and search databases for identifying and classifying viral RdRp sequences

### Installation

This version **requires Ubuntu** (or compatible) due to use of precompiled binaries.

1. Clone the repo, e.g. to `~/palm_annot`.

2. Add `~/palm_annot/bin` and ``~/palm_annot/py` to your `$PATH`.

### palm_nuc_search.py

Search nucleotide sequences, e.g contigs or genomes, for RdRp and RdRp-like sequences. Output is the subset of input
sequences predicted to have RdRp or RdRp-like hits, no trimming or annotation is performed. Trimming and annotation of the
matching sequences can by performed by (1) translating to aa by 6-frame translation or ORF-finding, then (2) running
palm_annot.py on the aa sequences.
<pre>
palm_nuc_search.py --input contigs.fna --output hits.fna

Optional arguments:
  --tmpdir TMPDIR       Directory for temporary files
  --evalue EVALUE       Max E-value for HMM and diamond search (default 1e-3)
  --dbsize DBSIZE       Effective db size for E-value (-Z option of hmmsearch, default 100000)
  --threads THREADS     Number of threads for HMM and diamond search (default relevant options not set)
  --sensitive {fast,midsensitive,more-sensitive,very-sensitive}
                        diamond sensitivity option
</pre>

### palm_annot.py

Classify amino acid sequences as RdRp or non-RdRp palm domain and trim to the domain by deleting non-palm flanking sequence
using '150pp150' trimming, i.e. allow no more than 150aa before the palmprint start and no more than 150aa after palmprint
end. See palm_nuc_search.py if you have nucleotide sequence such as contigs or genomes.

<pre>
palm_annot.py --input sixframe.faa --fev hits.fev --rdrp rdrp.faa -xdxp xdxp.fa

Options:
  --fev FEV             Annotation output file (tab-separated text in field=value format)
  --fasta FASTA         FASTA output trimmed sequences (150pp150)
  --rdrp RDRP           FASTA output RdRp sequences (150pp150)
  --xdxp XDXP           FASTA output non-RdRp palm domain sequences (150pp150)
  --maxscorexdxp MAXSCOREXDXP
                        Maximum RdRp score for --xdxp FASTA output (0 to 100, default 25)
  --minscorerdrp MINSCORERDRP
                        Minimum score for --rdrp FASTA output (0 to 100, default 75)
  --threads THREADS     Number of threads (default depends on invoked script or binary)
  --tmpdir TMPDIR       Directory for temporary files (default /tmp)
</pre>
</pre>
