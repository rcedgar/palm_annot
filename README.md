## palm_annot
Scripts, HMMs and search databases for identifying and classifying viral RdRp sequences

### Installation

This version **requires Ubuntu** (or compatible) due to use of precompiled binaries.

1. Clone the repo, e.g. to `~/palm_annot`.

2. Add `~/palm_annot/bin` and `~/palm_annot/py` to your `$PATH`.

3. Make sure execute permission is set on `~/palm_annot/bin/*` and `~/palm_annot/py/*`.

### palm_annot.py

Find RdRp and trim non-RdRp flanking sequence using '150pp150' trimming, i.e.
allow no more than 150aa before the palmprint start and no more than 150aa after palmprint end, partial
palmprints are allowed. This is also called 'palmcore' trimming.

<pre>
usage: palm_annot.py [-h] --input INPUT --seqtype {nt,aa} [--fev FEV] [--rdrp RDRP] [--fullnt FULLNT]
                     [--minscore MINSCORE] [--threads THREADS] [--tmpdir TMPDIR] [--keeptmp {no,yes}]
                     [--white {truncate,replace}] [--whitestr WHITESTR] [--dupes {delete,relabel}]
                     [--pattern PATTERN] [--framestr FRAMESTR] [--minpssmscore MINPSSMSCORE]

optional arguments:
  -h, --help            show this help message and exit
  --input INPUT         Input FASTA
  --seqtype {nt,aa}     Sequence type (nt or aa)
  --fev FEV             Annotation output file (tab-separated text in field=value format)
  --rdrp RDRP           Trimmed RdRp aa sequences (FASTA)
  --fullnt FULLNT       Full-length nt input sequences where RdRp found (FASTA)
  --minscore MINSCORE   Minimum score for RdRp (0 to 100, default 75)
  --threads THREADS     Number of threads (default don't set thread options)
  --tmpdir TMPDIR       Directory for temporary files (default /tmp)
  --keeptmp {no,yes}    Keep tmp files for trouble-shooting (default no)
  --white {truncate,replace}
                        Remove label white space by truncating or replacing (default truncate)
  --whitestr WHITESTR   Replacement string for --white replace (default '_')
  --dupes {delete,relabel}
                        Eliminate duplicate labels by deleting sequence or relabeling (default delete)
  --pattern PATTERN     Pattern to append to duplicate label, must contain @ which is replaced by 1,2...
                        (default '/dupe@')
  --framestr FRAMESTR   Append this string to translated sequence labels followed by frame -3 .. +3 (default
                        _frame=)
  --minpssmscore MINPSSMSCORE
                        Minimum PSSM score to report hit, 20 is high confidence (default 10.0)
</pre>

### fev2tsv.py

Convert tab-separated text name=value (fev) format to tab-separated text with
values only and optional header with field names.

<pre>
  --input INPUT         Input file in fev format (required)
  --output OUTPUT       Output file in tsv format (default stdout)
  --fields FIELDS       Comma-separated list of field names
                          (default: all fields found in fev file)
  --nullvalue NULLVALUE
                        String to use if value not specified (default empty string)
  --header {no,yes}     Include tsv header with field names as first line, yes
                        or no (default yes)
</pre>
