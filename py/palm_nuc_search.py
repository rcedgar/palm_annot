#!/usr/bin/env python3

import sys
import argparse
import os
import random
import fasta

Usage = \
(
"Search nucleotide sequences, e.g contigs or genomes, for RdRp and RdRp-like sequences."
" Output is the subset of input sequences predicted to have RdRp or RdRp-like hits,"
" no trimming or annotation is performed."
" Trimming and annotation of the matching sequences can be done by"
" (1) translating to aa by 6-frame translation or ORF-finding, then"
" (2) running palm_annot.py on the aa sequences."
)

# RepoDir = path name of repository
# Can be set in environment variable PALM_ANNOT_DIR, if not
#   set is assumed to be parent directory of this script.
RepoDir = os.environ.get("PALM_ANNOT_DIR", None)
if RepoDir is None:
	RepoDir = os.path.dirname(os.path.realpath(__file__ + "/.."))
if not RepoDir.endswith("/"):
	RepoDir += "/"

AP = argparse.ArgumentParser(description = Usage)

AP.add_argument("--input",
  required=True,
  help="Input FASTA")

AP.add_argument("--output",
  required=True,
  help="Sequences with RdRp or RdRp-like hits (FASTA format)")

AP.add_argument("--tmpdir",
  required=False,
  default="/tmp",
  help="Directory for temporary files")

AP.add_argument("--evalue",
  required=False,
  default=1e-3,
  help="Max E-value for HMM and diamond search (default 1e-3)")

AP.add_argument("--dbsize",
  required=False,
  default=100000,
  help="Effective db size for E-value (-Z option of hmmsearch, default 100000)")

AP.add_argument("--threads",
  type=int,
  required=False,
  help="Number of threads for HMM and diamond search (default relevant options not set)")

AP.add_argument("--sensitive",
  required=False,
  default="very-sensitive",
  choices=[ "fast", "midsensitive", "more-sensitive", "very-sensitive"],
  help="diamond sensitivity option")

AP.add_argument("--keeptmp",
  required=False,
  default="no",
  choices=[ "no", "yes" ],
  help="Keep tmp files for trouble-shooting (default no)")


Args = AP.parse_args()
def Exec(CmdLine):
	Code = os.system(CmdLine)
	if Code != 0:
		sys.stderr.write("\n")
		sys.stderr.write(CmdLine + "\n")
		sys.stderr.write("\n")
		sys.stderr.write("Error code %d\n" % Code)
		assert False

TmpDir = Args.tmpdir
if not TmpDir.endswith("/"):
	TmpDir += "/"
if not os.path.isdir(TmpDir):
	sys.stderr.write("Not a directory: --tmpdir %s\n" % TmpDir)
	assert False

pid = os.getpid()
r = random.randint(0, 999999)
TmpPrefix = TmpDir + "palm_nuc_search.%d.%d." % (pid, r)
sys.stderr.write("TmpPrefix = %s\n" % TmpPrefix)

fOut = open(Args.output, "w")

TmpFa = TmpPrefix + "xlat.fa"

CmdLine = "palmscan2"
CmdLine += " -fasta_xlat " + Args.input
CmdLine += " -trunclabels"
CmdLine += " -sep '|frame='"
CmdLine += " -fastaout " + TmpFa

sys.stderr.write("Six-frame translation...")
sys.stderr.flush()
Exec(CmdLine)
sys.stderr.write(" done.\n")

HitsFN = TmpPrefix + "hits"
HMMDb = RepoDir + "hmmdbs/palm"

CmdLine = "hmmsearch"
CmdLine += " -o /dev/null"
CmdLine += " --domtbl /dev/stdout"
if not Args.threads is None:
	CmdLine += " --cpu %d" % Args.threads
CmdLine += " -E %.3g" % Args.evalue
CmdLine += " -Z %d" % Args.dbsize
CmdLine += "  " + HMMDb
CmdLine += "  " + TmpFa
CmdLine += " | grep -v ^#"
CmdLine += " | sed '-es/[| ].*//'"
CmdLine += " > " + HitsFN

sys.stderr.write("Running hmmsearch...")
sys.stderr.flush()
Exec(CmdLine)
sys.stderr.write(" done.\n")

RefDb = RepoDir + "diamond_refdbs/rdrp_plus_abc"

CmdLine = "diamond blastx"
CmdLine += " --query " + Args.input
CmdLine += " --evalue %.3g" % Args.evalue
CmdLine += "  --max-target-seqs 1"
CmdLine += " --%s" % Args.sensitive
CmdLine += " --db " + RefDb
if not Args.threads is None:
	CmdLine += " --threads %d" % Args.threads
CmdLine += " --outfmt 6 qseqid"
CmdLine += " >> " + HitsFN

sys.stderr.write("Running diamond ...")
sys.stderr.flush()
Exec(CmdLine)
sys.stderr.write(" done.\n")

sys.stderr.write("Unique hits ...")
sys.stderr.flush()

def StripLabel(Label):
	Label = Label.strip()
	Label = Label.split()[0]
	Label = Label.split('|')[0]
	return Label

Labels = set()
for Line in open(HitsFN):
	Label = StripLabel(Line)
	Labels.add(Label)
N = len(Labels)
sys.stderr.write(" done, %d hits.\n" % N)

Nout = 0
Label0Set = set()
def OnSeq(Label, Seq):
	global Nout
	Label0 = StripLabel(Label)
	if Label0Set in Label0Set:
		sys.stderr.write("\n===ERROR===\nDuplicate label >%s\n\n" % Label0)
		sys.exit(1)

	if Label0 in Labels:
		Label0Set.add(Label)
		fasta.WriteSeq(fOut, Seq, Label)

fasta.ReadSeqsOnSeq(Args.input, OnSeq)
sys.stderr.write("%d sequences output\n" % Nout)
fOut.close()

if Args.keeptmp == "no":
	CmdLine = "rm -f " + TmpFa + " " + HitsFN
	Exec(CmdLine)
