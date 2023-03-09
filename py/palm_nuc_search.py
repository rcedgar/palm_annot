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
" Trimming and annotation of the matching sequences can by performed by"
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

DomTbl = TmpPrefix + ".domtbl"
HMMDb = RepoDir + "hmmdbs/palm"

CmdLine = "hmmsearch"
CmdLine += " --domtbl " + DomTbl
if not Args.threads is None:
	CmdLine += " --cpu %d" % Args.threads
CmdLine += " -E %.3g" % Args.evalue
CmdLine += " -Z %d" % Args.dbsize
CmdLine += "  " + HMMDb
CmdLine += "  " + Args.input
CmdLine += " > /dev/null"
