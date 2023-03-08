#!/usr/bin/env python3

import sys
import argparse
import os
import random
import fasta

Usage = \
(
"Useage: TODO\n"
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

AP.add_argument("--fev",
  required=True,
  help="Output file (fev format)")

AP.add_argument("--fasta",
  required=True,
  help="Trimmed sequences (150pp150, fasta format)")

AP.add_argument("--rdrp",
  required=False,
  help="Trimmed RdRp sequences (150pp150, fasta format)")

AP.add_argument("--xdxp",
  required=False,
  help="Trimmed non-RdRp palm domain sequences (150pp150, fasta format)")

AP.add_argument("--maxscorexdxp",
  required=False,
  type=float,
  default=25,
  help="Maximum RdRp score for --xdxp FASTA output (default 25)")

AP.add_argument("--minscorerdrp",
  required=False,
  type=float,
  default=75,
  help="Minimum score for --rdrp FASTA output (default 75)")

AP.add_argument("--threads",
  type=int,
  required=False,
  help="Number of threads (default depends on invoked script or binary)")

AP.add_argument("--tmpdir",
  required=False,
  default="/tmp",
  help="Directory for temporary files (default /tmp)")

Args = AP.parse_args()

TmpDir = Args.tmpdir
if not TmpDir.endswith("/"):
	TmpDir += "/"
if not os.path.isdir(TmpDir):
	sys.stderr.write("Not a directory: --tmpdir %s\n" % TmpDir)
	assert False

pid = os.getpid()
r = random.randint(0, 999999)
TmpPrefix = TmpDir + "palm_annot.pa.%d.%d." % (pid, r)
sys.stderr.write("TmpPrefix = %s\n" % TmpPrefix)

def Exec(CmdLine):
	Code = os.system(CmdLine)
	if Code != 0:
		sys.stderr.write("\n")
		sys.stderr.write(CmdLine + "\n")
		sys.stderr.write("\n")
		sys.stderr.write("Error code %d\n" % Code)
		assert False

PSSMModelFN = RepoDir + "pssms/palm.ppm"
PSSM_fev = TmpPrefix + "pssm.fev"
HMM_motif_fev = TmpPrefix + "hmm_motif.fev"
HMM_pm_fev = TmpPrefix + "hmm_pm.fev"
Dmnd_fev = TmpPrefix + "dmnd.fev"

CmdLine = RepoDir + "bin/palmscan2"
CmdLine += " -search_pssms " + Args.input
CmdLine += " -model " + PSSMModelFN
CmdLine += " -trunclabels"
CmdLine += "  -fev " + PSSM_fev
if not Args.threads is None:
	CmdLine += "  -threads " + Args.threads
Exec(CmdLine)

CmdLine = RepoDir + "py/palm_hmm_motif_search.py"
CmdLine += "  --input " + Args.input
CmdLine += " --fev " + HMM_motif_fev
if not Args.threads is None:
	CmdLine += "  --threads " + Args.threads
Exec(CmdLine)

CmdLine = RepoDir + "py/palm_hmm_search.py"
CmdLine += "  --input " + Args.input
CmdLine += "  --fev " + HMM_pm_fev
if not Args.threads is None:
	CmdLine += "  --threads " + Args.threads
Exec(CmdLine)

CmdLine = RepoDir + "py/palm_diamond_motif_search.py"
CmdLine += "  --input " + Args.input
CmdLine += " --fev " + Dmnd_fev
if not Args.threads is None:
	CmdLine += "  --threads " + Args.threads
Exec(CmdLine)

Labels = set()
LabelToFields = {}
for FN in [ PSSM_fev, HMM_motif_fev, HMM_pm_fev, Dmnd_fev ]:
	for Line in open(FN):
		Fields = Line[:-1].split('\t')
		Label = Fields[0]
		if not Label in Labels:
			LabelToFields[Label] = []
			Labels.add(Label)
		for Field in Fields[1:]:
			LabelToFields[Label].append(Field)

Tmp_fev = TmpPrefix + "merged.fev"
fTmpFev = open(Tmp_fev, "w")
LabelToL = {}
OutCount = 0
SeqCount = 0
def OnSeq(Label, Seq):
	global fTmpFev
	global OutCount
	global SeqCount
	SeqCount += 1
	if Label in Labels:
		L = len(Seq)
		Line = Label
		for Field in LabelToFields[Label]:
			Line += "\t" + Field
		Line += "\taaseq=" + Seq
		fTmpFev.write(Line + "\n")
		OutCount += 1

sys.stderr.write("Merging...")
fasta.ReadSeqsOnSeq(Args.input, OnSeq, True)
fTmpFev.close()
sys.stderr.write(" done.\n")

CmdLine = RepoDir + "bin/palmscan2"
CmdLine += " -pamerge " + Tmp_fev
CmdLine += " -fev " + Args.fev
CmdLine += " -fasta " + Args.fasta
CmdLine += " -minscore %.4g" % Args.minscorerdrp
CmdLine += " -maxscore %.4g" % Args.maxscorexdxp
if not Args.rdrp is None:
	CmdLine += " -fasta_rdrp " + Args.rdrp
if not Args.rdrp is None:
	CmdLine += " -fasta_xdxp " + Args.xdxp
Exec(CmdLine)

## for FN in [ PSSM_fev, HMM_motif_fev, HMM_pm_fev, Dmnd_fev, Tmp_fev ]:
	## Exec("rm -f " + FN)

Pct = 0
if SeqCount > 0:
	Pct = (100.0*OutCount)/SeqCount

sys.stderr.write("\n%u of %u (%.4g%%) palm sequences reported\n\n" % (SeqCount, OutCount, Pct))
