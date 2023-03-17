#!/usr/bin/env python3

import sys
import argparse
import os
import random
import fasta

Usage = \
(
"Find RdRp sequences in nt or aa sequences and trim non-RdRp "
" flanking sequence using '150pp150' trimming, i.e. allow no"
" more than 150aa before the palmprint start and no more than"
" 150aa after palmprint end, partial palmprints are allowed."
)

AP = argparse.ArgumentParser(description = Usage)

AP.add_argument("--input",
  required=True,
  help="Input FASTA")

AP.add_argument("--seqtype",
  required=True,
  choices=[ "nt", "aa" ],
  help="Sequence type (nt or aa)")

AP.add_argument("--fev",
  required=False,
  help="Annotation output file (tab-separated text in field=value format)")

AP.add_argument("--rdrp",
  required=False,
  help="Trimmed RdRp aa sequences (FASTA)")

AP.add_argument("--fullnt",
  required=False,
  help="Full-length nt input sequences where RdRp found (FASTA)")

AP.add_argument("--minscore",
  required=False,
  type=float,
  default=75,
  help="Minimum score for RdRp (0 to 100, default 75)")

AP.add_argument("--threads",
  type=int,
  required=False,
  help="Number of threads (default don't set thread options)")

AP.add_argument("--tmpdir",
  required=False,
  default="/tmp",
  help="Directory for temporary files (default /tmp)")

AP.add_argument("--keeptmp",
  required=False,
  default="no",
  choices=[ "no", "yes" ],
  help="Keep tmp files for trouble-shooting (default no)")

AP.add_argument("--white",
  required=False,
  choices=[ "truncate", "replace" ],
  default="truncate",
  help="Remove label white space by truncating or replacing (default truncate)")

AP.add_argument("--whitestr",
  required=False,
  default="_",
  help="Replacement string for --white replace (default '_')")

AP.add_argument("--dupes",
  required=False,
  choices=[ "delete", "relabel" ],
  default="delete",
  help="Eliminate duplicate labels by deleting sequence or relabeling (default delete)")

AP.add_argument("--pattern",
  required=False,
  default="/dupe@",
  help="Pattern to append to duplicate label, must contain @ which is replaced by 1,2... (default '/dupe@')")

AP.add_argument("--framestr",
  required=False,
  default="_frame=",
  help="Append this string to translated sequence labels followed by frame -3 .. +3 (default _frame=)")

AP.add_argument("--minpssmscore",
  required=False,
  type=float,
  default=10.0,
  help="Minimum PSSM score to report hit, 20 is high confidence (default 10.0)")

Args = AP.parse_args()

if not Args.fullnt is None and Args.seqtype != "nt":
	sys.stderr.write("\n\nERROR --fullnt not supported for aa input\n\n")
	sys.exit(1)

# Parent directory of this script is git repo, needed to find HMMs etc.
RepoDir = os.path.dirname(os.path.realpath(__file__ + "/..")) + "/"
assert os.path.isdir(RepoDir)

TmpDir = Args.tmpdir
if not TmpDir.endswith("/"):
	TmpDir += "/"
if not os.path.isdir(TmpDir):
	sys.stderr.write("Not a directory: --tmpdir %s\n" % TmpDir)
	assert False

pid = os.getpid()
r = random.randint(0, 999999)
TmpPrefix = TmpDir + "palm_annot.%d.%d." % (pid, r)

fLog = None
if Args.keeptmp == "yes":
	LogFN = TmpPrefix + "exec.log"
	fLog = open(LogFN, "w")

def Exec(CmdLine):
	if not fLog is None:
		fLog.write("exec %s\n" % CmdLine)
	Code = os.system(CmdLine)
	if not fLog is None:
		fLog.write("code %d\n" % Code)
	if Code != 0:
		sys.stderr.write("\n")
		sys.stderr.write(CmdLine + "\n")
		sys.stderr.write("\n")
		sys.stderr.write("Error code %d\n" % Code)
		assert False

TmpFiles = []
def MakeTmp(Name):
	global TmpFiles, TmpPrefix
	FN = TmpPrefix + Name
	TmpFiles.append(FN)
	return FN

CleanFN = MakeTmp("clean")

CmdLine = "fasta_clean.py"
CmdLine += " --input " + Args.input
CmdLine += " --output " + CleanFN
CmdLine += " --white " + Args.white
CmdLine += " --whitestr " + Args.whitestr
CmdLine += " --dupes " + Args.dupes
CmdLine += " --pattern " + Args.pattern
Exec(CmdLine)

if Args.seqtype == "nt":
	AaFN = MakeTmp("xlat")
	CmdLine = RepoDir + "bin/palmscan2"
	CmdLine += " -fasta_xlat " + CleanFN
	CmdLine += " -sep '%s'" % Args.framestr
	CmdLine += " -fastaout " + AaFN
	Exec(CmdLine)
	if Args.fullnt is None:
		Exec("rm -f " + CleanFN)
elif Args.seqtype == "aa":
	AaFN = CleanFN
else:
	assert False

PSSMModelFN = RepoDir + "pssms/palm.ppm"
PSSM_fev = MakeTmp("pssm.fev")
HMM_motif_fev = MakeTmp("hmm_motif.fev")
HMM_pm_fev = MakeTmp("hmm_pm.fev")
Dmnd_fev = MakeTmp("dmnd.fev")

CmdLine = RepoDir + "bin/palmscan2"
CmdLine += " -search_pssms " + AaFN
CmdLine += " -model " + PSSMModelFN
CmdLine += " -trunclabels"
CmdLine += " -min_palm_score %.1f" % Args.minpssmscore
CmdLine += " -fev " + PSSM_fev
if not Args.threads is None:
	CmdLine += "  -threads %d" % Args.threads
if Args.keeptmp == "yes":
	LogFN = MakeTmp("search_pssms.log")
	CmdLine += "  -log " + LogFN
Exec(CmdLine)

CmdLine = RepoDir + "py/palm_hmm_motif_search.py"
CmdLine += "  --input " + AaFN
CmdLine += " --fev " + HMM_motif_fev
if not Args.threads is None:
	CmdLine += "  --threads %d" % Args.threads
Exec(CmdLine)

CmdLine = RepoDir + "py/palm_hmm_search.py"
CmdLine += "  --input " + AaFN
CmdLine += "  --fev " + HMM_pm_fev
if not Args.threads is None:
	CmdLine += "  --threads %d" % Args.threads
Exec(CmdLine)

CmdLine = RepoDir + "py/palm_diamond_motif_search.py"
CmdLine += "  --input " + AaFN
CmdLine += " --fev " + Dmnd_fev
if not Args.threads is None:
	CmdLine += "  --threads %d" % Args.threads
Exec(CmdLine)

OutputLabels = set()
LabelToFields = {}
for FN in [ PSSM_fev, HMM_motif_fev, HMM_pm_fev, Dmnd_fev ]:
	Labels1 = set()
	for Line in open(FN):
		Fields = Line[:-1].split('\t')

		Label = Fields[0]
		if not Label in OutputLabels:
			LabelToFields[Label] = []
			OutputLabels.add(Label)
		for Field in Fields[1:]:
			LabelToFields[Label].append(Field)

Tmp_fev = MakeTmp("merged.fev")
fTmpFev = open(Tmp_fev, "w")
LabelToL = {}
OutCount = 0
SeqCount = 0
def OnSeq(Label, Seq):
	global fTmpFev
	global OutCount
	global SeqCount
	SeqCount += 1

	if Label in OutputLabels:
		L = len(Seq)
		Line = Label
		for Field in LabelToFields[Label]:
			Line += "\t" + Field
		Line += "\taaseq=" + Seq
		fTmpFev.write(Line + "\n")
		OutCount += 1

sys.stderr.write("Merging...")
fasta.ReadSeqsOnSeq(AaFN, OnSeq, True)
fTmpFev.close()
sys.stderr.write(" done.\n")

CmdLine = RepoDir + "bin/palmscan2"
CmdLine += " -pamerge " + Tmp_fev
CmdLine += " -minscore %.4g" % Args.minscore
if not Args.fev is None:
	CmdLine += " -fev " + Args.fev
if not Args.rdrp is None:
	CmdLine += " -fasta " + Args.rdrp
if Args.keeptmp == "no":
	LogFN = MakeTmp("pamerge.log")
	CmdLine += " -log " + LogFN
Exec(CmdLine)

def OnSeqFull(Label, Seq):
###	print("Label=" + Label)
	if Label in NtLabels:
		fasta.WriteSeq(fFull, Seq, Label)

if not Args.fullnt is None:
	if Args.seqtype != "nt":
		sys.stderr
	fFull = open(Args.fullnt, "w")
	NtLabels = set()
	for Label in OutputLabels:
		NtLabel = Label.split(Args.framestr)[0]
		NtLabels.add(NtLabel)
###	print(NtLabels)
	fasta.ReadSeqsOnSeq(CleanFN, OnSeqFull)
	fFull.close()

if Args.keeptmp == "no":
	for FN in TmpFiles:
		Exec("rm -f " + FN)

Pct = 0
if SeqCount > 0:
	Pct = (100.0*OutCount)/SeqCount

sys.stderr.write("\n%u of %u (%.4g%%) palm sequences reported\n\n" % (SeqCount, OutCount, Pct))
