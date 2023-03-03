#!/usr/bin/env python3

import sys
import argparse
import os
import random
import fasta
import polhummer

Usage = \
(
"Search sequences for palmprint A,B and C motifs\n"
"using PH/PN HMMs (motif match states are annotated).\n"
)

if len(sys.argv) == 0:
	print(Usage)

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

AP.add_argument("--palmcore",
  required=False,
  help="FASTA output with RdRp/RT sequences trimmed to palmcore with motif annotations")

AP.add_argument("--fev",
  required=False,
  help="Results in field-equals-value (fev) format")

AP.add_argument("--tmpdir",
  required=False,
  default="/tmp",
  help="Directory for temporary files")

AP.add_argument("--evalue",
  required=False,
  default=1e-3,
  help="Max E-value for HMM search (default and recommended 1e-3)")

AP.add_argument("--dbsize",
  required=False,
  default=100000,
  help="Effective db size for E-value (-Z option of hmmsearch, default 100000)")

AP.add_argument("--threads",
  required=False,
  help="Number of hmmsearch threads (default --cpu option of hmmsearch not set)")

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

sys.stderr.write("Reading input...")
sys.stderr.flush()
SeqDict = fasta.ReadSeqsDict(Args.input)
SeqCount = len(SeqDict.keys())
sys.stderr.write(" %d seqs.\n" % SeqCount)

pid = os.getpid()
r = random.randint(0, 999999)
TmpPrefix = TmpDir + "phms%d.%d." % (pid, r)
sys.stderr.write("TmpPrefix = %s\n" % TmpPrefix)

fCore = None
fFev = None
if not Args.palmcore is None:
	fCore = open(Args.palmcore, "w")
if not Args.fev is None:
	fFev = open(Args.fev , "w")

DomTbl = TmpPrefix + "ph.domtbl"
HMMDb = RepoDir + "hmmdbs/motifs"

CmdLine = "hmmsearch"
CmdLine += " --domtbl " + DomTbl
CmdLine += " -E %.3g" % Args.evalue
CmdLine += " -Z %d" % Args.dbsize
CmdLine += "  " + HMMDb
CmdLine += "  " + Args.input
if not Args.threads is None:
	CmdLine += " --cpu %d" % Args.threads
CmdLine += " > /dev/null"

sys.stderr.write("Running hmmsearch...")
sys.stderr.flush()
Exec(CmdLine)
sys.stderr.write(" done.\n")

HMMs = set()
Labels = set()
LabelToHMM = {}
LabelToE = {}
sys.stderr.write("Parsing domtbl...")
sys.stderr.flush()
for Line in open(DomTbl):
	if Line.startswith('#'):
		continue
	Fields = Line[:-1].split()
	Label = Fields[0]
	HMM = Fields[3].split('.')[0]
	E = float(Fields[11])
	if E > Args.evalue:
		continue
	HMMs.add(HMM)
	if Label not in Labels:
		Labels.add(Label)
		LabelToHMM[Label] = HMM
		LabelToE[Label] = E
	elif E < LabelToE[Label]:
		LabelToHMM[Label] = HMM
		LabelToE[Label] = E
sys.stderr.write("%d HMM hits\n" % (len(Labels)))

Exec("rm -f " + DomTbl)

HMMToLabels = {}
for HMM in HMMs:
	HMMToLabels[HMM] = []
for Label in Labels:
	HMM = LabelToHMM[Label]
	HMMToLabels[HMM].append(Label)

def FindMotif(Seq, Motif):
	if Motif == "" or Motif == "-" or Motif == ".":
		return 0

	n = Seq.find(Motif)
	r = Seq.rfind(Motif)
	if n != r:
		return -1
	return n

def GetTrimSeq(Seq, PosA, PosB, PosC):
	if PosA > 0 and PosC > 0:
		Lo = PosA - 149
		Hi = PosC + 150
	elif PosC > 0:
		Lo = PosC = 200
		Hi = PosC + 150
	elif PosB > 0:
		Lo = PosB - 150
		Hi = PosC - 150
	elif PosA > 0:
		Lo = PosA - 149
		Hi = PosA + 200
	else:
		return None
	if Lo < 0:
		Lo = 0
	L = len(Seq)
	if Hi > L:
		Hi = L
	return Seq[Lo:Hi]

def OnChunkSeq(Label, Seq):
	SeqA, SeqB, SeqC = polhummer.GetMotifs(Seq)
	if SeqA == "" and SeqB == "" and SeqC == "":
		return
	FullSeq = ""
	for c in Seq:
		if c == "-" or c == ".":
			continue
		FullSeq += c.upper()
	E = LabelToE[Label]
	Label0 = Label.split()[0]
	FullPosA = FindMotif(FullSeq, SeqA)
	FullPosB = FindMotif(FullSeq, SeqB)
	FullPosC = FindMotif(FullSeq, SeqC)
	TrimSeq = GetTrimSeq(FullSeq, FullPosA, FullPosB, FullPosC)

	TrimPosA = FindMotif(TrimSeq, SeqA)
	TrimPosB = FindMotif(TrimSeq, SeqB)
	TrimPosC = FindMotif(TrimSeq, SeqC)

	ABCAnnot = "A:%d:%s" % (TrimPosA+1, SeqA)
	ABCAnnot += " B:%d:%s" % (TrimPosB+1, SeqB)
	ABCAnnot += " C:%d:%s" % (TrimPosC+1, SeqC)
	NewLabel = Label0 + " " + ABCAnnot + " " + "HMM:" + HMM
	if not fCore is None:
		fasta.WriteSeq(fCore, TrimSeq, NewLabel)

	FevRec = Label0
	FevRec += "\tmotif_hmm=%s/%.3g" % (HMM, E)
	FevRec += "\thmm_A=%s,%d" % (SeqA, FullPosA)
	FevRec += "\thmm_B=%s,%d" % (SeqB, FullPosB)
	FevRec += "\thmm_C=%s,%d" % (SeqC, FullPosC)
	if not fFev is None:
		fFev.write(FevRec + "\n")

for HMM in HMMs:
	sys.stderr.write(HMM + "\n")

	HMMFN = RepoDir + "motif_hmms/" + HMM + ".hmm"
	StoFN = RepoDir + "motif_stos/" + HMM + ".sto"
	if not os.path.isfile(HMMFN):
		sys.stderr.write("HMM not found %s\n" % HMMFN)
		assert False
	if not os.path.isfile(StoFN):
		sys.stderr.write("Stockholm file not found %s\n" % StoFN)
		assert False
	polhummer.ReadSto(StoFN)

	HitsFaFN = TmpPrefix + "hits.fa"
	HMMLabels = HMMToLabels[HMM]
	N = len(HMMLabels)
	M = 5000
	Lo = 0
	while True:
		if Lo >= N:
			break
		Hi = Lo + M
		if Hi > N:
			Hi = N
		f = open(HitsFaFN, "w")
		for i in range(Lo, Hi):
			Label = HMMLabels[i]
			Seq = SeqDict[Label]
			fasta.WriteSeq(f, Seq, Label)
		f.close()
		Lo = Hi

		ChunkA2mFN = TmpPrefix + "chunk.a2m"
		CmdLine = "hmmalign "
		CmdLine += " --outformat a2m"
		CmdLine += " " + HMMFN
		CmdLine += " " + HitsFaFN
		CmdLine += " > " + ChunkA2mFN
		Exec(CmdLine)

		fasta.ReadSeqsOnSeq(ChunkA2mFN, OnChunkSeq)

		Exec("rm -f " + ChunkA2mFN)
		Exec("rm -f " + HitsFaFN)

if not fCore is None:
	fCore.close()
if not fFev is None:
	fFev.close()
