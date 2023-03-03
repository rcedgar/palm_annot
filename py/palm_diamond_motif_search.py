#!/usr/bin/env python3

import sys
import argparse
import os
import random
import cigar

Usage = \
(
"Classify sequences by diamond search against RdRp search database with motif annotations.\n"
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

AP.add_argument("--fev",
  required=True,
  help="Results in field-equals-value (fev) format")

AP.add_argument("--evalue",
  required=False,
  default=1e-9,
  help="Max E-value for HMM search (default 1e-9)")

AP.add_argument("--threads",
  required=False,
  help="Number of hmmsearch threads (default --threads option of diamond not set)")

AP.add_argument("--tmpdir",
  required=False,
  default="/tmp",
  help="Directory for temporary files")

AP.add_argument("--sensitive",
  required=False,
  default="very-sensitive",
  choices=[ "fast", "midsensitive", "more-sensitive", "very-sensitive"],
  help="diamond sensitivity option")

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
TmpPrefix = TmpDir + "pdms%d.%d." % (pid, r)
sys.stderr.write("TmpPrefix = %s\n" % TmpPrefix)

TsvFN = TmpPrefix + ".tsv"

RefDb = RepoDir + "diamond_refdbs/rdrp_plus_abc"

fFev = open(Args.fev , "w")

CmdLine = "diamond blastp"
CmdLine += " --query " + Args.input
CmdLine += " --evalue %.3g" % Args.evalue
CmdLine += "  --max-target-seqs 1"
CmdLine += " --%s" % Args.sensitive
CmdLine += " --db " + RefDb
if not Args.threads is None:
	CmdLine += " --cpu %d" % Args.threads
CmdLine += " --outfmt 6 qseqid qstart qend sseqid evalue qseq sseq cigar"
CmdLine += " > " + TsvFN

sys.stderr.write("Running diamond ...")
sys.stderr.flush()
Exec(CmdLine)
sys.stderr.write(" done.\n")

sys.stderr.write("Parsing hits...")
sys.stderr.flush()

# rdva.ERR1301806_NODE_722_length_1097_cov_42.330595;A=IALDLNKSDTSL;B=SGSMFTNIITTIAT;C=TYGDNMYV
def GetABC(Label):
	A = None
	B = None
	C = None
	Fields = Label.split(";")
	for Field in Fields:
		if Field.startswith("A="):
			A = Field[2:]
			assert len(A) == 12
		elif Field.startswith("B="):
			B = Field[2:]
			assert len(B) == 14
		elif Field.startswith("C="):
			C = Field[2:]
			assert len(C) == 8
	if A is None or B is None or C is None:
		sys.stderr.write("Missing motif(s) in db label >%s\n" % Label)
		assert False
	return A, B, C

def GetMotif(QRow, TRow, TPos, n):
	Motif = ""
	j = 0
	ColCount = len(QRow)
	assert len(TRow) == ColCount
	for Col in range(0, ColCount):
		q = QRow[Col]
		t = TRow[Col]
		if j >= TPos and q != '-':
			Motif += q
			if len(Motif) == n:
				return Motif
		if t != '-':
			j += 1
	return None

for Line in open(TsvFN):
#      0      1    2      3      4    5    6     7
# qseqid qstart qend sseqid evalue qseq sseq cigar
	Fields = Line[:-1].split()
	assert len(Fields) == 8
	E = float(Fields[4])
	if E > Args.evalue:
		continue
	QueryLabel = Fields[0]
	Lo = int(Fields[1])
	Hi = int(Fields[2])
	TargetLabel = Fields[3]
	QSeg = Fields[5]
	TSeg = Fields[6]
	CIGAR = Fields[7]

	cigar.ParseCigar(CIGAR)
	QRow = cigar.GetReadRow(QSeg)
	TRow = cigar.GetTargetRow(TSeg)

	TA, TB, TC = GetABC(TargetLabel)

	SeqA = None
	SeqB = None
	SeqC = None

	TPosA = TSeg.find(TA)
	TPosB = TSeg.find(TB)
	TPosC = TSeg.find(TC)
	if TPosA >= 0:
		SeqA = GetMotif(QRow, TRow, TPosA, 10)
	if TPosB >= 0:
		SeqB = GetMotif(QRow, TRow, TPosB, 12)
	if TPosC >= 0:
		SeqC = GetMotif(QRow, TRow, TPosC, 8)

	FevRec = QueryLabel
	FevRec += "\tdmnd_evalue=%.3g" % E
	FevRec += "\tdmnd_lo=%d" % Lo
	FevRec += "\tdmnd_hi=%d" % Hi
	if not SeqA is None:
		FevRec += "\tdmnd_A=%s" % SeqA
	if not SeqB is None:
		FevRec += "\tdmnd_B=%s" % SeqB
	if not SeqC is None:
		FevRec += "\tdmnd_C=%s" % SeqC

	fFev.write(FevRec + "\n")

fFev.close()
