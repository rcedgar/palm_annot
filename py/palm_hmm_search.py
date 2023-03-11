#!/usr/bin/env python3

import sys
import argparse
import os
import random
import fasta
import polhummer

Usage = \
(
"Classify sequences using two HMM libraries (1) RdRp and (2) non-RdRp homologs e.g. RT\n"
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

AP.add_argument("--tmpdir",
  required=False,
  default="/tmp",
  help="Directory for temporary files")

AP.add_argument("--evalue",
  required=False,
  default=1e-6,
  help="Max E-value for HMM search (default 1e-6)")

AP.add_argument("--dbsize",
  required=False,
  default=100000,
  help="Effective db size for E-value (-Z option of hmmsearch, default 100000)")

AP.add_argument("--threads",
  type=int,
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

pid = os.getpid()
r = random.randint(0, 999999)
TmpPrefix = TmpDir + "palm_annot.pmhmm.%d.%d." % (pid, r)
sys.stderr.write("TmpPrefix = %s\n" % TmpPrefix)

fFev = open(Args.fev , "w")

Labels = set()
LabelToHMMp = {}
LabelToHMMm = {}
LabelToLoHim = {}
LabelToLoHip = {}
LabelToEp = {}
LabelToEm = {}
DomTbl = TmpPrefix + "ph.domtbl"
for rdrp_xxx in [ "rdrp_plus", "rdrp_minus" ]:
	HMMDb = RepoDir + "hmmdbs/" + rdrp_xxx

	CmdLine = "hmmsearch"
	CmdLine += " --domtbl " + DomTbl
	if not Args.threads is None:
		CmdLine += " --cpu %d" % Args.threads
	CmdLine += " -E %.3g" % Args.evalue
	CmdLine += " -Z %d" % Args.dbsize
	CmdLine += "  " + HMMDb
	CmdLine += "  " + Args.input
	CmdLine += " > /dev/null"

	sys.stderr.write("Running hmmsearch %s..." % rdrp_xxx)
	sys.stderr.flush()
	Exec(CmdLine)
	sys.stderr.write(" done.\n")

	sys.stderr.write("Parsing domtbl...")
	sys.stderr.flush()
	UniqueLabels = set()
	for Line in open(DomTbl):
		if Line.startswith('#'):
			continue
		Fields = Line[:-1].split()
		Label = Fields[0].split()[0]
		HMM = Fields[3].split('.')[0]
		E = float(Fields[11])
		if E > Args.evalue:
			continue
		Lo = int(Fields[19])
		Hi = int(Fields[20])

		UniqueLabels.add(Label)
		if Label not in Labels:
			Labels.add(Label)
			LabelToHMMp[Label] = "."
			LabelToHMMm[Label] = "."
			LabelToEp[Label] = 999
			LabelToEm[Label] = 999
			LabelToLoHip[Label] = (0, 0)
			LabelToLoHim[Label] = (0, 0)

		if rdrp_xxx == "rdrp_plus":
			if E < LabelToEp[Label]:
				LabelToEp[Label] = E
				LabelToHMMp[Label] = HMM
				LabelToLoHip[Label] = (Lo, Hi)
		elif rdrp_xxx == "rdrp_minus":
			if E < LabelToEm[Label]:
				LabelToEm[Label] = E
				LabelToHMMm[Label] = HMM
				LabelToLoHim[Label] = (Lo, Hi)
		else:
			assert False
	sys.stderr.write("%d %s hits\n" % (len(UniqueLabels), rdrp_xxx))

	Exec("rm -f " + DomTbl)

for Label in Labels:
	HMMp = LabelToHMMp[Label]
	HMMm = LabelToHMMm[Label]
	Ep = LabelToEp[Label]
	Em = LabelToEm[Label]
	(Lop, Hip) = LabelToLoHip[Label]
	(Lom, Him) = LabelToLoHim[Label]

	if Em > Args.evalue and Ep > Args.evalue:
		continue

	FevRec = Label
	if Ep <= Args.evalue:
		FevRec += "\thmm_rdrp_plus=%s" % (HMMp)
		FevRec += "\thmm_rdrp_plus_evalue=%.3g" % (Ep)
		FevRec += "\thmm_rdrp_plus_lo=%d" % (Lop)
		FevRec += "\thmm_rdrp_plus_hi=%d" % (Hip)
	if Em <= Args.evalue:
		FevRec += "\thmm_rdrp_minus=%s" % (HMMm)
		FevRec += "\thmm_rdrp_minus_evalue=%.3g" % (Em)
		FevRec += "\thmm_rdrp_minus_lo=%d" % (Lom)
		FevRec += "\thmm_rdrp_minus_hi=%d" % (Him)
	fFev.write(FevRec + "\n")

fFev.close()
