#!/usr/bin/env python3

########################################################
########################################################
##	Cmd += " -exclude RT+RT2"
########################################################
########################################################

import sys
import argparse
import os
from time import gmtime, strftime

Usage = \
(
"TODO"
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
  help="Input RdRp/RT sequences (FASTA format)")

AP.add_argument("--maxhmms",
  required=False,
  type=int,
  default=32,
  help="Max HMMs to generate (default 32)")

AP.add_argument("--pssms",
  required=False,
  help="Initial PSSMs (default pssms/palm.ppm)")

AP.add_argument("--outdir",
  required=True,
  help="Output directory (must not exist)")

AP.add_argument("--evalue",
  required=False,
  default=0.001,
  help="Max E-value (default 0.001)")

AP.add_argument("--PH",
  required=False,
  default="PH",
  choices = [ "PH", "PN" ],
  help="Name prefix, PH for RdRp PN for non-RdRp (default PH)")

AP.add_argument("--minseqs",
  required=False,
  type=int,
  default=3,
  help="Min centroids at 75%% id to build HMM (default 3)")

AP.add_argument("--minhits",
  required=False,
  type=int,
  default=10,
  help="Min hits to continue iterating (default 3)")

AP.add_argument("--minhitpct",
  required=False,
  type=float,
  default=5,
  help="Min hit percent to continue iterating (default 3)")

AP.add_argument("--maxseed",
  required=False,
  type=int,
  default=500,
  help="Max seqs per seed alignment")

AP.add_argument("--maxflank",
  required=False,
  type=int,
  default=150,
  help="Max letters to include in pp flanks (default 150)")

AP.add_argument("--minflank",
  required=False,
  type=int,
  default=75,
  help="Min letters to include in pp flanks (default 75)")

AP.add_argument("--mcs",
  required=False,
  type=float,
  default=2.0,
  help="Minimum cluster score (default 2.0)")

AP.add_argument("--mps",
  required=False,
  type=float,
  default=20.0,
  help="Minimum palm score (default 20.0)")

AP.add_argument("--threads",
  type=int,
  required=False,
  help="Number of threads for hmmsearch --cpu option (default not set)")

AP.add_argument("--include",
  required=False,
  help="Include only these PSSMs, exclude others (e.g. RT+RT2, default use all")

AP.add_argument("--exclude",
  required=False,
  help="Exclude these PSSMs (e.g. RT+RT2, default use all")

Args = AP.parse_args()

if os.path.isdir(Args.outdir):
	sys.stderr.write("\n\nERROR: output directory already exists %s\n\n" % Args.outdir)
	sys.exit(1)

if (not Args.include is None) and (not Args.exclude is None):
	sys.stderr.write("\n\nERROR: --includes and --excludes\n\n")
	sys.exit(1)

# Need -p in case creating multiple levels
os.system("mkdir -p " + Args.outdir)

PSSMFN = RepoDir + "pssms/palm.ppm"
if not Args.pssms is None:
	PSSMFN = Args.pssms

PH = Args.PH

IterIndex = 0
fLog = None
LogFN = Args.outdir + "/build_hmms.log"
StderrFN = Args.outdir + "/stderr.txt"
def Log(s):
	global fLog
	ts = strftime("%Y-%m-%d %H:%M:%S", gmtime())
	if fLog is None:
		fLog = open(LogFN, "w")
		s = ts
		for Arg in sys.argv:
			s += " " + Arg
		fLog.write(s + "\n")
	ml = s.replace(' -', '\n        -')
	sys.stdout.write("[" + str(IterIndex) + "] " + ml + "\n")
	fLog.write(ts + " [" + str(IterIndex) + "] " + s + "\n")

def Exec(Cmd):
	assert Cmd[0] != '/'
	Log(Cmd)
	Redir = " 2>> %s" % StderrFN
	Code = os.system(Cmd + Redir)
	if Code != 0:
		Log("Exec failed rc=%d\n" % Code)
		sys.stderr.write("\n")
		sys.stderr.write(Cmd + "\n")
		sys.stderr.write("\n")
		sys.stderr.write("Error code %d\n" % Code)
		assert False

def GetFileSize(FN):
	try:
		Bytes = os.path.getsize(FN)
	except:
		Bytes = -1
	return Bytes

IterDir = Args.outdir + "/iters/"
Iter0Dir = IterDir + "/iter0/"
StoDir = Args.outdir + "/sto/"
HmmDir = Args.outdir + "/hmm/"

Exec("mkdir " + StoDir)
Exec("mkdir " + HmmDir)
Exec("mkdir " + IterDir)
Exec("mkdir " + Iter0Dir)

Exec("cp %s %s/miss.fa" % (Args.input, Iter0Dir))

def GetSeqCount(FN):
	try:
		f = open(FN)
	except:
		return 0
	n = 0
	for Line in f:
		if Line.startswith(">"):
			n += 1
	f.close()
	return n

def MissingOrEmptyFile(FN):
	Bytes = GetFileSize(FN)
	if Bytes <= 0:
		return True
	return False

def Iter(IterIndex):
	IterN_1Dir = IterDir + "iter" + str(IterIndex-1) + "/"
	IterNDir = IterDir + "iter" + str(IterIndex) + "/"
	Exec("mkdir " + IterNDir)

	IterInputFN = "%s/miss.fa" % IterN_1Dir
	InputSeqCount = GetSeqCount(IterInputFN)
	SeedSize = InputSeqCount
	if SeedSize > Args.maxseed:
		SeedSize = Args.maxseed

	Log("%d input seqs from %s" % (InputSeqCount, IterInputFN))
	Log("Seed size %d" % SeedSize)

	StoFN = StoDir + PH + str(IterIndex) + ".sto"
	HmmFN = HmmDir + PH + str(IterIndex) + ".hmm"
	DomFN = IterNDir + "hmmsearch.domtbl"
	LabelsFN = IterNDir + "labels.txt"

###############################################################
# Cluster motifs, report palmprint segments for largest cluster
###############################################################
	Cmd= "palmscan2"
	Cmd += " -cluster_motifs_greedy " + IterInputFN
	Cmd += " -model " + PSSMFN
	Cmd += " -motif_cluster_minscore %.1f" % Args.mcs
	Cmd += " -cluster_tsv " + IterNDir + "cluster.tsv"
	Cmd += " -min_palm_score %.1f" % Args.mps
	Cmd +="  -topn 1"
	Cmd += " -seg_fasta_prefix " + IterNDir + "seg_"
	if not Args.include is None:
		Cmd += " -include %s" % Args.include
	if not Args.exclude is None:
		Cmd += " -exclude %s" % Args.exclude
	Cmd += " -quiet"
	Exec(Cmd)

##############################################
# Cluster at 75% id to choose ABC or CAB order
##############################################
	seg_ABC_PP = IterNDir + "seg_ABC_PP"
	seg_CAB_PP = IterNDir + "seg_CAB_PP"
	NABC_PP = 0
	NCAB_PP = 0
	if GetFileSize(seg_ABC_PP) > 0:
		Cmd = "usearch"
		Cmd += " -cluster_fast " + seg_ABC_PP
		Cmd += " -id 0.75"
		Cmd += " -centroids " + seg_ABC_PP + ".id75"
		Cmd += " -quiet"
		Exec(Cmd)
		
		NABC_PP = GetSeqCount(seg_ABC_PP + ".id75")

	if GetFileSize(seg_CAB_PP) > 0:
		Cmd = "usearch"
		Cmd += " -cluster_fast " + seg_CAB_PP
		Cmd += " -id 0.75"
		Cmd += " -centroids " + seg_CAB_PP + ".id75"
		Cmd += " -quiet"
		Exec(Cmd)
		
		NCAB_PP = GetSeqCount(seg_CAB_PP + ".id75")

	if NABC_PP < Args.minseqs and NCAB_PP < Args.minseqs:
		Log("Iter %d terminating, NABC_PP=%d, NABC_PP=%d\n" \
		  % (IterIndex, NABC_PP, NCAB_PP))
		return False

	if NABC_PP >= NCAB_PP:
		ABC = "ABC"
	else:
		ABC = "CAB"

####################################################
# Motifs already aligned
#	(*) Uniques => X.uniques
#	(*) Subsample => X.sub
####################################################
	for X in [ "A", "B", "C" ]:
		SegFN = IterNDir + "seg_" + ABC + "_" + X
		UniquesFN = IterNDir + X + ".uniques"
		SubFN = IterNDir + X + ".sub"

		Cmd = "usearch"
		Cmd += " -fastx_uniques " + SegFN
		Cmd += " -fastaout " + UniquesFN
		Cmd += " -quiet"
		Exec(Cmd)

		Cmd = "usearch"
		Cmd += " -fastx_subsample " + UniquesFN
		Cmd += " -sample_size %d" % SeedSize
		Cmd += " -fastaout " + SubFN
		Cmd += " -quiet"
		Exec(Cmd)

####################################################
# Palmprint segments V1 and V2
# No padding required, no empty files/sequences
#	(*) Uniques => X.uniques
#	(*) Subsample => X.sub
#	(*) Align => X.afa
####################################################
	for X in [ "V1", "V2" ]:
		SegFN = IterNDir + "seg_" + ABC + "_" + X
		UniquesFN = IterNDir + X + ".uniques"
		PadxFN = IterNDir + X + ".padx"
		SubFN = IterNDir + X + ".sub"
		AfaFN = IterNDir + X + ".afa"

		Cmd = "usearch"
		Cmd += " -fastx_uniques " + SegFN
		Cmd += " -fastaout " + UniquesFN
		Cmd += " -quiet"
		Exec(Cmd)

		Cmd = "usearch"
		Cmd += " -fastx_subsample " + UniquesFN
		Cmd += " -sample_size %d" % SeedSize
		Cmd += " -fastaout " + SubFN
		Cmd += " -quiet"
		Exec(Cmd)

		Cmd = "muscle"
		Cmd += " -super5 " + SubFN
		Cmd += " -output " + AfaFN
		Cmd += " -quiet"
		Exec(Cmd)

		SegFN = None
		UniquesFN = None
		PadxFN = None
		SubFN = None
		AfaFN = None

###################################################################
# Left / right palmprint flanking sequences
# Pad with Xs to encourage HMM match states rather than inserts
#	(*) Uniques => X.uniques
#	(*) Subsample => X.sub
#	(*) Pad with Xs => X.padx to desired flank length
#	(*) Trim 5' and 3' until >=3 sequences with non-X
#	(*) Align => X.afa (may be empty)
###################################################################
	for X in [ "Left", "Right" ]:
		SegFN = IterNDir + "seg_" + ABC + "_" + X
		UniquesFN = IterNDir + X + ".uniques"
		SubFN = IterNDir + X + ".sub"
		PadxFN = IterNDir + X + ".padx"
		AfaFN = IterNDir + X + ".afa"
		if MissingOrEmptyFile(SegFN):
			Log("Missing or empty " + SegFN)
			Exec("touch " + AfaFN)
			continue

		Cmd = "usearch"
		Cmd += " -fastx_uniques " + SegFN
		Cmd += " -fastaout " + UniquesFN
		Cmd += " -quiet"
		Exec(Cmd)

		n = GetSeqCount(UniquesFN)
		if n > SeedSize:
			n = SeedSize

	# Extend to desired flank length (150)
	# Discard if flank to short (<75)
		Cmd = "polhummer_fasta_subsample_and_padx.py"
		Cmd += " --input " + UniquesFN
		Cmd += " --minlength %d" % Args.minflank
		Cmd += " --trimlength %d" % Args.maxflank
		Cmd += " --maxn " + str(n)
		Cmd += " --end " + X
		Cmd += " --output " + PadxFN
		Exec(Cmd)

	# Trim 5' and 3' until >=3 sequences with non-X
		Cmd = "polhummer_fasta_trimx.py"
		Cmd += " --input " + PadxFN
		Cmd += " --end " + X
		Cmd += " --minseqs 3"
		Cmd += " --output " + SubFN
		Exec(Cmd)

		if GetSeqCount(SubFN) == 0:
			Exec("touch " + AfaFN)
		else:
			Cmd = "muscle"
			Cmd += " -super5 " + SubFN
			Cmd += " -output " + AfaFN
			Cmd += " -quiet"
			Exec(Cmd)

		SegFN = None
		UniquesFN = None
		SubFN = None
		PadxFN = None

###################################################################
# Assemble segment alignments into seed alignment
#	Motifs are are aligned to PSSMs, fixed length.
#	Other segments aligned by Super5.
#	Stockholm format annotates motif columns.
###################################################################
	AsubFN = IterNDir + "A.sub"
	BsubFN = IterNDir + "B.sub"
	CsubFN = IterNDir + "C.sub"
	V1afaFN = IterNDir + "V1.afa"
	V2afaFN = IterNDir + "V2.afa"
	LeftafaFN = IterNDir + "Left.afa"
	RightafaFN = IterNDir + "Right.afa"
	Cmd = "polhummer_assemble.py"
	Cmd += " --name %s%d" % (PH, IterIndex)
	Cmd += " --afa_A " + AsubFN
	Cmd += " --afa_B " + BsubFN
	Cmd += " --afa_C " + CsubFN
	Cmd += " --afa_V1 " + V1afaFN
	Cmd += " --afa_V2 " + V2afaFN
	Cmd += " --afa_Left " + LeftafaFN
	Cmd += " --afa_Right " + RightafaFN
	Cmd += " --order " + ABC
	Cmd += " --motif_msa_prefix " + IterNDir + "motif"
	Cmd += " --output " + StoFN
	Exec(Cmd)

###################################################################
# Build HMM from seed
###################################################################
	Cmd = "hmmbuild"
	Cmd += " --hand"
	Cmd += " -n %s%d" % (PH, IterIndex)
	Cmd += " " + HmmFN
	Cmd += " " + StoFN
	Cmd += " > /dev/null"
	Exec(Cmd)

###################################################################
# Search input sequences with new HMM
###################################################################
	Cmd = "hmmsearch"
	Cmd += " --domtblout " + DomFN
	Cmd += " -E %.4g" % Args.evalue
	Cmd += " -Z 1000000"
	if not Args.threads is None:
		Cmd += " --cpu %d" % Args.threads
	Cmd += " " + HmmFN
	Cmd += " " + IterInputFN
	Cmd += " > /dev/null"
	Exec(Cmd)

###################################################################
# Extract hit labels from HMM output
###################################################################
	fLabels = open(LabelsFN, "w")
	HitCount = 0
	for Line in open(DomFN):
		if Line.startswith("#"):
			continue
		Fields = Line[:-1].split()
		E = float(Fields[12])
		if E > Args.evalue:
			continue
		Label = Fields[0]
		fLabels.write(Label + "\n")
		HitCount += 1
	fLabels.close()
	Log("%d hits" % HitCount)

####################################################################
# Extract sequences which don't match HMM = input for next iteration
####################################################################
	HitsFN = IterNDir + "hits.fa"
	MissFN = IterNDir + "miss.fa"
	Cmd = "usearch"
	Cmd += " -fastx_getseqs " + IterInputFN
	Cmd += " -labels " + LabelsFN
	Cmd += " -fastaout " + HitsFN
	Cmd += " -notmatched " + MissFN
	Cmd += " -quiet"
	Exec(Cmd)

	Log("%d input, %d hits\n" % (InputSeqCount, HitCount))
	HitSeqCount = GetSeqCount(HitsFN)
	if HitSeqCount < Args.minhits:
		Log("Terminating, hits %d < min hits %d" % (HitSeqCount, Args.minhits))
		return False

	HitPct = HitSeqCount*100.0/InputSeqCount
	if HitPct < Args.minhitpct:
		Log("Terminating, hits %.3g%% < min pct %.3g%%" % (HitPct, Args.minhitpct))
		return False

	return True

for IterIndex in range(1, Args.maxhmms):
	Ok = Iter(IterIndex)
	if not Ok:
		break

Log("Done, %d HMMs created" % (IterIndex-1))
