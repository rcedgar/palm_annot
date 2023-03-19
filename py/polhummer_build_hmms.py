#!/usr/bin/env python3

########################################################
########################################################
##	Cmd += " -exclude RT,RT2"
########################################################
########################################################

import sys
import argparse
import os
import random
import fasta
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
  default=1e-3,
  help="Max E-value (default 0.001)")

AP.add_argument("--minseqs",
  required=False,
  type=int,
  default=3,
  help="Min centroids at 75%% id to build HMM (default 3)")

AP.add_argument("--maxseed",
  required=False,
  type=int,
  default=500,
  help="Max seqs per seed alignment")

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

Args = AP.parse_args()

if os.path.isdir(Args.outdir):
	sys.stderr.write("\n\nERROR: output directory already exists %s\n\n" % Args.outdir)
	sys.exit(1)

os.system("mkdir -p " + Args.outdir)

PSSMFN = RepoDir + "pssms/palm.ppm"
if not Args.pssms is None:
	PSSMFN = Args.pssms

fLog = None
LogFN = Args.outdir + "/build_hmms.log"
def Log(s):
	global fLog
	ts = strftime("%Y-%m-%d %H:%M:%S", gmtime())
	if fLog is None:
		fLog = open(LogFN, "w")
		s = ts
		for Arg in sys.argv:
			s += " " + Arg
		fLog.write(s + "\n")
	sys.stdout.write(s + "\n")
	fLog.write(ts + " " + s + "\n")

def Exec(Cmd):
	assert Cmd[0] != '/'
	Log(Cmd)
	Code = os.system(Cmd)
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
A2mDir = Args.outdir + "/a2m/"

Exec("mkdir " + StoDir)
Exec("mkdir " + A2mDir)
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
	Log("Iter %d, %d seqs" % (IterIndex, InputSeqCount))

	SeedSize = InputSeqCount
	if SeedSize > Args.maxseed:
		SeedSize = Args.maxseed

	StoFN = IterNDir + "PH" + str(IterIndex) + ".sto"
	HmmFN = IterNDir + "PH" + str(IterIndex) + ".hmm"
	DomFN = IterNDir + "hmmsearch.domtbl"
	LabelsFN = IterNDir + "labels.txt"

	Cmd= "palmscan2"
	Cmd += " -cluster_motifs_greedy " + IterInputFN
	Cmd += " -model " + PSSMFN
	Cmd += " -motif_cluster_minscore %.1f" % Args.mcs
	Cmd += " -cluster_tsv " + IterNDir + "cluster.tsv"
	Cmd += " -min_palm_score %.1f" % Args.mps
	Cmd +="  -topn 1"
	Cmd += " -seg_fasta_prefix " + IterNDir + "seg_"
	Cmd += " -exclude RT,RT2" ########################################################
	Exec(Cmd)

	seg_ABC_PP = IterNDir + "seg_ABC_PP"
	seg_CAB_PP = IterNDir + "seg_CAB_PP"
	NABC_PP = 0
	NCAB_PP = 0

	if GetFileSize(seg_ABC_PP) > 0:
		Cmd = "usearch"
		Cmd += " -cluster_fast " + seg_ABC_PP
		Cmd += " -id 0.75"
		Cmd += " -centroids " + seg_ABC_PP + ".id75"
		Exec(Cmd)
		
		NABC_PP = GetSeqCount(seg_ABC_PP + ".id75")

	if GetFileSize(seg_CAB_PP) > 0:
		Cmd = "usearch"
		Cmd += " -cluster_fast " + seg_CAB_PP
		Cmd += " -id 0.75"
		Cmd += " -centroids " + seg_CAB_PP + ".id75"
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

	for X in [ "A", "B", "C", "V1", "V2" ]:
		FN = IterNDir + "seg_" + ABC + "_" + X
		Cmd = "usearch"
		Cmd += " -fastx_uniques " + FN
		Cmd += " -fastaout " + FN + ".uniques"
		Exec(Cmd)

		Cmd = "usearch"
		Cmd += " -fastx_subsample " + FN + ".uniques"
		Cmd += " -sample_size %d" % SeedSize
		Cmd += " -fastaout " + FN + ".subx"
		Exec(Cmd)

	for Left in [ "Left", "Right" ]:
		FN = IterNDir + "seg_" + ABC + "_" + Left
		if MissingOrEmptyFile(FN):
			Exec("touch " + FN)
			continue

		Cmd = "usearch"
		Cmd += " -fastx_uniques " + FN
		Cmd += " -fastaout " + FN + ".uniques"
		Exec(Cmd)

		n = GetSeqCount(FN)
		if n > SeedSize:
			n = SeedSize

		SubFN = FN + ".sub"
		SubxFN = FN + ".subx"

		Cmd = "fasta_subsample_and_padx.py"
		Cmd += " --input " + FN + ".uniques"
		Cmd += "  --minlength 75"
		Cmd += " --trimlength 150"
		Cmd += " --maxn " + str(n)
		Cmd += " --end " + Left
		Cmd += " > " + SubFN
		Exec(Cmd)

		Cmd = "fasta_trimx.py"
		Cmd += " --input " + SubxFN
		Cmd += " --end " + Left
		Cmd += " --minseqs 3"
		Cmd += " > " + SubFN
		Exec(Cmd)

	for X in [ "Left", "Right", "V1", "V2" ]:
		SubxFN = "seg_" + ABC + "_" + X + ".subx"
		SubFN = "seg_" + ABC + "_" + X + ".sub"
		AfaFN = "seg_" + ABC + "_" + X + ".afa"
		if GetSeqCount(SubxFN) == 0:
			Exec("touch " + SubFN)
			continue

		Cmd = "muscle"
		Cmd += " -super5 " + SubFN
		Cmd += " -output " + AfaFN
		Exec(Cmd)

	Cmd = "polhummer_assemble_segafas.py"
	Cmd += " --name PH%d" % IterIndex
	Cmd += " --afa_A seg_" + ABC + "_A.sub"
	Cmd += " --afa_B seg_" + ABC + "_B.sub"
	Cmd += " --afa_C seg_" + ABC + "_C.sub"
	Cmd += " --afa_V1 seg_" + ABC + "_V1.afa"
	Cmd += " --afa_V2 seg_" + ABC + "_V2.afa"
	Cmd += " --afa_Left seg_" + ABC + "_Left.afa"
	Cmd += " --afa_Right seg_" + ABC + "_Right.afa"
	Cmd += " --order " + ABC
	Cmd += " > " + StoFN
	Exec(Cmd)

	Cmd = "hmmbuild"
	Cmd += " --hand"
	Cmd += " -n PH" + str(IterIndex)
	Cmd += " " + HmmFN
	Cmd += " " + StoFN
	Exec(Cmd)

	Cmd = "hmmsearch"
	Cmd += " --domtblout " + DomFN
	Cmd += " --evalue %.4g" % Args.evalue
	Cmd += " --Z 1000000"
	if not Args.threads is None:
		Cmd += " --cpu %d" % Args.threads
	Cmd += " " + IterInputFN
	Cmd += " " + HmmFN
	Cmd += " > /dev/null"
	Exec(Cmd)

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

	Cmd = "usearch"
	Cmd += " -fastx_getseqs " + IterInputFN
	Cmd += " -labels " + LabelsFN
	Cmd += " -notmatched miss.fa"
	Exec(Cmd)

for IterIndex in range(1, Args.maxhmms):
	Ok = Iter(IterIndex)
	if not Ok:
		break
