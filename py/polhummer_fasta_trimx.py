#!/usr/bin/python3

import sys
import argparse
import fasta

AP = argparse.ArgumentParser()
AP.add_argument("--input", required=True, help="Input FASTA")
AP.add_argument("--output", required=True, help="Output FASTA")
AP.add_argument("--minseqs", required=False, type=int, default=1, help="Minimum number of sequences with non-X")
AP.add_argument("--end", required=True, choices=[ "Left", "Right"], help="Left or right")
Args = AP.parse_args()

fOut = open(Args.output, "w")

MinNX = None
MaxNX = None
def GetNX(Seq):
	global MinNX, MaxNX
	NX = 0
	if Args.end == "Left":
		for c in Seq:
			if c != 'X':
				break
			NX += 1
	elif Args.end == "Right":
		for c in Seq[::-1]:
			if c != 'X':
				break
			NX += 1
	else:
		assert False
	if MinNX == None:
		MinNX = NX
		MaxNX = NX
	else:
		MinNX = min(NX, MinNX)
		MaxNX = max(NX, MaxNX)
	return NX

Labels = []
Seqs = []
NXs = []
def OnSeq(Label, Seq):
	NX = GetNX(Seq)
	Labels.append(Label)
	Seqs.append(Seq)
	NXs.append(NX)
fasta.ReadSeqsOnSeq(Args.input, OnSeq)
N = len(Labels)

# Number of Xs to trim
def GetM():
	if N == 0:
		return 0
	SortedNXs = NXs[:]
	SortedNXs.sort()
	n = Args.minseqs
	n = min(n, len(SortedNXs))
	M = SortedNXs[n-1]
	return M

M = GetM()

ND = 0
NOut = 0
for i in range(N):
	Label = Labels[i]
	Seq = Seqs[i]

	for i in range(M):
		if Seq[-1] == 'X':
			Seq = Seq[:-1]
	if len(Seq) == 0:
		ND += 1
		continue
	NOut += 1
	fasta.WriteSeq(fOut, Seq, Label)

if N == 0:
	sys.stderr.write("polhummer_fasta_trimx.py, empty input %s\n" % Args.input)
else:
	sys.stderr.write("NX %d .. %d, M %d, input %d, deleted %d, output %d\n" % (MinNX, MaxNX, M, N, ND, NOut))
