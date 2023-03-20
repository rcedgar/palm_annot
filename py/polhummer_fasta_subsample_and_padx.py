#!/usr/bin/python3

import sys
import argparse
import fasta

AP = argparse.ArgumentParser()
AP.add_argument("--input", required=True, help="Input FASTA")
AP.add_argument("--output", required=True, help="Output FASTA")
AP.add_argument("--minlength", required=True, type=int, help="Minimum seq length, shorter discarded")
AP.add_argument("--trimlength", required=True, type=int, help="Longer truncated, shorter padded with Xs")
AP.add_argument("--end", required=True, choices=[ "Left", "Right"], help="Left or right")
AP.add_argument("--maxn", required=True, type=int, help="Maximum number to output, longer first")
Args = AP.parse_args()

Triples = []

fOut = open(Args.output, "w")

def OnSeq(Label, Seq):
	L = len(Seq)
	if L < Args.minlength:
		return
	if L > Args.trimlength:
		L = Args.trimlength
		Seq = Seq[:L]
		assert len(Seq) == L
	Pad = Args.trimlength - L
	if Pad > 0:
		Xs = "X"*Pad
		if Args.end == "Left":
			Seq = Xs + Seq
		elif Args.end == "Right":
			Seq = Seq + Xs
		else:
			assert False

	Triple = (Label, Seq, L)
	Triples.append(Triple)

fasta.ReadSeqsOnSeq(Args.input, OnSeq)

SortedTriples = sorted(Triples, key=lambda Triple: Triple[2])
M = len(SortedTriples)
N = min(Args.maxn, M)
for i in range(N):
	Label, Seq, L = SortedTriples[M-i-1]
	fasta.WriteSeq(fOut, Seq, Label)
