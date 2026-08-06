#!/usr/bin/python3

import sys
import argparse
import fasta
import sortdict3

AL = 12
BL = 14
CL = 8

AP = argparse.ArgumentParser()
AP.add_argument("--afa_A", required=True)
AP.add_argument("--afa_B", required=True)
AP.add_argument("--afa_C", required=True)
AP.add_argument("--afa_V1", required=True)
AP.add_argument("--afa_V2", required=True)
AP.add_argument("--afa_Left", required=True)
AP.add_argument("--afa_Right", required=True)
AP.add_argument("--order", choices=["ABC", "CAB"], required=True)
AP.add_argument("--name", default="Seq")
AP.add_argument("--motif_msa_prefix", default="motif")
AP.add_argument("--output", required=True)
Args = AP.parse_args()

fOut = open(Args.output, "w")

W = 16
BLOCK = 80

N = 0
def RSD(FN):
	global N
	Seqs = []
	D = fasta.ReadSeqsDict(FN)
	Labels = list(D.keys())
	for Label in Labels:
		Seq = D[Label]
		Seqs.append(Seq)
	n = len(Seqs)
	N = max(n, N)
	if n == 0:
		n = 1
		Seqs.append("")
	return n, Seqs

N_A, Seqs_A = RSD(Args.afa_A)
N_B, Seqs_B = RSD(Args.afa_B)
N_C, Seqs_C = RSD(Args.afa_C)
N_V1, Seqs_V1 = RSD(Args.afa_V1)
N_V2, Seqs_V2 = RSD(Args.afa_V2)
N_Left, Seqs_Left = RSD(Args.afa_Left)
N_Right, Seqs_Right = RSD(Args.afa_Right)

sys.stderr.write("N = %d\n" % N)

def UngapL(Seq):
	n = 0
	for c in Seq:
		if c != '-' and c != '.':
			n += 1
	return n

Labels = []
As = []
Bs = []
Cs = []

MatchCount = None
ColCount = None
ColA = None
ColB = None
ColC = None

NLeft = None
NRight = None
NV1 = None
NV2 = None

Rows = []
for i in range(N):
	A = Seqs_A[i%N_A]
	B = Seqs_B[i%N_B]
	C = Seqs_C[i%N_C]
	V1 = Seqs_V1[i%N_V1]
	V2 = Seqs_V2[i%N_V2]
	Left = Seqs_Left[i%N_Left]
	Right = Seqs_Right[i%N_Right]

	assert len(A) == 12
	assert len(B) == 14
	assert len(C) == 8

	LeftUL = UngapL(Left)
	V1UL = UngapL(V1)
	V2UL = UngapL(V2)

	if Args.order == "ABC":
		PosA = LeftUL + 1
		PosB = PosA + AL + V1UL
		PosC = PosB + BL + V2UL
		Row = Left + A + V1 + B + V2 + C + Right
	elif Args.order == "CAB":
		PosC = LeftUL + 1
		PosA = PosC + CL + V1UL
		PosB = PosA + AL + V2UL
		Row = Left + C + V1 + A + V2 + B + Right
	else:
		assert False
	Rows.append(Row)

	if ColCount == None:
		ColCount = len(Row)
		NLeft = len(Left)
		NRight = len(Right)
		NV1 = len(V1)
		NV2 = len(V2)
		if Args.order == "ABC":
			ColA = NLeft
			ColB = ColA + 12 + NV1
			ColC = ColB + 14 + NV2
		elif Args.order == "CAB":
			ColC = NLeft
			ColA = ColC + 8 + NV1
			ColB = ColA + 12 + NV2
		else:
			assert False
	else:
		assert NLeft == len(Left)
		assert NRight == len(Right)
		assert NV1 == len(V1)
		assert NV2 == len(V2)
	assert ColCount == len(Row)
	assert ColCount == (NLeft + 12 + NV1 + 14 + NV2 + 8 + NRight)

	Label = "%s.%d" % (Args.name, i+1)
	Labels.append(Label)
	As.append(A)
	Bs.append(B)
	Cs.append(C)

N = len(Labels)
assert len(As) == N
assert len(Bs) == N
assert len(Cs) == N

f = open(Args.motif_msa_prefix + ".A", "w")
for i in range(N):
	f.write(">" + Labels[i] + "\n")
	f.write(As[i] + "\n")
f.close()

f = open(Args.motif_msa_prefix + ".B", "w")
for i in range(N):
	f.write(">" + Labels[i] + "\n")
	f.write(Bs[i] + "\n")
f.close()

f = open(Args.motif_msa_prefix + ".C", "w")
for i in range(N):
	f.write(">" + Labels[i] + "\n")
	f.write(Cs[i] + "\n")
f.close()

SeqCount = len(Rows)

def GetConsChar(Seqs, Col):
	N = len(Seqs)
	CharToCount = {}
	for i in range(N):
		c = Seqs[i][Col]
		sortdict3.IncCount(CharToCount, c)
	Order = sortdict3.GetOrder(CharToCount)
	Chars = list(CharToCount.keys())
	ConsChar = Chars[Order[0]]
	n = CharToCount[ConsChar]
	if n/N >= 0.5:
		return ConsChar
	elif n/N >= 0.25:
		return ConsChar.lower()
	return "x"

def GetConsSeq(Seqs):
	N = len(Seqs)
	if N == 0:
		return ""
	ConsSeq = ""
	L = len(Seqs[0])
	for i in range(L):
		c = GetConsChar(Seqs, i)
		ConsSeq += c
	return ConsSeq

ConsA = GetConsSeq(As)
ConsB = GetConsSeq(Bs)
ConsC = GetConsSeq(Cs)

ConsA = ConsA.replace("x", ".")
ConsB = ConsB.replace("x", ".")
ConsC = ConsC.replace("x", ".")

MatchCount = 0
IsMatch = []
ColToMatch = []
for Col in range(ColCount):
	GapCount = 0
	for i in range(N):
		c = Rows[i][Col]
		if c == '-':
			GapCount += 1
	Match = (float(GapCount)/N < 0.5)
	if Match:
		ColToMatch.append(MatchCount)
		MatchCount += 1
	else:
		ColToMatch.append(-1)
	IsMatch.append(Match)
assert len(ColToMatch) == ColCount

ConsRow = GetConsSeq(Rows)

'''
Feature   Description        Description
-------   -----------        --------------
RF        ReFerence annot.   Often the consensus RNA or protein sequence is used as a reference
                              Any non-gap character (e.g. x's) can indicate consensus/conserved/match columns
                              .'s or -'s indicate insert columns
                              ~'s indicate unaligned insertions
                              Upper and lower case can be used to discriminate strong and weakly conserved 
                              residues respectively

'''
RF = ""
assert len(ConsRow) == ColCount
assert len(IsMatch) == ColCount
for Col in range(ColCount):
	c = ConsRow[Col]
	if IsMatch[Col]:
		if not c.isalpha():
			c = 'x'
	else:
		c = '.'
	RF += c

MotifVec = ['_'] * ColCount
for Col in range(ColA, ColA+12):
	MotifVec[Col] = 'A'
for Col in range(ColB, ColB+14):
	MotifVec[Col] = 'B'
for Col in range(ColC, ColC+8):
	MotifVec[Col] = 'C'
MotifRow = ""
for c in MotifVec:
	MotifRow += c

fOut.write("# STOCKHOLM 1.0\n")

MatchA = ColToMatch[ColA]
MatchB = ColToMatch[ColB]
MatchC = ColToMatch[ColC]

assert MatchA >= 0
assert MatchB >= 0
assert MatchC >= 0

s = "#=CC PolHummer %s:%s:%d A:%d:%s B:%d:%s C:%d:%s" % \
  (Args.name, Args.order, MatchCount, MatchA + 1, ConsA, MatchB + 1, ConsB, MatchC + 1, ConsC)
fOut.write(s + "\n")

StartCol = 0
while 1:
	if StartCol >= ColCount:
		break
	fOut.write("\n")
	EndCol = StartCol + BLOCK - 1
	if EndCol >= ColCount:
		EndCol = ColCount - 1

	s = "%-*.*s  " % (W, W, "#=CC ABC")
	s += MotifRow[StartCol:EndCol+1]
	fOut.write(s + "\n")

	for SeqIndex in range(SeqCount):
		s = "%-*.*s  " % (W, W, Labels[SeqIndex])
		Row = Rows[SeqIndex]
		assert len(Row) == ColCount
		s += Row[StartCol:EndCol+1]
		fOut.write(s + "\n")
# Recommended placements: #=GF Above the alignment
# https://en.wikipedia.org/wiki/Stockholm_format#Recommended_features
	s = "%-*.*s  " % (W, W, "#=GC RF")
	s += RF[StartCol:EndCol+1]
	fOut.write(s + "\n")

	StartCol = EndCol + 1

fOut.write("//\n")
