#!/usr/bin/python3

# input is fasta with A:n:XXXXXX ... annotations.
# output those with intact palmcore 150pp150

import sys
import fasta

FLANK = 150
MAXL = 800

FileNames = sys.argv[1:]
if len(FileNames) == 0:
	FileNames = [ "/dev/stdin" ]

def OnSeq(Label, Seq):
	Fields = Label.split()

	A = None
	B = None
	C = None
	Fields = Label.split()
	for Field in Fields:
		if Field.startswith("A:"):
			A = Field.split(':')[2]
		elif Field.startswith("B:"):
			B = Field.split(':')[2]
		elif Field.startswith("C:"):
			C = Field.split(':')[2]
	if A is None or B is None or C is None:
		return

	Seq = Seq.upper()
	PosA = Seq.find(A)
	PosB = Seq.find(B)
	PosC = Seq.find(C)
	if PosA < 0 or PosC < 0 or PosB < 0:
		return
	ABC = (PosA < PosB and PosB < PosC)
	if not ABC:
		return
	PCStart = PosA - FLANK
	if PCStart < 0:
		return
	if PosC > len(Seq) + FLANK:
		return
	TrimSeq = Seq[PCStart:PosC+FLANK]
	if len(TrimSeq) > MAXL:
		return
	fasta.WriteSeq(sys.stdout, TrimSeq, Label)

for FaFN in FileNames:
	if FaFN.startswith("-FLANK"):
		FLANK = int(FaFN[2:])
		sys.stderr.write("FLANK = %d\n" % FLANK)
		continue
	sys.stderr.write(FaFN + "\n")
	fasta.ReadSeqsOnSeq(FaFN, OnSeq)
