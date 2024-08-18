#!/usr/bin/python3

# input is fasta with A:n:XXXXXX ... annotations.
# output those with intact palmprint (no gaps in any motif)

import sys
import fasta

Invert = False

FileNames = sys.argv[1:]
if len(FileNames) == 0:
	FileNames = [ "/dev/stdin" ]

def OnSeq(Label, Seq):
	Fields = Label.split()
	A = ""
	B = ""
	C = ""
	for Field in Fields:
		if Field.startswith("A:"):
			A = Field.split(':')[2].replace("-", "")
		elif Field.startswith("B:"):
			B = Field.split(':')[2].replace("-", "")
		elif Field.startswith("C:"):
			C = Field.split(':')[2].replace("-", "")

	Discard = len(A) < 10 or len(B) < 12 or len(C) < 6
	if Invert:
		Discard = not Discard
	if Discard:
		return
	fasta.WriteSeq(sys.stdout, Seq, Label)

for FaFN in FileNames:
	if FaFN == "--invert":
		Invert = True
		continue
	sys.stderr.write(FaFN + "\n")
	fasta.ReadSeqsOnSeq(FaFN, OnSeq)
