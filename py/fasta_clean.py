#!/usr/bin/env python3

Usage = \
(
"Clean FASTA file by removing white space from labels and eliminating duplicate labels."
)

import sys
import argparse
import fasta

AP = argparse.ArgumentParser(description = Usage)

AP.add_argument("--input",
  required=False,
  help="Input file in FASTA format (default stdin)")

AP.add_argument("--output",
  required=False,
  help="Output file in FASTA format (default stdout)")

AP.add_argument("--white",
  required=False,
  choices=[ "truncate", "replace" ],
  default="truncate",
  help="Remove label white space by truncating or replacing (default truncate)")

AP.add_argument("--whitestr",
  required=False,
  default="_",
  help="Replacement string for --white replace (default '_')")

AP.add_argument("--dupes",
  required=False,
  choices=[ "delete", "relabel" ],
  default="delete",
  help="Eliminate duplicate labels by deleting sequence or relabeling (default delete)")

AP.add_argument("--pattern",
  required=False,
  default="/dupe@",
  help="Pattern to append to duplicate label, must contain @ which is replaced by 1,2... (default '/dupe@')")

Args = AP.parse_args()

InputName = Args.input
if Args.input is None:
	InputFileName = "/dev/stdin"
else:
	InputFileName = Args.input

for c in Args.whitestr:
	if c.isspace():
		sys.stderr.write("\n\nERROR: White space in --whitestr\n\n")
		sys.exit(1)

Pattern = Args.pattern
if Pattern.find("%") >= 0:
	sys.stderr.write("\n\nERROR: % not allowed in pattern\n\n")
	sys.exit(1)

Pattern = Args.pattern.replace("@", "%d")
if Pattern.find("%d") < 0:
	sys.stderr.write("\n\nERROR: Must be @ in --pattern\n\n")
	sys.exit(1)

fOut = sys.stdout
if not Args.output is None:
	fOut = open(Args.output, "w")

SeqCount = 0
LabelsWithWhiteSpaceCount = 0
DupeLabelsCount = 0

def FixWhite(Label):
	global LabelsWithWhiteSpaceCount
	FixedLabel = ""
	Changes = False
	for c in Label:
		if c.isspace():
			if Args.white == "truncate":
				Changes = True
				break
			c = Args.whitestr
			Changes = True
		FixedLabel += c
	if Changes:
		LabelsWithWhiteSpaceCount += 1
	return FixedLabel

LabelSet = set()
def OnSeq(Label, Seq):
	global SeqCount
	global DupeLabelsCount

	SeqCount += 1
	FixedLabel = FixWhite(Label)
	if FixedLabel == "":
		return
	if FixedLabel in LabelSet:
		DupeLabelsCount += 1
		if Args.dupes == "delete":
			return
		for n in range(1, 10000):
			if n == 9999:
				sys.stderr.write("\n\nERROR too many duplicate labels >%s\n\n" % Label)
				sys.exit(1)
			Suffix = Pattern % n
			Labeln = FixedLabel + Suffix
			if not Labeln in LabelSet:
				FixedLabel = Labeln
				break
	assert FixedLabel not in LabelSet
	LabelSet.add(FixedLabel)
	fasta.WriteSeq(fOut, Seq, FixedLabel)

fasta.ReadSeqsOnSeq(InputFileName, OnSeq)

sys.stderr.write("%d seqs, %d labels with white space, %d dupe labels\n" % \
  (SeqCount, LabelsWithWhiteSpaceCount, DupeLabelsCount))
