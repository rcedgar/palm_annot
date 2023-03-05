#!/usr/bin/env python3

import sys
import fasta

FastaFileName = sys.argv[1]

def OnSeq(Label, Seq):
	Label0 = Label.split()[0]
	print("%s	seqlength=%d" % (Label0, len(Seq)))

fasta.ReadSeqsOnSeq(FastaFileName, OnSeq)
