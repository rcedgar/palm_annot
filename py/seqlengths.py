#!/usr/bin/env python3

import sys
import fasta

FastaFileName = sys.argv[1]

def OnSeq(Label, Seq):
	print("%s	seqlength=%d" % (Label, len(Seq)))

fasta.ReadSeqsOnSeq(FastaFileName, OnSeq)
