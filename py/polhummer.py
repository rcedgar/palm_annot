#!/usr/bin/python3

import sys
import fasta

'''
123456789012345
#=CC PolHummer PH14:ABC:382 A:301:n.sDySRmDGti B:583:SGspdTSl.NTlrn C:758:ygGDDGl.

#=GC RF           axelfnaLVe.......tFNxYGDxmxxxraVRiedidiWItQLNPNrNSGNPDYTPVSKEQAvNd.YWPVmRe......
'''

def ReadSto(FN):
	global M, RF, MA, MB, MC, ABC, MatchCount
	ABC = None
	M = None
	MA = None
	MB = None
	MC = None
	RF = ""
	for Line in open(FN):
		if Line.startswith("#=CC PolHummer"):
			assert M == None
			Fields = Line[15:].split()
			assert len(Fields) == 4
			sPH = Fields[0]
			sA = Fields[1]
			sB = Fields[2]
			sC = Fields[3]

			Fields2 = sPH.split(':')
			assert len(Fields2) == 3
			PH = Fields2[0][0:2]
			assert PH == "PH" or PH == "PN"
			ABC = Fields2[1]
			M = int(Fields2[2])
		
			assert ABC == "ABC" or ABC == "CAB"
			assert sA.startswith("A:")
			assert sB.startswith("B:")
			assert sC.startswith("C:")

			MA = int(sA.split(':')[1])
			MB = int(sB.split(':')[1])
			MC = int(sC.split(':')[1])

		elif Line.startswith("#=GC RF "):
			RF += Line[:-1].split()[2]

	ColCount = len(RF)
	MA -= 1
	MB -= 1
	MC -= 1
	MatchToCol = []
	ColToMatch = []
	MatchIndex = 0
	for Col in range(ColCount):
		c = RF[Col]
		if c.isalpha():
			MatchToCol.append(Col)
			ColToMatch.append(MatchIndex)
			MatchIndex += 1
		elif c == '.':
			ColToMatch.append(-1)
		else:
			assert False

	MatchCount = len(MatchToCol)

# Get subsequence covering n match states
#   starting at the m'th match state 
def GetMSeq(Seq, m, n):
	L = len(Seq)
	MatchIndex = 0
	InsertCount = 0
	MSeq = ""
	Pos = None
	SeqPos = 0
	for i in range(L):
		c = Seq[i]
		if c.islower():
			if MatchIndex > m and MatchIndex < m+n:
				InsertCount += 1
				MSeq += c.upper()
		elif c.isupper() or c == '-':
			if c != '-':
				SeqPos += 1
			if MatchIndex >= m and MatchIndex < m+n:
				if Pos == None:
					Pos = SeqPos
				MSeq += c
			MatchIndex += 1
		elif c == '.':
			pass
		else:
			assert False
	if MatchIndex != MatchCount:
		sys.stderr.write("MI=%d, MC=%d\n" % (MatchIndex, MatchCount))
		assert False
	if len(MSeq) != n + InsertCount:
		s = "\n"
		s += "RF =" + RF + "\n"
		s += "Seq=" + Seq + "\n"
		s += "len(MSeq)=%d n=%d InsertCount=%d\n" % (len(MSeq), n, InsertCount)
		s += "\n"
		sys.stderr.write(s)
		assert False

	UM = MSeq.replace('-', '')
	if len(UM) == 0:
		Pos = -1
	MSeq = MSeq[:n]
	return Pos, MSeq

def GetMotifs(Seq):
	PosA, SeqA = GetMSeq(Seq, MA, 12)
	PosB, SeqB = GetMSeq(Seq, MB, 14)
	PosC, SeqC = GetMSeq(Seq, MC, 8)

	return SeqA, SeqB, SeqC
