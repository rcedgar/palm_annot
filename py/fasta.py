def ReadSeqsOnSeq(FileName, OnSeq, trunclabels=False):
	Label = None
	Seq = ""
	for Line in open(FileName):
		Line = Line.strip()
		if len(Line) == 0:
			continue
		if Line[0] == ">":
			if len(Seq) > 0:
				if trunclabels:
					Label = Label.split()[0]
				OnSeq(Label, Seq)
			Label = Line[1:]
			Seq = ""
		else:
			Seq += Line.replace(" ", "")
	if len(Seq) > 0:
		OnSeq(Label, Seq)

def ReadSeqsDict(FileName, trunclabels=False):
	SeqDict = {}
	Label = None
	Seq = ""
	for Line in open(FileName):
		Line = Line.strip()
		if len(Line) == 0:
			continue
		if Line[0] == ">":
			if len(Seq) > 0:
				SeqDict[Label] = Seq
			Label = Line[1:]
			if trunclabels:
				Label = Label.split()[0]
			Seq = ""
		else:
			Seq += Line.replace(" ", "")
	if len(Seq) > 0:
		SeqDict[Label] = Seq
	return SeqDict

def WriteSeq(File, Seq, Label, BLOCKLENGTH = 80):
	if len(Seq) == 0:
		return
	if Label != "":
		File.write(">" + Label + "\n")
	if BLOCKLENGTH <= 0:
		File.write(Seq + "\n")
	else:
		SeqLength = len(Seq)
		BlockCount = (SeqLength + (BLOCKLENGTH-1))//BLOCKLENGTH
		for BlockIndex in range(BlockCount):
			Block = Seq[BlockIndex*BLOCKLENGTH:]
			Block = Block[:BLOCKLENGTH]
			File.write(Block + "\n")
