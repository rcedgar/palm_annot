Line = None
Dict = None
Label = None
Fields = None

def WriteDict(File):
	Line = Label
	for Name in list(Dict.keys()):
		Value = Dict[Name]
		Line += "\t" + Name + "=" + Value
	File.write(Line)
	File.write("\n")

def WriteRec(File):
	WriteDict(File)

def DeleteField(Name):
	Dict.pop(Name, None)

def SetField(Name, Value):
	if Value is None:
		DeleteField(Name)
		return
	Dict[Name] = Value

def GetValue(Name):
	try:
		Value = Dict[Name]
	except:
		Value = None
	return Value

def GetStrValue(Name):
	return GetValue(Name)

def GetIntValue(Name):
	try:
		Value = Dict[Name]
	except:
		return None
	return int(Value)

def GetFloatValue(Name):
	try:
		Value = Dict[Name]
	except:
		return None
	return float(Value)

def ParseLine(argLine):
	global Label
	global Dict
	global Line
	global Fields

	if argLine.endswith('\n'):
		Line = argLine[:-1];
	else:
		Line = argLine

	Fields = Line.split('\t')
	Label = Fields[0]
	if len(Label) == 0:
		return False
	Found = False
	Dict = {}
	for Field in Fields[1:]:
		n = Field.find('=')
		if n > 0:
			Name = Field[:n]
			Value = Field[n+1:]
			Dict[Name] = Value
			Found = True
	return Found

def ReadRecsOnRec(FileName, OnRec):
	global Line, Dict, Label

	f = open(FileName)
	while 1:
		Line = f.readline()
		if len(Line) == 0:
			f.close()
			return
		Any = ParseLine(Line)
		if Any:
			OnRec()
