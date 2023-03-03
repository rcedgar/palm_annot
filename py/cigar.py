import sys

# MIDNSHPX=

# M 0 alignment match (can be a sequence match or mismatch)
# I 1 insertion to the reference
# D 2 deletion from the reference
# N 3 skipped region from the reference
# S 4 soft clipping (clipped sequences present in SEQ)
# H 5 hard clipping (clipped sequences NOT present in SEQ)
# P 6 padding (silent deletion from padded reference)
# = 7 sequence match
# X 8 sequence mismatch

Ns = []
Letters = []

def ParseCigar(C):
	global Ns, Letters

	Ns = []
	Letters = []

	if C == "*":
		return

	n = len(C)
	assert n > 0
	if not C[0].isdigit():
		print(("not C[0].isdigit()", C))
		assert False
	assert C[n-1].isalpha()

	Ns = []
	Letters = []

	N = 0
	for c in C:
		if c.isdigit():
			N = N*10 + (ord(c) - ord('0'))
		elif c.isalpha() or c == '=':
			Letters.append(c)
			Ns.append(N)
			N = 0
		else:
			print(("not letter or digit", C))
			assert False
	return Ns, Letters

def GetLeftClip():
	if len(Letters) == 0:
		return 0
	if Letters[0] != 'S':
		return 0
	return Ns[0]

def GetRightClip():
	if len(Letters) < 2:
		return 0
	if Letters[-1] != 'S':
		return 0
	return Ns[-1]
	
def GetAlnLength():
	global Ns, Letters

	L = 0
	for i in range(0, len(Ns)):
		x = Letters[i]
		if x != 'S' and x != 'H':
			L += Ns[i]
	return L

def GetReadSegLength():
	global Ns, Letters

	L = 0
	for i in range(0, len(Ns)):
		x = Letters[i]
		if x == 'M' or x == 'D' or x == '=' or x == 'X':
			L += Ns[i]
	return L

def GetReadLength():
	global Ns, Letters

	L = 0
	for i in range(0, len(Ns)):
		x = Letters[i]
		if x == 'M' or x == 'I' or x == '=' or x == 'X' or x == 'S':
			L += Ns[i]
	return L

def GetRefSegLength():
	global Ns, Letters

	L = 0
	for i in range(0, len(Ns)):
		x = Letters[i]
		if x == 'M' or x == 'D' or x == '=' or x == 'X':
			L += Ns[i]
	return L

def GetSoftClipCount():
	L = 0
	for i in range(0, len(Ns)):
		x = Letters[i]
		if x == 'S':
			L += Ns[i]
	return L

def GetReadRow(SEQ):
	Row = ""
	Pos = 0
	C = len(Ns)
	for i in range(0, C):
		c = Letters[i]
		n = Ns[i]
		if c == 'M' or c == 'I':
			for k in range(0, n):
				Row += SEQ[Pos]
				Pos += 1
		elif c == 'D':
			for k in range(0, n):
				Row += '-'
		elif c == 'S':
			Pos += n
		else:
			assert False
	return Row

def GetTargetRow(T):
	Row = ""
	C = len(Ns)
	Pos = 0
	for i in range(0, C):
		c = Letters[i]
		n = Ns[i]
		if c == 'M' or c == 'D':
			for k in range(0, n):
				Row += T[Pos]
				Pos += 1
		elif c == 'I':
			for k in range(0, n):
				Row += '-'
		elif c == 'S':
			pass
		else:
			assert False
	return Row
