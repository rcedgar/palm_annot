d = {}
Order = []
Keys = []

def VectToDict(v):
	d = {}
	for x in v:
		try:
			d[x] += 1
		except:
			d[x] = 1
	return d

def CmpKey__(i):
	global d, Keys
	ki = Keys[i]
	ni = d[ki]
	return ni

def CmpKeyVec__(i):
	global Keys
	ki = Keys[i]
	return ki

def GetOrder(Dict):
	global d, Order, Keys

	d = Dict
	Keys = list(d.keys())
	N = len(Keys)
	Order = list(range(0, N))
	Order.sort(key=CmpKey__)
	Order.reverse()
	return Order

def IncCount(Dict, Key, n=1):
	try:
		Dict[Key] += n
	except:
		Dict[Key] = n
	return Dict[Key]

def IncCount2(Dict, Key1, Key2, n=1):
	try:
		x = Dict[Key1]
	except:
		Dict[Key1] = {}
	IncCount(Dict[Key1], Key2, n)

def GetCount(Dict, Key):
	try:
		return Dict[Key]
	except:
		return 0

def GetCount2(Dict, Key1, Key2):
	try:
		return Dict[Key1][Key2]
	except:
		return 0

def GetRanks(Dict):
	Order = GetOrder(Dict)
	Keys = list(Dict.keys())
	KeyToRank = {}
	N = len(Order)
	for Rank in range(0, N):
		i = Order[Rank]
		Key = Keys[i]
		KeyToRank[Key] = Rank
	return KeyToRank

def CmpVec__(i, j):
	global g_Vec
	if g_Vec[i] < g_Vec[j]:
		return 1
	elif g_Vec[i] > g_Vec[j]:
		return -1
	return 0

def GetOrderDesc(Vec):
	global Keys
	Keys = Vec
	Order = list(range(0, len(Vec)))
	Order.sort(key=CmpKeyVec__)
	Order.reverse()
	return Order

def GetOrderAsc(Vec):
	global Keys
	Keys = Vec
	Order = list(range(0, len(Vec)))
	Order.sort(key=CmpKeyVec__)
	return Order

def Test():
	v = [91, 27884, 25365, 3]
	OrderA = GetOrderAsc(v)
	OrderD = GetOrderDesc(v)
	print(v)
	print(("A:", OrderA))
	print(("D:", OrderD))

if __name__ == "__main__":
	Test()
