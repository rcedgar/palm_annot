#!/usr/bin/env python3

Usage = \
(
"Convert tab-separated text name=value (fev) format"
" to tab-separated text with values only and optional"
" header with field names."
)

import sys
import argparse

AP = argparse.ArgumentParser(description = Usage)

AP.add_argument("--input",
  required=True,
  help="Input file in fev format")

AP.add_argument("--output",
  required=False,
  help="Output file in tsv format (default stdout)")

AP.add_argument("--fields",
  required=False,
  help="Comma-separated list of field names (default: all fields found in fev file")

AP.add_argument("--nullvalue",
  required=False,
  default="",
  help="String to use if value not specified (default empty string")

AP.add_argument("--header",
  required=False,
  choices=[ "no", "yes" ],
  default="yes",
  help="Include tsv header with field names as first line, yes or no (default yes)")

Args = AP.parse_args()

def GetNamesFromFile(FN):
	NameSet = set()
	for Line in open(FN):
		Fields = Line[:-1].split('\t')
		for Field in Fields[1:]:
			n = Field.find('=')
			if n > 0:
				Name = Field[:n]
				NameSet.add(Name)
	Names = []
	for Name in NameSet:
		Names.append(Name)
	return Names

if Args.fields != None:
	Names = Args.fields.split(',')
else:
	Names = GetNamesFromFile(Args.input)

fOut = sys.stdout
if not Args.output is None:
	fOut = open(Args.output, "w")

if Args.header == "yes":
	Hdr = "Label"
	for Name in Names:
		Hdr += "\t" + Name
	fOut.write(Hdr + "\n")

NameSet = set()
for Name in Names:
	NameSet.add(Name)

for Line in open(Args.input):
	Fields = Line[:-1].split('\t')
	Label = Fields[0]
	NameToValue = {}
	for Name in Names:
		NameToValue[Name] = Args.nullvalue
	for Field in Fields[1:]:
		n = Field.find('=')
		if n > 0:
			Name = Field[:n]
			Value = Field[n+1:]
			NameToValue[Name] = Value
	s = Label
	for Name in Names:
		s += "\t" + NameToValue[Name]
	fOut.write(s + "\n")
