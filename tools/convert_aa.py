peptide = open('peptide.txt', 'r')
output = ''
for aa in peptide:
	if (aa.lower() == "ala"): output += "A"
	if (aa.lower() == "gly"): output += "G"
	if (aa.lower() == "ile"): output += "I"
	if (aa.lower() == "leu"): output += "L"
	if (aa.lower() == "pro"): output += "P"
	if (aa.lower() == "val"): output += "V"
	if (aa.lower() == "phe"): output += "F"
	if (aa.lower() == "trp"): output += "W"
	if (aa.lower() == "tyr"): output += "Y"
	if (aa.lower() == "asp"): output += "D"
	if (aa.lower() == "glu"): output += "E"
	if (aa.lower() == "arg"): output += "R"
	if (aa.lower() == "his"): output += "H"
	if (aa.lower() == "lys"): output += "K"
	if (aa.lower() == "ser"): output += "S"
	if (aa.lower() == "thr"): output += "T"
	if (aa.lower() == "cys"): output += "C"
	if (aa.lower() == "met"): output += "M"
	if (aa.lower() == "asn"): output += "N"
	if (aa.lower() == "gln"): output += "Q"

print output