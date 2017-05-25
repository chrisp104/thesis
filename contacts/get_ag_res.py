import os

# script to take contact files for array of contact names and determine aggregated residue contact points

# 1.
# getEpitopeRegion()
#
# ARGUMENTS
# 1. output: str - path to write output file to
# 2. directory: str - path to directory with all model files
# 3. abs: array - containing Abs to aggregate contacts for, leave empty if want all
#
# RETURNS 
# 	outputs one file containing the Ag positions that are present across all 
# 	docking models
def getEpitopeRegion(output, directory, abs):
	
	# whether to do all Abs or just the ones in the array
	allAbs = False
	if len(abs) == 0:
		allAbs = True

	residues = []

	os.chdir(directory)
	for fn in os.listdir(directory):
		if allAbs and fn[0] != 'D':
			continue
		elif not allAbs and not fn[:9] in abs:
			continue


		# loop through file and determine which residues have at least one contact in Ab
		file = open(fn, 'r')
		for line in file.readlines():
			if line[:3] == "Ag#" or line[0] == "T" or line[0] == "O":
				continue

			if len(line) > 9 and not line[:3] in residues:
				residues.append(line[:3])

		file.close()

	out = open(output, 'w')
	for residue in sorted(residues):
		out.write(residue+"\n")
	out.close()




getEpitopeRegion("/Users/Chris/GitHub/thesis/contacts/D206/iter_6.txt",
	"/Users/Chris/GitHub/thesis/contacts/D206/n2",
	[
"D206n2d00",
"D206n2d17",
	])





