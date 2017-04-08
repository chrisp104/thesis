# FUNCTIONS
# 1. find_region()
# 2. find_bin()
#

# 1.
# find_region()
# for the three chains in the Ag, find the number of contacts for each
# peptide and cumulate
#
# ARGUMENTS
# file: string - file name of first pdb
# chains: array containing characters of chains of Ag
#					e.g. ['O', 'R', 'T']
# dist: num - angstrom threshold for determining 'contact'
#
# RETURNS 
# 	1. the number of contacts for each chain in the 

def find_region(file, chains, dist):
	struct = open(file, 'r')
	lines = struct.readlines()

	# create the return dictionary
	chain_nums = {}
	for key in chains:
		chain_nums[key] = 0

	# loop through the Ag's residues and get contacts for each C beta
	for agen in lines:

		# skip if not the CB atom, CA and Gly, or BMET ** Comment out if wanna commpare all atoms
		if not agen[21] in chains:
			continue
		if not ((agen[17:20].strip() == "GLY" and agen[13:15].strip() == "CA") or 
			(agen[13:15].strip() == "CB")):
			continue
		if agen[16:20].strip() == "BMET":
			continue

		# loop through Ab residues (!Ag residues)
		for abody in lines:
			# skip if not the CB atom, CA and Gly, or if part of Ag
			if not ((abody[17:20].strip() == "GLY" and abody[13:15].strip() == "CA") or 
				(abody[13:15].strip() == "CB")):
				continue
			if (abody[16:20].strip() == "BMET" or abody[21] in chains):
				continue

			# else calculate the rmsd between the Ag residue and Ab residue
			xsq = (float(agen[30:38].strip()) - float(abody[30:38].strip()))**2
			ysq = (float(agen[38:46].strip()) - float(abody[38:46].strip()))**2
			zsq = (float(agen[46:54].strip()) - float(abody[46:54].strip()))**2
			squared = xsq + ysq + zsq
			distance = squared**(0.5)

			# count if less than $dist anstroms
			if (distance <= dist):
				chain = agen[21:22]
				chain_nums[chain] += 1

	struct.close()
	print chain_nums
	return chain_nums


## ******** TEST CASE **********
#find_region("example.pdb", ['O', 'R', 'T'], 8)


# 2.
# find_bin()
# bulk find_region for each docking model and find most contacted Ag chain
#
# ARGUMENTS
# model_pre: string - file name prefix (before number) of second model pdb files
# chains: array containing characters of chains of Ag
#					e.g. ['O', 'R', 'T']
# dist: num - angstrom threshold for determining 'contact'
# out: output text file name
#
# RETURNS 
# 	outputs a file and prints the total contact numbers for each chain
def find_bin(model_pre, chains, dist, out):
	
	output = open(out, 'w')
	totals = {}
	for chain in chains:
		totals[chain] = 0

	# Calculate contact counts for each Ag chain for each of the model files
	for i in range(40):
		if i < 10:
			i = "0" + str(i)
		name = model_pre + str(i) + ".pdb"
		try:
			file = open(name, 'r')
		except IOError as e:
			break
		result = find_region(name, chains, dist)
		output.write("Model #" + str(i)+'\n')
		output.write("Result: "+str(result)+'\n')

		# add all model results to find chain contact totals
		for ag_chain in result:
			totals[ag_chain] += result[ag_chain]

		file.close()

	# write the chain contact totals to the file
	output.write("TOTALS \n")
	for key in totals:
		output.write(key +": "+str(totals[key])+"\n")

	output.close()







