import os

# ******* use this when chain letters MUST be the same as well - for finding
# rmsd between models of the same Ag

# FUNCTIONS
# 1. rmsdPair
# 2. rmsdPairMultiple

# 1.
# rmsd()
# calculates rmsd between two pdb files, the first the crystal structure
# with the specific epitope and the second the docked Ab:Ag structure with 
# the whole Ag
#
# ****** CHAINS TO COMPARE DON'T HAVE TO BE SAME LABEL BUT MUST 
# ****** BE SAME AMINO ACID NUMBER!!
#
# ARGUMENTS
# f1: string - file name of first pdb
# f2: string - file name of second pdb
# chains: array - containing Ag chain letters in crystal structure
# NOTE: chain and ag should be same if just comparing two homologous structures
#
# RETURNS 
# 	1. the rmsd between the designated chain between the two structures
def rmsdPair(f1, f2, chains):

	## FILES

	# first structure file *** USUALLY THE CRYSTAL STRUCTURE if exists
	file_1 = open(f1, 'r')
	struct_1 = file_1.readlines()
	
	# second structure file
	file_2 = open(f2, 'r')
	struct_2 = file_2.readlines()
	

	# *** CALCULATING RMSD *** #
	n = 0
	count = 0		# for checking how many residues were compared

	# loop through the lines of each from the start of the Ag for both files
	rmsd_squared = 0
	for c_line in struct_1:

		# skip if not the CB atom, CA and Gly, or BMET ** Comment out if wanna compare all atoms
		if (c_line[0:3] == "TER"): continue
		if (c_line[0:1] == "\n" or c_line[0:3] == "END" or c_line[0:8] == "CDR_SASA"): break
		if not (c_line[21] in chains and c_line[0:4] == "ATOM"): continue
		if not ((c_line[17:20].strip() == "GLY" and c_line[13:15].strip() == "CA") or 
			(c_line[13:15].strip() == "CB")): continue
		if c_line[16:20].strip() == "BMET": continue

		for m_line in struct_2:
			# make sure that we are looking at the right AA in an Ag chain and find the right atom
			##********************** CHECK THIS
			if (m_line[0:3] == "TER"): continue
			if (m_line[0:1] == "\n" or m_line[0:3] == "END" or m_line[0:8] == "CDR_SASA"): break
			if not (m_line[21] in chains and m_line[0:4] == "ATOM"): continue
			if not ((m_line[17:20].strip() == "GLY" and m_line[13:15].strip() == "CA") or 
				(m_line[13:15].strip() == "CB")): continue
			if m_line[16:20].strip() == "BMET": continue

			if (c_line[23:27].strip() == m_line[23:27].strip() and c_line[17:20] == m_line[17:20]
				and c_line[21] == m_line[21]):
				# print c_line
				# print m_line

				xsq = (float(c_line[30:38].strip()) - float(m_line[30:38].strip()))**2
				ysq = (float(c_line[38:46].strip()) - float(m_line[38:46].strip()))**2
				zsq = (float(c_line[46:54].strip()) - float(m_line[46:54].strip()))**2
			
				squared = xsq + ysq + zsq
				rmsd_squared += squared
				# print squared
				# print rmsd_squared
				n += 1
				break

	rmsd_s_mean = rmsd_squared / n
	rmsd = rmsd_s_mean**(0.5)
	print f2 + ": " + str(rmsd)

	return rmsd

	# close files
	file_1.close()
	file_2.close()






# 2.
# rmsdPairMultiple()
# calculates rmsds between 10 models in a directory
#
# ARGUMENTS
# directory: string - directory with pdb files to pairwise compare
# chains: array - containing chain letters for the chains to compare
# out: string - path of the output file to write all rmsds to
#
# RESULT 
# 	1. prints rmsds
#		2. outputs file with rmsds
def rmsdPairMultiple(directory, chains, out):
	print directory

	output = open(out, 'w')
	rmsd_dict = {}
	
	for i in range(1, 11):
		for j in range(i+1, 11):

			f1 = directory+"model.0_"+str(i).zfill(4)+".pdb"
			f2 = directory+"model.0_"+str(j).zfill(4)+".pdb"

			try:
				file1 = open(f1, 'r')
			except IOError as e:
				continue

			try:
				file2 = open(f2, 'r')
			except IOError as e:
				continue

			result = rmsdPair(f1, f2, chains)

			rmsd_dict[result] = str(i)+ " & "+str(j)

	for key in sorted(rmsd_dict):
		output.write(rmsd_dict[key]+": "+str(key)+"\n\n")
	output.close()

