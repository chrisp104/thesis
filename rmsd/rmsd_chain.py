# FUNCTIONS
# 1. rmsd
# 2. rmsdMultiple

# rmsd()
# calculates rmsd between two pdb files BY ANTIGEN CHAIN
# ****** CHAIN TO COMPARE MUST BE THE SAME LABEL!!

# ARGUMENTS
# f1: string - file name of first pdb
# f2: string - file name of second pdb
# chain: string - letter of chain to compare

# RETURNS 
# 	1. the rmsd between the designated chain between the two structures

def rmsd(f1, f2, chain):

	## FILES

	# first structure file *** USUALLY THE CRYSTAL STRUCTURE if exists
	file_1 = open(f1, 'r')
	struct_1 = file_1.readlines()
	
	# second structure file
	file_2 = open(f2, 'r')
	struct_2 = file_2.readlines()


	# *** CALCULATING RMSD *** #
	n = 0

	# loop through the lines of each from the start of the Ag for both files
	rmsd_squared = 0
	for c_line in struct_1:

		# skip if not the CB atom, CA and Gly, or BMET ** Comment out if wanna compare all atoms
		if (c_line[0:1] == "\n" or c_line[0:3] == "END"): break
		if not (c_line[21] == chain and c_line[0:4] == "ATOM"): continue
		if not ((c_line[17:20].strip() == "GLY" and 
			c_line[13:15].strip() == "CA") or 
			(c_line[13:15].strip() == "CB")):
			continue
		if c_line[16:20].strip() == "BMET":
			continue

		n += 1

		for m_line in struct_2:
			# make sure that we are looking at the right AA and find the right atom
			if (c_line[23:26] == m_line[23:26]):
				if (c_line[13:15] == m_line[13:15]):
					xsq = (float(c_line[30:38].strip()) - float(m_line[30:38].strip()))**2
					ysq = (float(c_line[38:46].strip()) - float(m_line[38:46].strip()))**2
					zsq = (float(c_line[46:54].strip()) - float(m_line[46:54].strip()))**2
					#print xsq, ysq, zsq
					squared = xsq + ysq + zsq
					rmsd_squared += squared
					break

	rmsd_s_mean = rmsd_squared / n
	rmsd = rmsd_s_mean**(0.5)
	#print f2 + ": " + str(rmsd)
	return rmsd

	# close files
	file_1.close()
	file_2.close()

# rmsd("5d1q_aligned.pdb", "D206m2d00.pdb", "E")


# rmsdMultiple()
# calculates rmsds between a target file and many model files against it

# ARGUMENTS
# target: string - file name of target pdb
# model_pre: string - file name prefix (before number) of second model pdb files
# chain: string - letter of chain to compare
# out: string - name of the output file to write all rmsds to

# RESULT 
# 	1. prints rmsds
#		2. outputs file with rmsds

def rmsdMultiple(target, model_pre, chain, out):

	output = open(out, 'w')

	# Calculate rmsd for each of the model files
	for i in range(40):
		if i < 10:
			i = "0" + str(i)
		name = model_pre + str(i) + ".pdb"
		try:
			file = open(name, 'r')
		except IOError as e:
			break
		result = rmsd(target, name, chain)
		output.write("Model #" + str(i)+'\n')
		output.write("rmsd = "+str(result)+'\n')
		file.close()
	output.close()




