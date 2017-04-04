## rmsd.py
## calculates the rmsd of the ligand for each docking model to crystal s   cture
## FILES

# Crystal structure file
crystal_file = open('5d1z_aligned.pdb', 'r')
crystal = crystal_file.readlines()
output = open('rmsd_output.txt', 'w')

# Make dictionary to hold the model files
models = {}
for i in range(10):
	name = "D410m5d0" + str(i) + ".pdb"
	file = open(name, 'r')
	lines = file.readlines()
	models[str(i)] = lines
	file.close()

for i in range(10, 30):
	name = "D410m5d" + str(i) + ".pdb"
	file = open(name, 'r')
	lines = file.readlines()
	models[str(i)] = lines
	file.close()

# find the length of 5d1q chain E and line where it starts
c_start_found = False
c_start = 0
c_len = 0

for c_line in crystal:
	if c_line[0:3] == "": break
	if (c_line[21] == 'I' and c_line[0:4] == "ATOM"):

		# if first time seeing E then set the start line to this line
		if c_start_found == False:
			c_start_found = True

		c_len += 1

	if c_start_found == False:
		c_start += 1


for key in models:
	print key
	# find the length of the model ligand and line where it starts
	m_start_found = False
	m_start = 0
	m_len = 0

	model = models[key]
	for m_line in model:
		if m_line[0:3] == "": break
		if (m_line[21] == 'I' and m_line[0:4] == "ATOM"):

			if m_start_found == False:
				m_start_found = True

			m_len += 1

		if m_start_found == False:
			m_start += 1

	#print c_start, c_len, c_start_found, m_start, m_len, m_start_found


	# *** CALCULATING RMSD *** #
	output.write("Model #" + str(key)+'\n')
	n = 0

	# loop through the lines of each from the start of the Ag for both files
	rmsd_squared = 0
	for i in range(977):

		# skip if not the CB atom, CA and Gly, or BMET ** Comment out if wanna   mpare all atoms
		if not ((crystal[c_start+i][17:20].strip() == "GLY" and 
			crystal[c_start+i][13:15].strip() == "CA") or 
			(crystal[c_start+i][13:15].strip() == "CB")):
			continue
		if crystal[c_start+i][16:20].strip() == "BMET":
			continue

		#print crystal[c_start+i][17:20]
		
		n += 1
		for j in range(1190):
			# make sure that we are looking at the right AA and find the right atom
			if (crystal[c_start+i][23:26] == model[m_start+j][23:26]):
				if (crystal[c_start+i][13:15] == model[m_start+j][13:15]):
					xsq = (float(crystal[c_start+i][30:38].strip()) - float(model[m_start+j][30:38].strip()))**2
					ysq = (float(crystal[c_start+i][38:46].strip()) - float(model[m_start+j][38:46].strip()))**2
					zsq = (float(crystal[c_start+i][46:54].strip()) - float(model[m_start+j][46:54].strip()))**2
					#print xsq, ysq, zsq
					squared = xsq + ysq + zsq
					rmsd_squared += squared
					break

	rmsd_s_mean = rmsd_squared / n
	rmsd = rmsd_s_mean**(0.5)
	output.write("rmsd = "+str(rmsd)+'\n')

# close files
crystal_file.close()

## ** For printing different parts of the PDB File
# print crystal[c_start][0:4]
# print crystal[c_start][6:11]
# print crystal[c_start][13:16]
# print crystal[c_start][17:20]
# print crystal[c_start][21:22]
# print crystal[c_start][22:26]
# print crystal[c_start][30:38]
# print crystal[c_start][38:46]
# print crystal[c_start][46:54]
# print crystal[c_start][76:78]
# break
