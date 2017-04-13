import os

#	Use if we want Ab contacts info for a given Ag as well
#
# FUNCTIONS
# 1. resContacts
# 2. writeResContacts
# 3. bulkContacts
# 4. bulkDirectory
# 5. percentContact
#
# functions to return an array with residues within dist angstroms
# for each residue's beta carbon 


# 1.
# contacts()
#
# ARGUMENTS
# pdb: string - pdb file name of the structure in question
# ag: array - the chains of the Ag
# dist: num - angstrom threshold for determining 'contact'
# res_num_only: boolean - set to true if key in return_contacts should just be residue number
#
# RETURNS 
# 	1. an array of tuples corresponding to the number of Fv residues in contact
#		with the nth residue in the Ag - (Ag chain, Ag res #, number of residues in contact)
#		2. an dictionary of arrays, key = Ag res # and key =  array of
#		tuples with (Ab contact residue #, AA abbreviation, chain letter)
#		3. a dictionary with the chains and number of contact residues for each chain
def resContacts(pdb, ag, dist, res_num_only=False):
	struct = open(pdb, 'r')
	lines = struct.readlines()

	return_contact_info = []		# first return value
	return_contacts = {}				# second return value
	return_chain_num = {}				# third return value

	# instantiate the desired chains
	for chain in ag:
		return_chain_num[chain] = 0

	# loop through the Ag's residues and get contacts for each C beta
	for agen in lines:

		num_contacts = 0		# the number of Ab residues this Ag residue is contacting
		res_contacts = []		# to hold the Ab contact residues for second return value

		# skip if not the CB atom, CA and Gly, or BMET ** Comment out if wanna commpare all atoms
		if not (agen[21] in ag):
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
			if (abody[16:20].strip() == "BMET" or abody[21] in ag):
				continue

			# else calculate the rmsd between the Ag residue and Ab residue
			xsq = (float(agen[30:38].strip()) - float(abody[30:38].strip()))**2
			ysq = (float(agen[38:46].strip()) - float(abody[38:46].strip()))**2
			zsq = (float(agen[46:54].strip()) - float(abody[46:54].strip()))**2
			squared = xsq + ysq + zsq
			distance = squared**(0.5)

			# count if less than $dist anstroms
			if (distance <= dist):
				num_contacts += 1

				# add to second return array's array
				ab_res = (abody[22:26].strip(), abody[17:20], abody[21:22])
				res_contacts.append(ab_res)

				# increment the chain's contact count
				chain = agen[21:22]
				return_chain_num[chain] += 1  


		## 2. Add information to return values for this Ag res
		return_contact_info.append((agen[21], agen[23:26], num_contacts))
		if not res_num_only:
			return_contacts[(agen[23:26], agen[17:20])] = (res_contacts)
		else:
			return_contacts[agen[23:26]] = (res_contacts)

	struct.close()
	print pdb + " finished."
	return	return_contact_info, return_contacts, return_chain_num




# 2. 
# writeResContacts()
#
# function to write and produce output of residue contact information into file
#
# ARUGMENTS
#	1. path: str - path to directory to write out files to
# 2. file_name: str - output file name
#	3. numbers: array - the first return value from resContacts()
#	4. residues: dictionary - the second return value from resContacts()
# 5. chain_nums: dictionary - the third
def writeResContacts(path, file_name, numbers, residues, chain_nums):
	output = open(path+file_name, 'w')
	for key in sorted(residues):

		# header information
		output.write("Ag res: "+str(key)+"\n")
		output.write("Contacts in Ab: \n")
		contacts = residues[key]
		
		# contact information
		for contact in contacts:
			output.write(contact[0]+' '+contact[1]+' '+contact[2]+"\n")
		output.write("\n")

	# chain contact numbers
	for key in chain_nums:
		output.write(key + ": " + str(chain_nums[key]) + "\n")




# 3. 
# bulkContacts()
#
# bulk writeResContacts for each docking model 
#
# ARGUMENTS
# 1. model_pre: string - file name prefix (before number) of second model pdb files
# 2. chains: array containing characters of chains of Ag
#					e.g. ['O', 'R', 'T']
# 3. dist: num - angstrom threshold for determining 'contact'
# 4. path: str - the directory to write the output files to
#
# RETURNS 
# 	outputs a file and prints the total contact numbers for each chain
def bulkContacts(model_pre, chains, dist, path):

	# Find contacts for each Ag chain for each of the model files
	for i in range(40):
		if i < 10:
			i = "0" + str(i)
		name = model_pre + str(i)
		pdb_name = name+".pdb"
		try:
			file = open(pdb_name, 'r')
		except IOError as e:
			break

		numbers, residues, chain_nums = resContacts(pdb_name, chains, dist)
		writeResContacts(path, name+"_contacts.txt", numbers, residues, chain_nums)

		file.close()





# 4. 
# bulkDirectory()
#
# bulk writeResContacts for each docking model in a given directory
#
# ARGUMENTS
# 1. chains: array containing characters of chains of Ag
#					e.g. ['O', 'R', 'T']
# 2. dist: num - angstrom threshold for determining 'contact'
# 3. path: str - the directory to write the output files to
#
# RETURNS 
# 	outputs a file and prints the total contact numbers for each chain
def bulkDirectory(chains, dist, path):

	# Find contacts for each Ag chain for each of the model files
	for fn in os.listdir('.'):
		if fn[0:1] == '.': continue
		try:
			file = open(fn, 'r')
		except IOError as e:
			break

		numbers, residues, chain_nums = resContacts(fn, chains, dist)
		writeResContacts(path, fn[:-7]+"_contacts.txt", numbers, residues, chain_nums)

		file.close()





# 5.
# percentContact()
#
# take two outputs of residue contacts (second return value of res_contacts first function)
# 1 of the correct contacts derived from the crystal structure and 1 from the mutants or 
# models to be compared (is a directory of models)
#
# PERCENT CORRECT CONTACT is calculated as number of contacts in non crystal
# structure docking model that are correctly within $dist angstroms / total number of contacts
# as determined by crystal structure
#
# INCORECT " is same idea???
#
# ARGUMENTS
# 1. crystal: str - crystal structure pdb
# 2. chains: array containing characters of chains of Ag
#					e.g. ['O', 'R', 'T']
# 3. dist: num - angstrom threshold for determining 'contact'
# 4. output: str - path to write outpu file to
#
# RETURNS 
# 	outputs one file containing percent correct contact for each non-crystal model
def percentContact(crystal, chains, dist, output):
	out = open(output, 'w')

	# GET CRYSTAL CORRECT CONTACTS
	c_numbers, correct_contacts, c_chain_nums = resContacts(crystal, chains, dist, True)

	# Find contacts for each Ag chain for each of the model files
	for fn in os.listdir('.'):
		print fn
		if fn[0:1] == '.' or fn[-3:] != "pdb": continue
		try:
			file = open(fn, 'r')
		except IOError as e:
			break

		numbers, residues, chain_nums = resContacts(fn, chains, dist, True)
		file.close()


		# ** COMPARE against the correct contacts determined from crystal structure
		num_correct = 0
		total = 0
		for resi in sorted(correct_contacts):

			c_con = correct_contacts[resi]
			m_con = residues[resi]

			# loop through the correct Ab contacts in the crystal contact array
			for contact in c_con:
				total += 1
				if contact in m_con:
					num_correct += 1

		percent_correct = round(float(num_correct)/total, 3)

		print num_correct, total
		print str(percent_correct) + "\n"

		out.write(fn+": "+ str(percent_correct*100)+"% correct\n\n")













