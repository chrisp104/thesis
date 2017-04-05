import os

#	Use if we want antigen RESIDUE-specific contact information

# FUNCTIONS
# 1. resContacts
# 2. writeResContacts

# functions to return an array with residues within dist angstroms
# for each residue's beta carbon 

# *** function contacts()
#
# ARGUMENTS
# pdb: string - pdb file name of the structure in question
# ag: array - the chains of the Ag
# dist: num - angstrom threshold for determining 'contact'
#
# RETURNS 
# 	1. an array of numbers corresponding to the number of Fv residues in contact
#		with the nth residue in the Ag - (Ag res #, number of residues in contact)
#		2. an array of arrays, each array representing one Ag residue and containing
#		tuples with (Ag res #, res # of Ab contact residue, AA abbreviation, chain letter)
#		3. a dictionary with the chains and number of contact residues for each chain

def resContacts(pdb, ag, dist):
	struct = open(pdb, 'r')
	lines = struct.readlines()

	return_contact_info = []		# first return value
	return_contacts = []				# second return value
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

		# ## 1. returns with Ag res #s
		#	return_contact_info.append((agen[22:26].strip(), num_contacts))
		# return_contacts.append((agen[22:26].strip(), res_contacts))

		## 2. returns without Ag res #s
		return_contact_info.append((agen[21], agen[23:26], num_contacts))
		return_contacts.append((res_contacts))

	struct.close()
	# print	return_contact_info
	# print return_contacts
	return	return_contact_info, return_contacts, return_chain_num




# function to write and produce output of residue contact information into file
# ARUGMENTS
#	1. the output file name
#	2. the first return value from resContacts()
#	3. the second return value from resContacts()

def writeResContacts(file_name, numbers, residues, chain_nums):
	if not os.path.exists("./contact_output"):
		os.makedirs("./contact_output", 0777)
	output = open("./contact_output/"+file_name, 'w')
	for i in range(len(numbers)):
		# don't write anything if Ag residue has no contacts
		if len(residues[i]) == 0:
			continue

		# header information
		output.write("Ag Chain: "+str(numbers[i][0])+", "
			"res #: "+str(numbers[i][1])+"\n")
		output.write("# contacts: "+str(numbers[i][2])+"\n")
		output.write("Contacts in Ab: \n")
		
		# contact information
		for residue in residues[i]:
			output.write(residue[0]+' '+residue[1]+' '+residue[2]+"\n")
		output.write("\n")

	# chain contact numbers
	for key in chain_nums:
		output.write(key + ": " + str(chain_nums[key]) + "\n")




# bulkContacts()
# bulk writeResContacts for each docking model
#
# ARGUMENTS
# model_pre: string - file name prefix (before number) of second model pdb files
# chains: array containing characters of chains of Ag
#					e.g. ['O', 'R', 'T']
# dist: num - angstrom threshold for determining 'contact'
#
# RETURNS 
# 	outputs a file and prints the total contact numbers for each chain
def bulkContacts(model_pre, chains, dist):

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
		writeResContacts(name+"_contacts.txt", numbers, residues, chain_nums)

		file.close()












