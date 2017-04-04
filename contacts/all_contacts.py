#	Use if we want all contact information that is NOT Ag residue specific

# FUNCTIONS
# 1. allContacts
# 2. writeAllContacts

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
# 	1. an array of arrays - each entry is [Ab residue #, # contacts]

def allContacts(pdb, ag, dist):
	struct = open(pdb, 'r')
	lines = struct.readlines()

	contacts = []		# return value

	# loop through the Ag's residues and get contacts for each C beta
	for agen in lines:

		num_contacts = 0		# the number of Ab residues this Ag residue is contacting

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

		# ******* STRIP if you don't want leading white spaces
		ag_number = agen[23:26]
		contacts.append([ag_number, num_contacts])

	struct.close()
	# print	return_contact_info
	# print return_contacts
	return	contacts



# function to write and produce output of residue contact information into file
# ARUGMENTS
#	1. the output file name
#	2. the return value from allContacts()

def writeAllContacts(file_name, numbers):
	output = open(file_name, 'w')

	for i in range(len(numbers)):
		# don't write anything if Ag residue has no contacts
		if numbers[i][1] == 0:
			continue

		# otherwise display the Ag res number and number of contacts
		output.write("Res #" + numbers[i][0]+": "+str(numbers[i][1])+"\n")


numbers = allContacts("example.pdb", ['O', 'R', 'T'], 8)
writeAllContacts("all_num_out", numbers)













