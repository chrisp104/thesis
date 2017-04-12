import os

#	Use if we want all contact information that is NOT Ab residue specific
# for each Ag residue
#
# FUNCTIONS
# 1. allContacts
# 2. writeAllContacts
# 3. bulkAllContacts
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
# binary: boolean - if true, will just set contact number as 1, else count
#
# RETURNS 
# 	1. an array of arrays - each entry is [Ab residue #, # contacts]

def allContacts(pdb, ag, dist, binary=True):
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

		# good for if we don't need to normalize and want absolute contact numbers
		# across all docks
		if not binary:
			contacts.append([ag_number, num_contacts])

		# if we just want num_contacts to be 1.
		if binary:
		 	if num_contacts > 0:
				contacts.append([ag_number, 1])
			else:
				contacts.append([ag_number, 0])

	struct.close()
	return contacts




# 2.
#	writeAllContacts()
#
# function to write and produce output of residue contact information into file
#
# ARUGMENTS
#	1. path: str - path to directory to write out files to
# 2. file_name: str - output file name
# 3. numbers: array - the return value from allContacts()
def writeAllContacts(path, file_name, numbers):
	output = open(path+file_name, 'w')

	for i in range(len(numbers)):
		# don't write anything if Ag residue has no contacts
		if numbers[i][1] == 0:
			continue

		# otherwise display the Ag res number and number of contacts
		output.write("Res #" + numbers[i][0]+": "+str(numbers[i][1])+"\n")




# 3.
# bulkAllContacts()
#
# bulk writeAllContacts for each docking model
#
# ARGUMENTS
# 1. path: str - the directory to write the output files to
# 2. model_pre: string - file name prefix (before number) of second model pdb files
# 3. chains: array containing characters of chains of Ag
#					e.g. ['O', 'R', 'T']
# 4. dist: num - angstrom threshold for determining 'contact'
# 5. normalized: boolean - set to True to account for varying numbers of dock models for a given 
# 												 Ab. Will set absolute parameter to True as well for allContacts()
#
# RETURNS 
# 	outputs a file for each docking model and prints the total number of Ab 
#			residues contacting each Ag res
#		ouptuts one file for Ag contacts accumulated over all docking models
def bulkAllContacts(path, model_pre, chains, dist, normalized=True):

	if not os.path.exists("./contact_num_output"):
		os.makedirs("./contact_num_output", 0777)

	# dictionary to hold Ag res contact totals over all models
	contact_totals = {}

	# Find contacts for each Ag chain for each of the model files
	counter = 0		# to keep track of how many dock files there are to normalize totals
	for i in range(40):
		if i < 10:
			i = "0" + str(i)
		name = model_pre + str(i)
		pdb_name = name+".pdb"
		try:
			file = open(pdb_name, 'r')
		except IOError as e:
			break

		counter += 1

		numbers = allContacts(pdb_name, chains, dist, normalized)
		writeAllContacts(path, name+"_contacts.txt", numbers)

		# increment Ag contacts in contact_totals
		for i in range(len(numbers)):

			res_num = numbers[i][0]
			num_contacts = numbers[i][1]
			if res_num in contact_totals:
				contact_totals[res_num] += num_contacts
			else:
				contact_totals[res_num] = num_contacts

		file.close()

	# now contact_totals holds NON NORMALIZED contact totals for each Ag res

	# write the contact totals to one file
	if model_pre[4:7] == "mut":
		output = open("./contact_num_output/totals_"+model_pre[0:10]+".txt", 'w')
	else:
		output = open("./contact_num_output/totals_"+model_pre[0:6]+".txt", 'w')
	for key in sorted(contact_totals):
		num_contacts = contact_totals[key]

		# normalize if necessary 
		if normalized:
			num_contacts = round(float(num_contacts)/counter, 2)

		output.write("Res #" + key+": "+str(num_contacts)+"\n")










