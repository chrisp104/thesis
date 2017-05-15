import os

# script to create sets of mutations that disrupt binding between the Ab and Ag
# FUNCTIONS
# 1. findMutations - will find optimal mutations for a given file
# 2. writeMutations - writes out result from findMutations to file
# 3. bulkDirMutations - executes functions 1 and 2 for all files in a given directory



# 1.
# findMutations()
#
# take contact output file from res_contacts.py and find optimal Ag mutations
# to disrupt binding
#
# ARGUMENTS
# 1. matrix: dictionary - mutation disruption dictionary with
#			key = "RES RES", value = int(score)
# 2. file: str - file name of resi contact information
#
# RETURNS 
# 	DICTIONARY representing an Ag residue and disruptive mutations
#		with key = Ag residue number, 
#				 value = ARRAY of TUPLES of (mutation resi, score)
def findMutations(matrix, file):
	AA = ['ARG', 'LYS', 'ASN', 'ASP', 'GLN', 'GLU', 'HIS', 'PRO', 'TYR', 'TRP',
	'SER', 'THR', 'GLY', 'ALA', 'MET', 'CYS', 'PHE', 'LEU', 'VAL', 'ILE']

	THRESHOLD = 0		# net disruption score for Ag residue above which counts as good

	f = open(file, 'r')
	data = f.readlines()
	f.close()

	mut_data = {}		# return data

	start_found = False
	non_lines = ['A', 'O', 'R', 'T', 'I', 'E']
	for line in data:

		# # only loop through resi info
		# if line[:-1] == "START":
		# 	start_found = True
		# 	continue
		# if line[:-1] == "END":
		# 	break
		# if not start_found:
		# 	continue
		if line[0] in non_lines:
			continue

		# loop through the Ag resi lines
		split_line = line.split(':')
		ag = split_line[0]
		ab = split_line[1]

		ag_split = ag.split(',')
		ag_num = ag_split[0]
		ag_res = ag_split[1]

		ab_resis = ab.split('-')

		# if no contacts then just set to null, effectively
		if len(line) == 9: 
			mut_data[ag_num] = 0
			continue

		# loop through all 19 possible mutations for Ag res
		aa_scores = []
		for aa in AA:
			if ag_res == aa: continue

			net_score = 0		# to keep track of mutation score over all contacting Ab resis

			# loop through ab_resis 
			for ab_res_info in ab_resis:
				if ab_res_info == '\n': continue
				info_split = ab_res_info.split(',')
				ab_num = info_split[0]
				ab_res = info_split[1]
				ab_chain = info_split[2]

				current_pair = ag_res+':'+ab_res	# actual Ag:Ab contact
				mutated_pair = aa+':'+ab_res

				cur_score = float(matrix[current_pair])
				mut_score = float(matrix[mutated_pair])

				difference = mut_score - cur_score

				# print current_pair, cur_score
				# print mutated_pair, mut_score
				# print difference
				net_score += difference

			# now write in data for this mutation
			# if net_score > THRESHOLD:		# only append if score is net disruptive
			aa_scores.append((aa, net_score))

		# now write in data for all mutations for this Ag res
		aa_scores = sorted(aa_scores, key=lambda tup: tup[1], reverse=True)
		mut_data[ag_num] = aa_scores

	return mut_data






# 2.
# writeMutations()
#
# write out return value from findMutations into a file
#
# ARGUMENTS
# 1. data: return value from findMutations (1.)
# 2. out_dir: str - directory path to write output file to
# 3. out_name: str - file name of output
#
# RETURNS
#		an output file with the best mutations for that Ag resi, starting with 
#		most disruptive to least
def writeMutations(data, out_dir, out_name):
	if not os.path.exists(out_dir):
		os.makedirs(out_dir)
	output = open(out_dir+out_name, 'w')

	for key in sorted(data):
		output.write(key+": ")

		# if no contacts
		if data[key] == 0:
			output.write('\n')
			continue

		first = True
		for mutation in data[key]:
			if first:
				output.write(mutation[0]+'|'+str(mutation[1]))
				first = False
				continue
			output.write(','+mutation[0]+'|'+str(mutation[1]))
		output.write('\n')

	output.close()








# 3.
# bulkDirMutations()
#
# executes functions 1 and 2 for all files in a given directory
#
# ARGUMENTS
# 1. out_dir: str - directory path to write output file to
#
# RETURNS
#		an output file with the best mutations for that Ag resi, starting with 
#		most disruptive to least
def bulkDirMutations(out_dir):

	# Read in resi pair disruption score matrix file
	matrix_file = open("/Users/Chris/GitHub/thesis/mutagenesis/sipper_i.mtx", 'r')
	lines = matrix_file.readlines()
	matrix_file.close()

	matrix = {}

	for line in lines:
		pair = line[:3]+':'+line[4:7]
		score = line[8:-1]
		matrix[pair] = score


	# Find contacts for each Ag chain for each of the model files
	for fn in os.listdir('.'):
		if fn[0:1] == '.' or fn[-12:] != "contacts.txt": continue
		try:
			file = open(fn, 'r')
		except IOError as e:
			break

		print fn

		mut_data = findMutations(matrix, fn)
		print mut_data

		writeMutations(mut_data, out_dir, fn[:-12]+"_mutations.txt")

		file.close()



os.chdir("/Users/Chris/GitHub/thesis/contacts/crystals/")
bulkDirMutations("/Users/Chris/GitHub/thesis/contacts/crystals/")









