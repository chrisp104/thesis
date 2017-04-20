import os

# script to take contact_total files for a given docking model and compare those

# 1.
# percentContactFromNums()
#
# take two outputs of residue contacts (second return value of res_contacts first function)
# to be compared for % contacts similarity
#
# PERCENT SAME CONTACT is calculated as the following numerator / denominator
# numerator = sum - smaller number of Ab residues in contact from either dock
# denominator = sum - larger number
#
# ARGUMENTS
# 1. output: str - path to write outpu file to
# 2. order_file: str - file that contains the order of the Abs we should sort the matrix in
# 3. ignore_models: bool - will aggregate data across models for a given Ab
#
# RETURNS 
# 	outputs one file containing percent correct contact for each model pair
# 	at the end after the *** it outputs the order of the Ab models for array maker
def percentContactFromNums(output, order_file="NONE", ignore_models=False):
	out = open(output, 'w')
	data = []

	# read in all of the files and store into data first
	for fn in os.listdir("./"):
		if fn[0] != "D": continue

		file = open(fn, 'r')
		info = file.readlines()
		data.append((fn[:-4], info))
		file.close()


	# sort the data array if specified in params
	sorted_data = []
	if order_file != "NONE":

		order = open(order_file, 'r')
		for name in order.readlines():
			ab = name[0:2] + name[3:5]
			for model in data:
				if model[0][:-2] == ab:
					sorted_data.append(model)

		data = sorted_data
		order.close()

	# aggregate data for each Ab if specified in params
	if ignore_models:

		# find number of models for each Ab
		num_docks = 0
		first_ab = data[0][0][:-2]
		print first_ab
		for model in data:
			if model[0][:-2] == first_ab:
				num_docks += 1
			else:
				break
		print num_docks

		# now iterate through and add 
		temp = []		# to hold the new information - analogous to sorted_data
		current = ''
		for i in range(0, len(data), num_docks):
			ab_name = data[i][0][:-2]
			aggr_nums = {}
			for line in data[i][1]:
				aggr_nums[line[:10]] = 0
			ab_contacts = []

			# loop through that antibody's models
			for j in range(num_docks):
				print i, j
				cur_model = data[i+j]
				for line in cur_model[1]:
					aggr_nums[line[:10]] += int(line[10:])

			for key in sorted(aggr_nums):
				ab_contacts.append(key+str(aggr_nums[key]))
			
			temp.append((ab_name, ab_contacts))

		data = temp
		print data

	# compute pairwise percent similarities
	for m1 in data:
		for m2 in data:
			# print m1[0], m2[0]
			cont1 = m1[1]
			cont2 = m2[1]

			numerator = 0
			denominator = 0
			for i in range(len(cont1)):
				line1 = cont1[i]
				line2 = cont2[i]
				if line1[5:8] != line2[5:8]:
					print "mismatch"

				num1 = int(line1[10:])
				num2 = int(line2[10:])
				numerator += min(num1, num2)
				denominator += max(num1, num2)
			
			if denominator == 0:
				percent = 0
			else:
				percent = float(numerator)/denominator
			out.write(m1[0]+" & "+m2[0]+": "+str(percent)+"\n")

	# now write the sorted names of the Ab models
	out.write("\n***\n")
	for model in data:
		print model[0]
		out.write(model[0]+"\n")

	out.close()


os.chdir("/Users/Chris/GitHub/thesis/contacts/all_models/")
percentContactFromNums("/Users/Chris/GitHub/thesis/contacts/all_models/percents_aggr_ordered.txt", 
	"/Users/Chris/GitHub/thesis/contacts/ab_order.txt", True)




