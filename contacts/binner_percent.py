import os
import random

# script to take individual contact files and compare the 30 docking contact files from one
# Ab pairwise with each of the others' 30 docking contact files
#
# *** WILL NEED TO CHOOSE VARIOUS STARTING ANTIBODIES AND REPEAT FOR EACH STARTING POINT ****

# 1.
# binByPercent()
#
# PERCENT SAME CONTACT is calculated as the following numerator / denominator
# numerator = sum - smaller number of Ab residues in contact from either dock
# denominator = sum - larger number
#
# ARGUMENTS
# 1. bin_dir: str - path to directory with the contact files organized in Ab dir --> model dir
# 2. output: str - path to write output file to
#
# RETURNS 
# 	outputs one file containing the sets of docking models for each Ab that creates
#		a set with the others in the set that "cross-block" each other at least wrt contacts
# 	within 8 angstroms
#		
#		returns an array of tuples: (model names that "cross-block", total summed % similarity score)
def binByPercent(bin_dir, output):
	THRESHOLD = 0.35
	out = open(output, 'w')
	data = {}		# to just hold file information
	antibodies = []		# to hold the antibody model names, e.g., D206m1

	# read in all of the files and store into data first
	for ab_dir in os.listdir(bin_dir):
		if ab_dir[0] != 'D': continue
		os.chdir(bin_dir+"/"+ab_dir)
		for model_dir in os.listdir("./"):
			if model_dir[0] != 'm': continue
			antibodies.append(ab_dir+model_dir)
			os.chdir(bin_dir+"/"+ab_dir+"/"+model_dir)
			for fn in os.listdir("./"):

				file = open(fn, 'r')
				info = file.readlines()
				data[fn[:-13]]=info
				file.close()

	# the array of result sets to return
	results = []
	# populate with one model's contact files
	start = random.choice(data.keys())

	antibodies.remove(start[:-3])
	for key in data:
		if key[:-3] == start[:-3]:
			results.append(([key], 0))

	print results
	print antibodies


	# **** ALGORITHM ****

	# compute pairwise percent similarities
	random.shuffle(antibodies)
	for current_ab in antibodies:
		temp_results = []
		for model in data:
			if model[:-3] != current_ab: 
				continue
			print model
			model_data = data[model]	# contact information for current docking model
			for ab_set in results:
				# now we are comparing a given docking model from the current Ab tier
				# to all of the sets of Abs in results

				# to break out of this loop if one of the pairs was under threshold
				set_removed = False

				new_ab_list = []
				new_ab_list.append(model)
				score = ab_set[1]

				# loop through Abs in the set
				for ab in ab_set[0]:

					new_ab_list.append(ab)
					ab_data = data[ab]
					# model_data from above

					numerator = 0
					denominator = 0
					for i in range(len(ab_data)):
						line1 = ab_data[i]
						line2 = model_data[i]
						if line1[5:8] != line2[5:8]:
							print "mismatch"

						num1 = int(line1[10:])
						num2 = int(line2[10:])

						numerator += min(num1, num2)
						denominator += (max(num1, num2))
					
					if denominator == 0:
						percent = 0
					else:
						percent = float(numerator)/denominator

					# if pairwise score is at any point below a threshold then remove that set
					if percent < THRESHOLD:
						set_removed = True
						break

					print percent
					score += percent

				# don't add if one of the scores was too low
				if set_removed: continue

				# now to amend results to accomodate this Ab's models into the sets
				temp_results.append((new_ab_list, score))

		results = temp_results


	# now write the sorted names of the Ab models
	for ab_set in results:
		for ab in ab_set[0]:
			out.write(ab+"-")
		out.write(str(ab_set[1])+"\n")

	out.close()


os.chdir("/Users/Chris/GitHub/thesis/contacts/bin_c/")
binByPercent("/Users/Chris/GitHub/thesis/contacts/bin_c/", 
	"/Users/Chris/GitHub/thesis/contacts/bin_c/output.txt")




