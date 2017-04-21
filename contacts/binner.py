import os

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
def binByPercent(bin_dir, output):
	out = open(output, 'w')
	data = {}

	# read in all of the files and store into data first
	for ab_dir in os.listdir(bin_dir):
		if ab_dir[0] != 'D': continue
		os.chdir(bin_dir+"/"+ab_dir)
		for model_dir in os.listdir("./"):
			if model_dir[0] != 'm': continue
			os.chdir(model_dir)
			for fn in os.listdir("./"):

				file = open(fn, 'r')
				info = file.readlines()
				data[fn[:-4]]=info
				file.close()


	# compute pairwise percent similarities
	for m1 in data:
		for m2 in data:
			print m1, m2
			if m1[0:4] == m2[0:4]: continue
			cont1 = data[m1]
			cont2 = data[m2]

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


os.chdir("/Users/Chris/GitHub/thesis/contacts/bin_c/")
binByPercent("/Users/Chris/GitHub/thesis/contacts/bin_c/", 
	"/Users/Chris/GitHub/thesis/contacts/bin_c/output.txt")




