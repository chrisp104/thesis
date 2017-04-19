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
#
# RETURNS 
# 	outputs one file containing percent correct contact for each model pair
# 	at the end after the *** it outputs the order of the Ab models for array maker
def percentContactFromNums(output):
	out = open(output, 'w')
	data = {}

	# read in all of the files and store into data first
	for fn in os.listdir("./"):
		if fn[0] != "D": continue

		file = open(fn, 'r')
		info = file.readlines()
		data[fn[:-4]] = info


	# compute pairwise percent similarities
	for m1 in sorted(data):
		for m2 in sorted(data):
			print m1, m2
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
			out.write(m1+" & "+m2+": "+str(percent)+"\n")

	# now write the sorted names of the Ab models
	out.write("\n***\n")
	for name in sorted(data):
		out.write(name+"\n")

	out.close()


os.chdir("/Users/Chris/GitHub/thesis/contacts/all_models/")
percentContactFromNums("/Users/Chris/GitHub/thesis/contacts/all_models/percents.txt")




