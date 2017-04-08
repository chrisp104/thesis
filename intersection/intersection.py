import os

#	To use with the results of all_contacts (bulk_all_contacts especially)
#
#	1.
# intersection()
# 
# takes in an array of file names, where each file has the number of Ab contacts
# for each Ag residue, and finds percentage intersection between them
#
# ARGUMENTS
# 1. files = array: array of file names with the results
#	2. threshold = how many contacts / percent determines a "contact" across 
# 3. out = output file name where results will be written
#
# RETURNS 
# 	1. array of tuples, tuple = (str - model comparison pair, int - percentage intersection)
#
#	OUTPUTS
#		a file with the docking model pairs and their overlapping contact percentages
def intersection(files, threshold, out):

	# stores all files
	data = {}

	# return value
	results = []

	# read in all files and store
	for file in files:
		f = open(file, 'r')
		lines = f.readlines()
		data[file] = lines
		f.close()

	# nested loop to go through all pairs of files
	for i in range(len(files)):
		# for j in range(i+1, len(files)): 	# if we don't want redundancy
		for j in range(len(files)):
			print i, j
			file_one = files[i]
			file_two = files[j]

			# retrieve the file data
			one = data[file_one]
			two = data[file_two]

			# loop through Ag residues and add up total intersection
			numerator = 0
			denominator = 0
			for k in range(len(one)):
				line_one = one[k]
				line_two = two[k]

				# **** RESIDUE NUMBERS MUST BE IN EXACT SAME ORDER ACROSS ALL FILES
				res_one = line_one[5:8]
				res_two = line_two[5:8]
				total_one = float(line_one[10:].strip())
				total_two = float(line_two[10:].strip())

				# numerator incremented if both are above threshold
				# denominator incremented if one is above threshold
				if total_one >= threshold or total_two >= threshold:
					denominator += 1
				if total_one >= threshold and total_two >= threshold:
					numerator += 1

			# write to return array
			percent_intersect = round(numerator/float(denominator), 3)
			result = (file_one[7:13] + " & " + file_two[7:13], percent_intersect)
			results.append(result)

	output = open(out, 'w')
	for result in results:
		output.write(result[0]+": "+str(result[1])+"\n")

	output.close()

	print results
	return results

# all models
intersection(["totals_D206m2.txt", "totals_D206m9.txt", "totals_D206mR.txt",
	"totals_D410m1.txt", "totals_D410m5.txt", "totals_D410mR.txt",
	"totals_D430m7.txt", "totals_D430m8.txt", "totals_D430mR.txt"], 0.05, "output_all.txt", )

# just crystal structure models
intersection(["totals_D206mR.txt", "totals_D410mR.txt", "totals_D430mR.txt"], 0.05, "output_crystals.txt", )












