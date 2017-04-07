import os

#	To use with the results of all_contacts (bulk_all_contacts especially)
#
#	1. intersection()
# 
# takes in an array of file names, where each file has the number of Ab contacts
# for each Ag residue, and finds percentage intersection between them
#
# ARGUMENTS
# files = array: array of file names with the results
# out = output file name where results will be written
#
# RETURNS 
# 	1. array of tuples, tuple = (model comparison pair, percentage intersection)

def intersection(files, out):

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
		for j in range(i+1, len(files)):
			file_one = files[i]
			file_two = files[j]

			# retrieve the file data
			one = data[file_one]
			two = data[file_two]

			# loop through Ag residues and add up total intersection
			intersection = 0
			for k in range(len(one)):
				line_one = one[k]
				line_two = two[k]

				# **** RESIDUE NUMBERS MUST BE IN EXACT SAME ORDER ACROSS ALL FILES
				res_one = line_one[5:8]
				res_two = line_two[5:8]
				total_one = int(line_one[10:].strip())
				total_two = int(line_two[10:].strip())

				common = abs(total_one - total_two)
				intersection += common

			# write to return array
			result = (file_one + " & " + file_two, intersection)
			results.append(result)

		output = open(out, 'w')
		for result in results:
			output.write(result[0][7:13]+": "+str(result[1])+"\n")

		output.close()

		print results
		return results

intersection(["totals_D206m2.txt", "totals_D206m9.txt", "totals_D206mR.txt",
	"totals_D410m1.txt", "totals_D410m5.txt", "totals_D410mR.txt",
	"totals_D430m7.txt", "totals_D430m8.txt", "totals_D430mR.txt"], "output.txt")














