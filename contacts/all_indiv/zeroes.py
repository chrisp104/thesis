# program to enter in 0s for Ag residues if I didn't print those lines
# ARGUMENTS
# 1. bin_dir: str - path to directory with the contact files organized in Ab dir --> model dir
import os
def zeroes(bin_dir):
	# read in all of the files and store into data first
	for ab_dir in os.listdir(bin_dir):
		if ab_dir[0] != 'D': continue
		os.chdir(bin_dir+"/"+ab_dir)
		for model_dir in os.listdir("./"):
			if model_dir[0] != 'm': continue
			os.chdir(bin_dir+"/"+ab_dir+"/"+model_dir)
			for fn in os.listdir("./"):

				# read in file
				file = open(fn, 'r')
				data = file.readlines()
				file.close()

				new_file = open(fn, 'w')
				for i in range(119, 457, 1):
					line_found = False

					# try to find the line first
					for line in data:
						if line[5:8] == str(i):
							new_file.write(line)
							line_found = True
							break

					# if not then just write out 0
					if not line_found:
						new_file.write("Res #"+str(i)+": 0\n")

				new_file.close()

zeroes("/Users/Chris/GitHub/thesis/contacts/all_indiv/")