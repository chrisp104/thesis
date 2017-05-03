import os
import shutil

for antibody in os.listdir("/Users/Chris/GitHub/thesis/antibodies/"):
	if antibody[0] != 'D': continue
	os.chdir("/Users/Chris/GitHub/thesis/antibodies/"+antibody+"/decoys")
	scores = open('scores.txt', 'r')
	lines = scores.readlines()
	scores.close()
	if len(lines) == 0:
		continue
	first = lines[0]
	best = first[11:12]
	model = 'm'+best
	os.makedirs("/Users/Chris/GitHub/thesis/mutagenesis/one_models/"+antibody+"/")
	# os.makedirs("/Users/Chris/GitHub/thesis/mutagenesis/one_models/"+antibody+"/"+model+"/")
	shutil.copytree("/Users/Chris/GitHub/thesis/mutagenesis/two_models/"+antibody+"/"+model+"/", 
		"/Users/Chris/GitHub/thesis/mutagenesis/one_models/"+antibody+"/"+model+"/")