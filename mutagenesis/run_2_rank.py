from analyze import *
import os


for directory in os.listdir("/Users/Chris/GitHub/thesis/mutagenesis/mutations"):
	if directory[0] != 'D': continue
	os.chdir("/Users/Chris/GitHub/thesis/mutagenesis/mutations/"+directory)
	for model in os.listdir("."):
		if model[0] != 'm': continue
		print directory
		print model
		os.chdir("/Users/Chris/GitHub/thesis/mutagenesis/models/"+directory+"/"+model)
		rankMutations("/Users/Chris/GitHub/thesis/mutagenesis/mutations/"+directory+'/'+model+'/',
			"/Users/Chris/GitHub/thesis/mutagenesis/ranked_mutations/"+directory+model+".txt")