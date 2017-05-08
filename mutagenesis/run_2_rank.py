from rank_mutations import *
import os


for directory in os.listdir("/Users/Chris/GitHub/thesis/mutagenesis/one_mutations"):
	if directory[0] != 'D': continue
	os.chdir("/Users/Chris/GitHub/thesis/mutagenesis/one_mutations/"+directory)
	for model in os.listdir("."):
		if model[0] != 'm': continue
		print directory
		print model
		os.chdir("/Users/Chris/GitHub/thesis/mutagenesis/one_models/"+directory+"/"+model)
		rankForAb("/Users/Chris/GitHub/thesis/mutagenesis/one_mutations/"+directory+'/'+model+'/',
			"/Users/Chris/GitHub/thesis/mutagenesis/one_ranked_mutations/"+directory+model+".txt", 2, 1)

rankForAll("/Users/Chris/GitHub/thesis/mutagenesis/one_mutations/",
			"/Users/Chris/GitHub/thesis/mutagenesis/one_ranked_mutations/all.txt", 2, 1)