from mutate import *
import os


# for directory in os.listdir("/Users/Chris/GitHub/thesis/mutagenesis/one_models"):
# 	if directory[0] != 'D': continue
# 	os.chdir("/Users/Chris/GitHub/thesis/mutagenesis/one_models/"+directory)
# 	for model in os.listdir("."):
# 		if model[0] != 'm': continue
# 		print directory
# 		print model
# 		os.chdir("/Users/Chris/GitHub/thesis/mutagenesis/one_models/"+directory+"/"+model)
# 		bulkDirMutations("/Users/Chris/GitHub/thesis/mutagenesis/one_mutations/"+directory+'/'+model+'/')


for directory in os.listdir("/Users/Chris/GitHub/thesis/mutagenesis/sep_models"):
	if directory[0] != 'D': continue
	os.chdir("/Users/Chris/GitHub/thesis/mutagenesis/sep_models/"+directory)
	for model in os.listdir("."):
		if model[0] != 'n': continue
		print directory
		print model
		os.chdir("/Users/Chris/GitHub/thesis/mutagenesis/sep_models/"+directory+"/"+model)
		bulkDirMutations("/Users/Chris/GitHub/thesis/mutagenesis/sep_mutations/"+directory+'/'+model+'/')