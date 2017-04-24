from mutate import *
import os


# for directory in os.listdir("/home/anthill/cpark/thesis/mutagenesis"):
	# if directory[0] != 'D': continue
	# os.chdir("/home/anthill/cpark/thesis/antibodies/"+directory)
for model in os.listdir("."):
	if model[0] != 'm': continue
	os.chdir("/Users/Chris/GitHub/thesis/mutagenesis/"+model)
	bulkDirMutations("/Users/Chris/GitHub/thesis/mutagenesis/results/")