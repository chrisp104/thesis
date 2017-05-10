from rank_mutations import *
from make_exclusions import *
import os


# ********** WITHOUT EXCLUSIONS
for directory in os.listdir("/Users/Chris/GitHub/thesis/mutagenesis/one_mutations"):
	if directory[0] != 'D': continue
	os.chdir("/Users/Chris/GitHub/thesis/mutagenesis/one_mutations/"+directory)
	for model in os.listdir("."):
		if model[0] != 'm': continue
		os.chdir("/Users/Chris/GitHub/thesis/mutagenesis/one_models/"+directory+"/"+model)
		rankForAb("/Users/Chris/GitHub/thesis/mutagenesis/one_mutations/"+directory+'/'+model+'/',
			"/Users/Chris/GitHub/thesis/mutagenesis/run_1/ranked_mutations/"+directory+model+".txt", 2, 1)

rankForAll("/Users/Chris/GitHub/thesis/mutagenesis/one_mutations/",
	"/Users/Chris/GitHub/thesis/mutagenesis/run_1/ranked_mutations/all.txt", 2, 1)



# ************* WITH EXCLUSIONS

# exclusions = makeExclusions("/Users/Chris/GitHub/thesis/mutagenesis/final_clusters.txt",
# 	"/Users/Chris/GitHub/thesis/mutagenesis/bins.txt", 
# 	"/Users/Chris/GitHub/thesis/mutagenesis/all.txt",
# 	"/Users/Chris/GitHub/thesis/mutagenesis/confirmed_mutations.txt",
# 	"/Users/Chris/GitHub/thesis/mutagenesis/merged_isdb.pdb")

# for directory in os.listdir("/Users/Chris/GitHub/thesis/mutagenesis/one_mutations"):
# 	if directory[0] != 'D': continue
# 	os.chdir("/Users/Chris/GitHub/thesis/mutagenesis/one_mutations/"+directory)
# 	for model in os.listdir("."):
# 		if model[0] != 'm': continue
# 		os.chdir("/Users/Chris/GitHub/thesis/mutagenesis/one_models/"+directory+"/"+model)
# 		rankForAb("/Users/Chris/GitHub/thesis/mutagenesis/one_mutations/"+directory+'/'+model+'/',
# 			"/Users/Chris/GitHub/thesis/mutagenesis/one_ranked_mutations_2/"+directory+model+".txt", 2, 1, exclusions)

# rankForAll("/Users/Chris/GitHub/thesis/mutagenesis/one_mutations/",
# 			"/Users/Chris/GitHub/thesis/mutagenesis/one_ranked_mutations_2/all.txt", 2, 1, exclusions)