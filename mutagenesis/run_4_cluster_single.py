from cluster import *
import os


os.chdir("/Users/Chris/GitHub/thesis/mutagenesis/ranked_mutations/")
for model in os.listdir("."):
	if model[0] != 'D': continue
	print model
	clusters = kMedCluster("/Users/Chris/GitHub/thesis/mutagenesis/ranked_mutations/"+model, 
		"/Users/Chris/GitHub/thesis/mutagenesis/merged_isdb.pdb", 20, 
		"/Users/Chris/GitHub/thesis/mutagenesis/clusters/"+model[:-4]+"_cluster.txt")
	print checkCoverage(clusters, "/Users/Chris/GitHub/thesis/mutagenesis/ranked_mutations/"+model)