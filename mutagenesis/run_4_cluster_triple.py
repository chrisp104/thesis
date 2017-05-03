# clusters each Ab's mutation variants, then finds best set of variants that covers all
# docking models with the highest score
from cluster_triple import *
import os


os.chdir("/Users/Chris/GitHub/thesis/mutagenesis/one_ranked_mutations/")
for model in os.listdir("."):
	if model[0] != 'D': continue
	print model
	clusters = clusterVariants("/Users/Chris/GitHub/thesis/mutagenesis/one_variants/"+model, 
		"/Users/Chris/GitHub/thesis/mutagenesis/merged_isdb.pdb", 4, 
		"/Users/Chris/GitHub/thesis/mutagenesis/one_clusters/"+model[:-4]+"_cluster.txt")
	print checkCoverage(clusters.keys(), "/Users/Chris/GitHub/thesis/mutagenesis/one_ranked_mutations/"+model)