# clusters each Ab's mutation variants, then finds best set of variants that covers all
# docking models with the highest score
from cluster_triple import *
import os


os.chdir("/Users/Chris/GitHub/thesis/mutagenesis/ranked_mutations_one/")
for model in os.listdir("."):
	if model[0] != 'D': continue
	print model
	clusters = clusterVariants("/Users/Chris/GitHub/thesis/mutagenesis/variants_one/"+model, 
		"/Users/Chris/GitHub/thesis/mutagenesis/merged_isdb.pdb", 4, 
		"/Users/Chris/GitHub/thesis/mutagenesis/triple_clusters_one/"+model[:-4]+"_cluster.txt")
	print checkCoverage(clusters.keys(), "/Users/Chris/GitHub/thesis/mutagenesis/ranked_mutations_one/"+model)