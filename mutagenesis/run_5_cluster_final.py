# clusters each Ab's mutation variants, then finds best set of variants that covers all
# docking models with the highest score
from cluster_variants import *
import os

setFound = False
tries = 0
while not setFound:
	tries += 1
	
	# FIND MEDOIDS AND CLUSTERS
	clusters = clusterBestVariants("/Users/Chris/GitHub/thesis/mutagenesis/run_1/medoids_2.txt", 
		"/Users/Chris/GitHub/thesis/mutagenesis/merged_isdb.pdb", 3, 
		"/Users/Chris/GitHub/thesis/mutagenesis/run_1/final_clusters_2.txt")
	
	# FIND BEST CANDIDATE
	# if best_candidate == 0:
	# 	continue
	setFound = True
	if tries == 3:
		break









