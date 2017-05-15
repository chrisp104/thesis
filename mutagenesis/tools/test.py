from cluster_variants import *


final_clusters = clusterBestVariants("/Users/Chris/GitHub/thesis/mutagenesis/variants.txt", 
	"/Users/Chris/GitHub/thesis/mutagenesis/merged_isdb.pdb", 3, 
	"/Users/Chris/GitHub/thesis/mutagenesis/test_out.txt")

findBestFinalVariants(final_clusters, "/Users/Chris/GitHub/thesis/mutagenesis/all.txt")