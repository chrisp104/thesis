# clusters each Ab's mutation variants, then finds best set of variants that covers all
# docking models with the highest score
from cluster_triple import *
import os

out_medoids = open("/Users/Chris/GitHub/thesis/mutagenesis/medoids_3_2.txt", 'w')
out_best_variants = open("/Users/Chris/GitHub/thesis/mutagenesis/variants_3_2.txt", 'w')

os.chdir("/Users/Chris/GitHub/thesis/mutagenesis/one_ranked_mutations/")
for model in os.listdir("."):
	if model[0] != 'D': continue
	print model

	setFound = False
	tries = 0
	while not setFound:
		tries += 1
		
		# FIND MEDOIDS AND CLUSTERS
		clusters = clusterVariants("/Users/Chris/GitHub/thesis/mutagenesis/one_variants_2/"+model, 
			"/Users/Chris/GitHub/thesis/mutagenesis/merged_isdb.pdb", 3, 
			"/Users/Chris/GitHub/thesis/mutagenesis/one_clusters_k3_2/"+model[:-4]+"_cluster.txt")
		
		# FIND BEST CANDIDATE
		best_candidate, max_covered = findBestVariants(clusters, "/Users/Chris/GitHub/thesis/mutagenesis/one_ranked_mutations_2/"+model)
		if best_candidate == 0:
			continue
		else:
			setFound = True
		if tries == 3:
			break


		# write out the medoids to file
		out_medoids.write(model+",0, | ")
		for medoid in clusters:
			for m in medoid:
				out_medoids.write(m+' ')
			out_medoids.write('| ')
		out_medoids.write(',')
		out_medoids.write("Score: "+str(0)+"\n")
		out_medoids.write('\n')


		#candidates = [best_candidate, second, third]
		candidates = [best_candidate]
		out_best_variants.write(model+","+str(max_covered)+', | ')
		print candidates
		for candidate in candidates:
			# write out_best_variants information to file
			variants = candidate[0]
			score = candidate[1]
			for x in variants:
				for m in x:
					out_best_variants.write(m+' ')
				out_best_variants.write('| ')
			out_best_variants.write(',')
			out_best_variants.write("Score: "+str(score)+"\n")
		out_best_variants.write('\n')

out_medoids.close()
out_best_variants.close()








