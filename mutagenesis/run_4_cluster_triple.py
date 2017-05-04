# clusters each Ab's mutation variants, then finds best set of variants that covers all
# docking models with the highest score
from cluster_triple import *
import os

out = open("/Users/Chris/GitHub/thesis/mutagenesis/bestVariants.txt", 'w')

os.chdir("/Users/Chris/GitHub/thesis/mutagenesis/one_ranked_mutations/")
for model in os.listdir("."):
	if model[0] != 'D': continue
	print model

	setFound = False
	tries = 0
	while not setFound:
		tries += 1
		clusters = clusterVariants("/Users/Chris/GitHub/thesis/mutagenesis/one_variants/"+model, 
			"/Users/Chris/GitHub/thesis/mutagenesis/merged_isdb.pdb", 4, 
			"/Users/Chris/GitHub/thesis/mutagenesis/one_clusters/"+model[:-4]+"_cluster.txt")
		best_candidate, max_covered = findBestVariants(clusters, "/Users/Chris/GitHub/thesis/mutagenesis/one_ranked_mutations/"+model)
		print best_candidate
		if best_candidate != 0:
			setFound = True
		if tries == 3:
			break

		#candidates = [best_candidate, second, third]
		candidates = [best_candidate]
		out.write(model+","+str(max_covered)+', | ')
		for candidate in candidates:
			# write out information to file
			variants = candidate[0]
			score = candidate[1]
			for x in variants:
				for m in x:
					out.write(m+' ')
				out.write('| ')
			out.write(',')
			out.write("Score: "+str(score)+"\n")
		out.write('\n')

out.close()