import os
from rank_mutations import *
from make_exclusions import *
from create_variants import *
from cluster_variants import *


# Vary
# - n - number of mutations / variant
# - k - num clusters for all variants
# - k - num clusters for final variants
# - score and # models in ranking
# - whether i final cluster based on medoids or best variants


# run all scripts in loop to narrow down appropriate docking models
# the equivalent of running run_2 through run_5 in this directory but in a loop

# NOTE: run_1_mutate.py should be run separate first since this does not change
# NOTE: run_6 and 7 to make the heat map can be incorporated in later



# 1.
# runAll()
# 
# ARGUMENTS
#		iteration
#		exclusions
#		out_dir: str - directory all output from this run will go into
# 	num_affected: int - for ranking, minimum number of models that need to be affected with + disruption score
# 	cutoff_score: float - for ranking, minimum
#		n: int - number of mutations per variant
#		k1: int - k for k medoids clustering of docking model variants
#		k2: int - k for k medoids clustering of final variants
#
# RETURNS 
# 	
def runAll(iteration, exclusions, out_dir, num_affected, cutoff_score, n, k1, k2):
	if not os.path.exists(out_dir):
		os.makedirs(out_dir)

	# ************* RUN RANKER ************

	for directory in os.listdir("/Users/Chris/GitHub/thesis/mutagenesis/one_mutations"):
		if directory[0] != 'D': continue
		os.chdir("/Users/Chris/GitHub/thesis/mutagenesis/one_mutations/"+directory)
		for model in os.listdir("."):
			if model[0] != 'm': continue
			os.chdir("/Users/Chris/GitHub/thesis/mutagenesis/one_models/"+directory+"/"+model)
			rankForAb("/Users/Chris/GitHub/thesis/mutagenesis/one_mutations/"+directory+'/'+model+'/',
				out_dir+str(iteration)+'_'+"ranked_mutations/"+directory+model+".txt",
				num_affected, cutoff_score, exclusions)

	rankForAll("/Users/Chris/GitHub/thesis/mutagenesis/one_mutations/",
		out_dir+str(iteration)+'_'+"ranked_mutations/all.txt",
		num_affected, cutoff_score, exclusions)


	# ************** CREATE VARIANTS ******************

	for fn in os.listdir(out_dir+str(iteration)+'_'+"ranked_mutations/"):
		if fn[0] != 'D': continue
		print fn
		createVariant(out_dir+str(iteration)+'_'+"ranked_mutations/"+fn, 
		"/Users/Chris/GitHub/thesis/mutagenesis/merged_isdb.pdb", n,
		out_dir+str(iteration)+'_variants_n'+str(n)+'/'+fn)


	# ************** CLUSTER DOCKING MODEL VARIANTS ************

	out_medoids = open(out_dir+str(iteration)+"_medoids_k"+str(k1)+".txt", 'w')
	out_best_variants = open(out_dir+str(iteration)+"_variants_k"+str(k1)+".txt", 'w')

	os.chdir(out_dir+str(iteration)+'_'+"ranked_mutations/")
	for model in os.listdir("."):
		if model[0] != 'D': continue
		print model

		setFound = False
		tries = 0
		while not setFound:
			tries += 1
			
			# FIND MEDOIDS AND CLUSTERS
			clusters = clusterVariants(out_dir+str(iteration)+'_variants_n'+str(n)+'/'+model, 
				"/Users/Chris/GitHub/thesis/mutagenesis/merged_isdb.pdb", k1, 
				out_dir+str(iteration)+'_clusters_k'+str(k1)+'/'+model[:-4]+"_cluster.txt")
			
			# FIND BEST CANDIDATE
			best_candidate, max_covered = findBestVariants(clusters, out_dir+str(iteration)+'_'+"ranked_mutations/"+model)
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


	# ************** CLUSTER FINAL VARIANTS ************

	clusters = clusterBestVariants(out_dir+str(iteration)+"_medoids_k"+str(k1)+".txt", 
		"/Users/Chris/GitHub/thesis/mutagenesis/merged_isdb.pdb", k2, 
		out_dir+str(iteration)+"_final_clusters.txt")




iteration = 2

while iteration < 4:

	log_file = open(out_dir+str(iteration)+"_log.txt", 'w')

	exclusions = []
	if not iteration == 1:

		exclusions = makeExclusions("/Users/Chris/GitHub/thesis/mutagenesis/run_1/"+str(iteration-1)+"_final_clusters.txt",
			"/Users/Chris/GitHub/thesis/mutagenesis/bins.txt", 
			"/Users/Chris/GitHub/thesis/mutagenesis/all.txt",
			"/Users/Chris/GitHub/thesis/mutagenesis/confirmed_mutations.txt",
			"/Users/Chris/GitHub/thesis/mutagenesis/merged_isdb.pdb")


	out = open("/Users/Chris/GitHub/thesis/mutagenesis/run_1/"+str(iteration-1)+"_exclusions.txt", 'w')
	for e in exclusions:
		out.write(e+'\n')
	out.close()

	runAll(iteration=iteration, exclusions=exclusions,
		out_dir="/Users/Chris/GitHub/thesis/mutagenesis/run_1/", 
		num_affected=2, cutoff_score=1, 
		n=3, k1=3, k2=3)

	iteration += 1







