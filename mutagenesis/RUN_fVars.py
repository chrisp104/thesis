import os
from rank_mutations import *
from random import randint
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
# the equivalent of running run_2 thru run_5 in directory but in a loop

# NOTE: run_1_mutate.py should be run separate first since this does not change
# NOTE: run_6 to make the heat map can be incorporated in later



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
# 	exclude: array - Ag residues numbers to exclude 
#
# RETURNS 
# 	final clusters
#		array of excluded docking models
def runAll(iteration, exclusions, out_dir, num_affected, cutoff_score, n, k1, k2, exclude):

	log_file = open(out_dir+"log.txt", 'a')
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
		1, 0, exclusions, 'm')


	# ************** CREATE VARIANTS ******************

	for fn in os.listdir(out_dir+str(iteration)+'_'+"ranked_mutations/"):
		if fn[0] != 'D': continue
		print fn
		createVariant(out_dir+str(iteration)+'_'+"ranked_mutations/"+fn, 
		"/Users/Chris/GitHub/thesis/mutagenesis/merged_isdb.pdb", n,
		out_dir+str(iteration)+'_variants_n'+str(n)+'/'+fn, exclude)


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

	final_clusters = clusterBestVariants(out_dir+str(iteration)+"_variants_k"+str(k1)+".txt", 
		"/Users/Chris/GitHub/thesis/mutagenesis/merged_isdb.pdb", k1, 
		out_dir+str(iteration)+"_final_clusters.txt")


	# find the mutations that cover the most models across all Abs
	final_clusters_covered = findBestFinalVariants(final_clusters, out_dir+str(iteration)+'_'+"ranked_mutations/all.txt")
	print final_clusters_covered

	log_file.write("The best variants for this round were\n")
	return_clusters = []
	for final_cluster in final_clusters_covered:
		return_clusters.append(list(final_cluster))
		for m in final_cluster:
			log_file.write(m+"\n")

			# add to exclude
			exclude.append(m)

	log_file.close()
		
	print "exclude: "
	print exclude

	print return_clusters
	return return_clusters, exclude






# ******************** THE ACTUAL RUNNING OF IT ALL **********************


# 1.
# literallyRunAll()
# Runs runAll iteratively with exclusions after each variant design 
#
# ARGUMENTS
#		iteration
#		exclusions
#		out_dir: str - directory all output from this run will go into
# 	num_affected: int - for ranking, minimum number of models that need to be affected with + disruption score
# 	cutoff_score: float - for ranking, minimum
#		n: int - number of mutations per variant
#		k1: int - k for k medoids clustering of docking model variants
#		num_binC: int - number of antibodies to "experimentally" test each iteration (not > 9)
#		num_binP: int - number of antibodies to "experimentally" test each iteration (not > 9)
# 	exclude: array - Ag residues numbers to exclude 
#
# RETURNS 
# 	final clusters
#		array of excluded docking models
def literallyRunAll(num_affected, cutoff_score, n, k1, num_binC, num_binP, out_directory):
	if not os.path.exists(out_directory):
		os.makedirs(out_directory)

	# decide which are the random Abs we will experimentally test
	testC = []
	testP = []
	binC = ["D229","D110","D214","D216","D228","D312","D320","D323","D324","D331","D410","D413"]
	binP = ["D106","D204","D206","D301","D302","D305","D307","D318","D430","D431"]

	while len(testC) < num_binC:
		choose = (randint(0,len(binC)-1))
		if not binC[choose] in testC:
			testC.append(binC[choose])

	while len(testP) < num_binP:
		choose = (randint(0,len(binP)-1))
		if not binP[choose] in testP:
			testP.append(binP[choose])

	# combine the lists
	final_abs = []
	for c in testC:
		final_abs.append(c)
	for p in testP:
		final_abs.append(p)


	# START HERE
	log_file = open("/Users/Chris/GitHub/thesis/mutagenesis/"+out_directory+"/log.txt", 'w')
	log_file.write("Experimentally Tested Abs:\n")
	for z in final_abs:
		log_file.write(z+" ")
	log_file.write('\n\n')

	iteration = 1
	exclusions = []
	excludeMutations = []

	while iteration < 7:
		
		nothing_changed = True

		log_file = open("/Users/Chris/GitHub/thesis/mutagenesis/"+out_directory+"/log.txt", 'a')
		log_file.write("Iteration: "+str(iteration)+"\n")
		log_file.write("num_affected: "+str(num_affected)+"\n")
		log_file.write("cutoff_score: "+str(cutoff_score)+"\n")
		log_file.write("num mutations per variant: "+str(n)+"\n")
		log_file.write("k: "+str(k1)+"\n")
		
		final_resis, excludeMutations = runAll(iteration=iteration, exclusions=exclusions,
			out_dir="/Users/Chris/GitHub/thesis/mutagenesis/"+out_directory+"/", 
			num_affected=num_affected, cutoff_score=cutoff_score, 
			n=n, k1=k1, k2=k1, exclude=excludeMutations)

		cur_exclusions, cur_remaining = makeExclusions(final_resis,
			"/Users/Chris/GitHub/thesis/mutagenesis/bins.txt", 
			"/Users/Chris/GitHub/thesis/mutagenesis/all.txt",
			"/Users/Chris/GitHub/thesis/mutagenesis/confirmed_mutations.txt",
			"/Users/Chris/GitHub/thesis/mutagenesis/merged_isdb.pdb",
			final_abs,
			'n')

		new_exclusions = []
		for new_e in cur_exclusions:
			if not new_e in exclusions:
				new_exclusions.append(new_e)
				nothing_changed = False
				exclusions.append(new_e)

		log_file.write("Excluded after this round: "+str(len(new_exclusions))+"\n\n")

		out = open("/Users/Chris/GitHub/thesis/mutagenesis/"+out_directory+"/"+str(iteration)+"_remaining.txt", 'w')
		for ab in sorted(cur_remaining):
			models = cur_remaining[ab]
			for model in models:
				if not model in exclusions:
					out.write(model+'\n')
			out.write('\n')
		out.write("Eliminated "+str(len(new_exclusions))+" models\n")
		out.close()

		# if after the first iteration and we aren't eliminating any models, 
		# check to see if we have honed in on NEAT2 binding domain and if not then increase resolution maybe?

		# ************** THIS IS WHERE TO MANIPULATE PARAMS TO HONE IN ********************
		# if len(new_exclusions) < 20:
			# *********** K
			# k1 += 1

		# ********** NUM AFFECTED
		if num_affected > 1:
			num_affected -= 1

		iteration += 1

	log_file.close()


for i in range(10):
	os.chdir("/Users/Chris/GitHub/thesis/mutagenesis/")
	literallyRunAll(num_affected=3, cutoff_score=1, n=3, k1=3, num_binC=5, num_binP=8, 
		out_directory="run_n3_k3_"+str(i))

