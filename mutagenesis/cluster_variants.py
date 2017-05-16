import os
import random
import itertools
import time

# scripts to take output of triple_mutation.py (ranked k-mutation variants that disrupt binding in all
# affected docking models) and run K-MEDOIDS (k-dimensional) clustering on them
#
# FUNCTIONS
# 1. clusterVariants() - cluster Ab mutations using k-medoids algorithm
# 2. clusterBestVariants() - cluster returned best variants or medoids representative of all docking
#														 docking models for that Ab
# 3. kMedoids() - the k medoids algorithm
# 4. findBestVariants() - takes the clusters from function 1 and then selects a set of variants
#			 										covers all docking models and is most disruptive
# 5. findBestVariants() - takes the clusters from function 1 and then selects a set of variants
#			 										covers all docking models and is most disruptive
# 6. createHeatMapData()
# 7. checkCoverage() - check to see if medoid mutations cover all docking models
# 8. rmsdFromPDB()
# 9. hausdorff()





# 1.
# clusterVariants()
#
# take a variants file file that covers mutations aggregated for all 30 docking models
# and cluster them into k group based on k-medoids
# 
# ARGUMENTS
# 1. variants: str - path to the variants file for this Ab
# 2. isdb: str - path of isdb Ag pdb file
# 3. k: num - number of clusters to use
# 4. out_path: str - output file path 
#
# RETURNS 
# 	clusters - returns set of k mutations representing the cluster medoids
def clusterVariants(variantsFile, isdb, k, out_path):
	clusters = {}
	out_dir = out_path[:-18]
	print out_dir
	if not os.path.exists(out_dir):
		os.makedirs(out_dir)

	# read in variants and store into an array
	vFile = open(variantsFile,'r')
	vlines = vFile.readlines()
	vFile.close()

	variants = []
	for line in vlines:
		mutations = line.strip().split(' ')[:3]
		variant = []
		for mutation in mutations:
			variant.append(mutation)
		tuple_variant = tuple(variant)
		variants.append(tuple_variant)

	if len(variants) < k:
		array = []
		for variant in variants:
			array.append(variant)
			clusters[variant] = [array, 0]
		return clusters

	# create distance matrix (dictionary)
	distances = rmsdFromPDB(isdb)
	# print sorted(mutations)

	clusters = kMedoids(variants, k, distances)

	# print the clusters for this one Ab to output
	out = open(out_path, 'w')
	for key in clusters:
		out.write("M = ")
		for m in key:
			out.write(m+' ')
		out.write("\n")
		for p in clusters[key][0]:
			for n in p:
				out.write(n+' ')
			out.write('| ')
		out.write("\n\n")
	out.close()


	# return medoid set
	print sorted(clusters.keys())
	return clusters







# 2.
# clusterBestVariants()
#
# takes the clustered medoids or best variants for all Abs and clusters those
# to see where the variants fall / cluster for ALL antibodies
# 
# ARGUMENTS
# 1. variants: str - path to the variants / medoids file for this Ab
# 2. isdb: str - path of isdb Ag pdb file
# 3. k: num - number of clusters to use
# 4. out_path: str - output file path 
#
# RETURNS 
# 	clusters - returns set of k mutations representing the medoids of the variants for all Abs
#		also writes stats to output
def clusterBestVariants(variantsFile, isdb, k, out_path):
	print variantsFile
	# dictionary to count how many times a mutation at a certain position occurs
	mutation_counts = {}
	clusters = {}

	# read in variants and store into an array
	vFile = open(variantsFile,'r')
	vlines = vFile.readlines()
	vFile.close()

	variants = []
	for line in vlines:
		if line == '\n': continue
		variants_line = line.strip().split(',')[2].strip().split('|')[1:-1]
		for triple_mutation in variants_line:
			variant = []
			mutation = triple_mutation.strip().split(' ')
			for m in mutation:
				# first, increment count
				if m in mutation_counts.keys():
					mutation_counts[m] += 1
				else:
					mutation_counts[m] = 1

				# add to variant array
				variant.append(m)
			tuple_variant = tuple(variant)
			variants.append(tuple_variant)

	if len(variants) < k:
		array = []
		for variant in variants:
			array.append(variant)
			clusters[variant] = [array,0]
		return clusters

	# create distance matrix (dictionary)
	distances = rmsdFromPDB(isdb)

	clusters = kMedoids(variants, k, distances)

	# print the clusters for this one Ab to output
	out = open(out_path, 'w')
	for key in clusters:
		out.write("M = ")
		for m in key:
			out.write(m+' ')
		out.write("\n")
		for p in clusters[key][0]:
			for n in p:
				out.write(n+' ')
			out.write('| ')
		out.write("\n\n")

	# write the frequencies
	out.write('FREQUENCIES\n')
	sorted_counts = sorted(mutation_counts, key=mutation_counts.get, reverse=True)
	for key in sorted_counts:
		out.write(key+": "+str(mutation_counts[key])+'\n')
	out.close()



	# return medoid set
	#print sorted(clusters.keys())
	return clusters








# 3.
# kMedoids()
#
# k-medoids algorithm
# 
# ARGUMENTS
# 1. variants: array - containing all variants
# 2. k
# 3. distances: dict - with the pairwise distances between IsdB atom positions 
#
# RETURNS 
# 	clusters - returns set of k mutations representing the cluster medoids
def kMedoids(variants, k, distances):

	# 1. Choose k entities to become the medoids
	# select k random medoids from the data
	clusters = {}
	while len(clusters) != k:
		randomVar = random.choice(variants)
		if not randomVar in clusters:
			clusters[randomVar] = [[], 0]

	# initialize some variables for the loop
	firstTime = True
	new_medoid_found = True
	

	# THE LOOP UNTIL MEDOIDS DO NOT CHANGE (CONVERGENCE)
	iteration = 1
	while new_medoid_found and iteration < 50:
		print clusters.keys()

		# 2. Assign every entity to its closest medoid
		# remove the medoids from the list
		non_medoids = []
		removed = []			# so we don't remove a variant twice just because it's in there mutliple times
		for variant in variants:
			if variant in clusters.keys() and not variant in removed:
				removed.append(variant)
				continue
			non_medoids.append(variant)

		# assign each data point to a medoid using the Hausdorff distance
		for variant in non_medoids:
			closest = 0				# the closest medoid
			dist = 100000000	# the shortest distance to a medoid
			for medoid in clusters:
				hausD = hausdorff(variant, medoid, distances)
				if closest == 0 or hausD < dist:
					closest = medoid
					dist = hausD
			clusters[closest][0].append(variant)
			clusters[closest][1] += dist
		
		
		# 3. For each cluster search if any of the entities of the cluster lower the average 
		# distance total, if it does select the entity that lowers this coefficient the most 
		# as the medoid for this cluster
		new_medoid_found = False
		new_clusters = {}
		for m in clusters:

			# print clusters[m][1]
			prevCost = clusters[m][1]

			curCluster = list(clusters[m][0])		# all the data points in this cluster
			curCluster.append(m)

			# find if any points lower the average total distance
			best_point = m
			for o in curCluster:													# o is the new temp medoid
				curCost = 0
				placeholder = list(curCluster)
				placeholder.remove(o)												# the new cluster with the temp medoid removed

				# *** recompute the cost (sum of dist of points to the new medoid o)
				for p in placeholder:
					dist = hausdorff(p, o, distances)
					curCost += dist

				# if the total cost increased revert back, else do nothing
				if float(curCost) < (float(prevCost)-0.000005):	 	# to deal with identical costs
					# print curCost, prevCost
					best_point = o
					prevCost = curCost

			if best_point != m:
				new_medoid_found = True
				new_clusters[best_point] = [[], 0]
			else:
				new_clusters[m] = [[], 0]

		clusters = new_clusters
		iteration += 1



	# ***** CLUSTERING DONE *********


	# Assign every entity to its closest medoid
	# remove the medoids from the list
	non_medoids = []
	for variant in variants:
		if variant in clusters.keys():
			continue
		non_medoids.append(variant)

	# assign each data point to a medoid using the Hausdorff distance
	for variant in non_medoids:
		closest = 0				# the closest medoid
		dist = 100000000	# the shortest distance to a medoid
		for medoid in clusters:
			hausD = hausdorff(variant, medoid, distances)
			if closest == 0 or hausD < dist:
				closest = medoid
				dist = hausD
		clusters[closest][0].append(variant)
		clusters[closest][1] += dist

	return clusters






# 4.
# findBestVariants()
#
# takes the clusters from function 1 and then selects a set of variants that covers the maximum 
# number of docking models and is most disruptive
# 
# ARGUMENTS
# 1. clusters: dict - return value from function 1
# 2. rankedFile: str - path of the ranked mutations file
#
# RETURNS 
#		tuple (variant set, score)
# 	NOT variantSet: array - set of variants with highest disruption score that cover all docking models
def findBestVariants(clusters, rankedFile):
	# affected number of models
	affectedNumber = 0

	# read in rankedFile and store mutations into array
	rFile = open(rankedFile, 'r')
	rlines = rFile.readlines()
	rFile.close()

	# create data structure for each mutation point being
	# key = mutation point
	# value = array of tuples: (mutation resi, [affected models], total disruption score)
	covered = []
	mutations = {}
	for line in rlines:
		number = line[:3]
		data = line.split(' ')[1].split('|')

		mutRes = line[4:7]
		affected = data[0].split(',')
		score = float(data[2][:-1])
		entry = (mutRes, affected, score)

		# find how many models are affected
		for model in affected:
			if not model in covered:
				covered.append(model)
		affectedNumber = len(covered)

		if number in mutations.keys():
			mutations[number].append(entry)
		else:
			mutations[number] = [entry]

	print "affected: "+str(affectedNumber)
	for key in mutations:
		mutations[key] = sorted(mutations[key], key=lambda x: x[2], reverse=True)

	# add the medoids to the cluster arrays
	clusterArray = []						# create array to hold the clusters
	for medoid in clusters:
		points = clusters[medoid][0]
		points.append(medoid)
		clusters[medoid] = points
		clusterArray.append(points)

	
	# find all variant sets, one from each cluster, that CAN cover all docking models
	x_count = 0
	covered_count = 0
	max_covered = 0
	coveringSets = []
	for x in itertools.product(*clusterArray):
		x_count += 1

		# check if the mutations cover all docking models
		# make array of affected points
		covered = []
		points = []
		for variant in x:
			for p in variant:
				points.append(p)

		# for each mutation point in the set
		for p in points:
			info = mutations[p]
			# for info in mutationList:
			models = info[0][1]
			for model in models:
				if not model in covered:
					covered.append(model)

		# add to coveringSets if it covers all docking models
		if len(covered) > max_covered:
			max_covered = len(covered)
			coveringSets = [x]
		if len(covered) == max_covered:
			coveringSets.append(x)
			covered_count += 1

	print "Number of candidate sets: " + str(covered_count)
	print "Number of models covered: " + str(max_covered)
	#print x_count


	# *** FIND BEST SCORING MODEL ****
	# loop through sets that cover all Abs
	candidateSets = []		# to hold array of tuples (set, score)
	for variantSet in coveringSets:
		# make array of affected points
		covered = []
		points = []
		for variant in variantSet:
			for p in variant:
				points.append(p)

		total_score = 0
		for p in points:
			info = mutations[p]
			best = info[0]

			residue = best[0]
			models = best[1]
			score = best[2]

			total_score += score

			# check
			for model in models:
				if not model in covered:
					covered.append(model)

		if len(covered) == max_covered:
			candidateSets.append((variantSet, total_score))

	candidateSets = sorted(candidateSets, key=lambda x: x[1])
	if len(candidateSets) == 0:
		return 0
	bestCandidate = candidateSets[-1]
	print bestCandidate
	return bestCandidate, max_covered







# 5.
# findBestFinalVariants()
#
# takes the clusters from function 2 and then selects a set of variants that covers the maximum 
# number of docking models and is most disruptive
# 
# ARGUMENTS
# 1. clusters: dict - return value from function 2
# 2. rankedAllFile: str - path of the ranked mutations file for ALL docking models - "all.txt" in ranked mutations dir
#
# RETURNS 
#		tuple (variant set, score)
# 	NOT variantSet: array - set of variants with highest disruption score that cover all docking models
def findBestFinalVariants(clusters, rankedAllFile):
	# affected number of models
	affectedNumber = 0

	# read in rankedFile and store mutations into array
	rFile = open(rankedAllFile, 'r')
	rlines = rFile.readlines()
	rFile.close()

	# create data structure for each mutation point being
	# key = mutation point
	# value = array of tuples: (mutation resi, [affected models], total disruption score)
	covered = []
	mutations = {}
	for line in rlines:
		number = line[:3]
		data = line.split(' ')[1].split('|')

		mutRes = line[4:7]
		affectedNumber = data[0]
		score = float(data[1][:-1])
		entry = (mutRes, affectedNumber, score)

		if number in mutations.keys():
			mutations[number].append(entry)
		else:
			mutations[number] = [entry]

	print "affected: "+str(affectedNumber)
	for key in mutations:
		mutations[key] = sorted(mutations[key], key=lambda x: x[2], reverse=True)

	# make each cluster its own array and store in cluster_array
	cluster_array = []
	for cluster in clusters:
		cur_cluster = []
		cur_cluster.append(cluster)
		variants = clusters[cluster][0]
		for v in variants:
			cur_cluster.append(v)
		cluster_array.append(cur_cluster)

	# now for each cluster, find the variant with the mutations that disrupts the max docking models
	clusters_ranked = []
	for cluster in cluster_array:
		c_ranked = []
		for variant in cluster:
			variant_viable = True 	# to set to False if mutation was not disruptive
			total_affected = 0

			# check first that all mutations are present in mutations
			# if not, then across all Ab docking models, the total disruption score was negative
			for mutation in variant:
				if not mutation in mutations.keys():
					variant_viable = False
					continue

			# otherwise add 
			if variant_viable:
				for mutation in variant:
					number_affected = mutations[mutation][0][1]
					total_affected += int(number_affected)
				c_ranked.append((variant, total_affected))
		c_ranked = sorted(c_ranked, key=lambda x: x[1], reverse=True)
		clusters_ranked.append(c_ranked)
	print clusters_ranked


	final_variants = []
	for final_cluster in clusters_ranked:
		final_variants.append(final_cluster[0][0])
	print final_variants
	return final_variants







# 6.
# createHeatMapData()
#
# take written data output from bulk running function 2 and create data file for heat_mapper.py
# 
# ARGUMENTS
# 1. setFile: str - path of the file containing each model's variant set
# 2. orderFile: str - path of the file containing order of antibodies
# 3. isdb: str - path of isdb Ag pdb file
# 4. outFile: str - path of output file
#
# RETURNS 
# 	outputs a file containing all pairwise average Hausdorff distances for the pairwise variant sets
def createHeatMapData(setFile, orderFile, isdb, outFile):
	out = open(outFile, 'w')

	# create distance matrix (dictionary)
	distances = rmsdFromPDB(isdb)

	initial_data = {}
	sFile = open(setFile, 'r')
	for line in sFile:
		if line == "\n": continue
		ab = line[:4]
		print ab
		variants = line.split(',')[2].strip()[1:-1].split('|')
		new_variants = []
		for variant in variants:
			variant = variant.strip().split(' ')
			new_variants.append(variant)
		initial_data[ab] = new_variants

	# order the data according to the order file
	data = []		# will be 2 element array of Ab name and data
	order = open(orderFile, 'r')
	for name in order.readlines():
		data.append([name[:4], initial_data[name[:4]]])

	order.close()

	# loop through pairwise and compute average hausdorff distance and write out to file
	k = len(data[0][1])
	print k

	# pairwise variant sets
	for m1 in data:
		for m2 in data:
			# pairwise variants
			d1 = m1[1]
			d2 = m2[1]
			total_dist = 0
			used_j = []
			for i in range(len(d1)):
				low = 10000000		# to only sum the lowest hausdorff distance for a given variant
				cur_j = 0
				for j in range(len(d2)):
					v1 = d1[i]
					v2 = d2[j]
					dist = hausdorff(v1, v2, distances)
					if dist < low and not j in used_j:
						low = dist
						cur_j = j
					if m1 == m2:
						print v1
						print v2
						print low
						print cur_j
				used_j.append(cur_j)
				total_dist += low

			average_dist = total_dist / (k)
			print average_dist
			out.write(m1[0]+" & "+m2[0]+": "+str(average_dist)+"\n")

	# now write the sorted names of the Ab models
	out.write("\n***\n")
	for model in data:
		out.write(model[0]+"\n")

	out.close()







# 7.
# checkCoverage()
#
# take the medoids from function 1 and check to see if (they) mutations cover all docking models
# 
# ARGUMENTS
# 1. medoids: array - array containing the k Ag residue medoids
# 2. rankedFile: str - path of the ranked mutations file
#
# RETURNS 
# 	(True, 0) if covers all models, (False, [n1, n2, ...]) if some missing with an array of the 
# 	missing models
def checkCoverage(medoids, rankedFile):
	covered = []

	# read in rankedFile and store mutations into array
	rFile = open(rankedFile, 'r')
	rlines = rFile.readlines()
	rFile.close()

	# put single points into medoids
	single_medoids = []
	for m in medoids:
		for m1 in m:
			single_medoids.append(m1)
	print single_medoids	

	for medoid in single_medoids:
		for line in rlines:
			# if the this line corresponds to a medoid, then parse data
			if medoid != line[:3]: continue

			# if the mutation is one that happens on a medoid, then note which models it covers
			models = line.split(':')[1].strip().split('|')[0].split(',')
			for model in models:
				num = int(model)
				if not (num in covered):
					covered.append(num)

	not_covered = []
	for i in range(30):
		if not (i in covered):
			not_covered.append(i)

	if len(not_covered) == 0:
		return (True, 0)
	else:
		return (False, not_covered)







# 8.
# rmsdFromPDB()
#
# take pdb file lines and create distance matrix (dictionary) for pairwise CA atoms of every residue
# 
# ARGUMENTS
# 1. pdb: str - file path of pdb
#
# RETURNS 
# 	distances - dictionary: key = resi1:resi2, value = Euclidean distance
def rmsdFromPDB(pdb):
	isdbFile = open(pdb, 'r')
	ilines = isdbFile.readlines()
	isdbFile.close()

	# array containing CA atom residue lines
	ca = []
	for line in ilines:
		if line[13:15] == "CA":
			ca.append(line)

	# now loop through array and create distance matrix
	distances = {}
	for i in range(len(ca)):
		for j in range(len(ca)):
			l1 = ca[i]
			l2 = ca[j]

			res1 = l1[23:26]
			res2 = l2[23:26]
			key = res1+':'+res2

			xsq = (float(l1[30:38].strip()) - float(l2[30:38].strip()))**2
			ysq = (float(l1[38:46].strip()) - float(l2[38:46].strip()))**2
			zsq = (float(l1[46:54].strip()) - float(l2[46:54].strip()))**2

			squared = xsq + ysq + zsq
			rmsd = squared**(0.5)

			distances[key] = rmsd

	return distances







# 9.
# hausdorff()
#
# calculate the hausdorff distance between two medoids - basically the maximum distance two points
# in each data set can possibly be from each other
# 
# ARGUMENTS
# 1. v1: tuple - variant 1
# 2. v2: tuple 
# 3. distances: dict - with the pairwise distances between IsdB atom positions
#
# RETURNS 
# 	float - the hausdorff distance between the two
def hausdorff(v1, v2, distances):
	dist_array = []
	for p1 in v1:
		low = -1
		for p2 in v2:
			res1 = p1[:3]
			res2 = p2[:3]
			key = res1+":"+res2
			dist = float(distances[key])
			if low == -1:
				low = dist
			elif dist < low:
				low = dist
		dist_array.append(low)
	hausD = max(dist_array)
	return hausD














# ****** TEST CASES ***********

# clusters = kMedCluster("/Users/Chris/GitHub/thesis/mutagenesis/ranked_mutations/D110m2.txt", 
# 		"/Users/Chris/GitHub/thesis/mutagenesis/merged_isdb.pdb", 8, 
# 		"/Users/Chris/GitHub/thesis/mutagenesis/results.txt")

# print checkCoverage(clusters, "/Users/Chris/GitHub/thesis/mutagenesis/ranked_mutations/D110m2.txt")

# distances = rmsdFromPDB("/Users/Chris/GitHub/thesis/mutagenesis/merged_isdb.pdb")
# hausdorff(("440mGLU", "438mALA", "434mCYS"), ("293mGLY", "297mGLU", "281mASP"), distances)

# clusters = clusterVariants("/Users/Chris/GitHub/thesis/mutagenesis/variants_one/D102m2.txt", 
# 	"/Users/Chris/GitHub/thesis/mutagenesis/merged_isdb.pdb", 4, 
# 	"/Users/Chris/GitHub/thesis/mutagenesis/results.txt")
# findBestVariants(clusters, "/Users/Chris/GitHub/thesis/mutagenesis/ranked_mutations_one/D102m2.txt")

# print checkCoverage(clusters.keys(), "/Users/Chris/GitHub/thesis/mutagenesis/ranked_mutations/D102m0.txt")


# createHeatMapData("/Users/Chris/GitHub/thesis/mutagenesis/bestVariants.txt",
# 	"/Users/Chris/GitHub/thesis/mutagenesis/ab_order.txt",
# 	"/Users/Chris/GitHub/thesis/mutagenesis/merged_isdb.pdb",
# 	"/Users/Chris/GitHub/thesis/mutagenesis/heat_data.txt")






