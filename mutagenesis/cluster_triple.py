import os
import random
import itertools

# scripts to take output of triple_mutation.py (ranked k-mutation variants that disrupt binding in all
# affected docking models) and run K-MEDOIDS (k-dimensional) clustering on them
#
# FUNCTIONS
# 1. clusterVariants() - cluster mutations using k-medoids algorithm
# 2. findBestVariants() - takes the clusters from function 1 and then selects a set of variants
#			 										covers all docking models and is most disruptive
# 3. checkCoverage() - check to see if medoid mutations cover all docking models
# 4. rmsdFromPDB()
# 5. hausdorff()





# 1.
# clusterVariants()
#
# take a ranked mutation file that covers mutations aggregated for all 30 docking models
# and cluster them into k group based on k-medoids
# 
# ARGUMENTS
# 1. variants: str - path to the variants file for this Ab
# 2. rankedFile: str - path of the ranked mutations file
# 3. isdb: str - path of isdb Ag pdb file
# 4. k: num - number of clusters to use
# 5. out_path: str - output file path 
#
# RETURNS 
# 	clusters - returns set of k mutations representing the cluster medoids
def clusterVariants(variantsFile, isdb, k, out_path):
	# read in variants and store into an array
	vFile = open(variantsFile,'r')
	vlines = vFile.readlines()
	vFile.close()

	variants = []
	for line in vlines:
		mutations = line.strip().split(' ')
		variant = []
		for mutation in mutations:
			variant.append(mutation)
		tuple_variant = tuple(variant)
		variants.append(tuple_variant)

	# create distance matrix (dictionary)
	distances = rmsdFromPDB(isdb)
	# print sorted(mutations)



	# *** K-MEDOIDS ALGORITHM ***

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
	while new_medoid_found:
		print clusters.keys()

		# 2. Assign every entity to its closest medoid
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
		
		
		# 3. For each cluster search if any of the entities of the cluster lower the average 
		# distance total, if it does select the entity that lowers this coefficient the most 
		# as the medoid for this cluster
		new_medoid_found = False
		new_clusters = {}

		for m in clusters:
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

				if float(curCost) < float(prevCost)+0.000005:	 	# the addition is to deal with identical costs
					#print prevCost, curCost
					best_point = o
					prevCost = curCost

			if best_point != m:
				new_medoid_found = True
				new_clusters[best_point] = [[], prevCost]
			else:
				new_clusters[m] = [[], prevCost]

		clusters = new_clusters



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


	# print results to output
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


	# print just the medoids of all Abs to one file
	medoids_out = open("/Users/Chris/GitHub/thesis/mutagenesis/variant_clusters.txt", 'a')
	medoids_out.write(variantsFile[-10:] + '\n')
	for key in clusters:
		for m in key:
			medoids_out.write(m+',')
		medoids_out.write("|")
	medoids_out.write('\n')
	medoids_out.close()


	# return medoid set
	print sorted(clusters.keys())
	return clusters





# 2.
# findBestVariants()
#
# takes the clusters from function 1 and then selects a set of variants covers all 
# docking models and is most disruptive
# 
# ARGUMENTS
# 1. clusters: dict - return value from function 1
# 2. rankedFile: str - path of the ranked mutations file
#
# RETURNS 
#		tuple (variant set, score)
# 	NOT variantSet: array - set of variants with highest disruption score that cover all docking models
def findBestVariants(clusters, rankedFile):

	# read in rankedFile and store mutations into array
	rFile = open(rankedFile, 'r')
	rlines = rFile.readlines()
	rFile.close()

	# create data structure for each mutation point being
	# key = mutation point
	# value = array of tuples: (mutation resi, [affected models], total disruption score)
	mutations = {}
	for line in rlines:
		number = line[:3]
		data = line.split(' ')[1].split('|')

		mutRes = line[4:7]
		affected = data[0].split(',')
		score = float(data[2][:-1])
		entry = (mutRes, affected, score)

		if number in mutations.keys():
			mutations[number].append(entry)
		else:
			mutations[number] = [entry]

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
	# ******************* CHANGE HERE TO MAKE SURE ALL ARE THERE *****************
	x_count = 0
	covered_count = 0
	coveringSets = []
	for x in itertools.product(clusterArray[0], clusterArray[1], clusterArray[2], clusterArray[3]):
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
			# just take the models of the first possible mutation
			models = info[0][1]
			for model in models:
				if not model in covered:
					covered.append(model)

		# add to coveringSets if it covers all docking models
		if len(covered) == 30:
			coveringSets.append(x)
			covered_count += 1

	print covered_count
	print x_count


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

		if len(covered) == 30:
			candidateSets.append((variantSet, total_score))

	candidateSets = sorted(candidateSets, key=lambda x: x[1])
	bestCandidate = candidateSets[-1]
	print bestCandidate
	return bestCandidate








# 3.
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







# 4.
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







# 5.
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








