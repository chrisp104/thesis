import os
import random

# scripts to take output of triple_mutation.py (ranked k-mutation variants that disrupt binding in all
# affected docking models) and run K-MEDOIDS (k-dimensional) clustering on them
#
# FUNCTIONS
# 1. clusterVariants() - cluster mutations using k-medoids algorithm
# 2. checkCoverage() - check to see if medoid mutations cover all docking models
# 3. rmsdFromPDB()
# 4. hausdorff()



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
def clusterVariants(variantsFile, rankedFile, isdb, k, out_path):
	# read in variants and store into an array
	vFile = open(variantsFile,'r')
	vlines = vFile.readlines()
	vFile.close()

	variants = []
	for line in vlines:
		mutations = line.split('|')[0].strip().split(' ')
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
		# print clusters.keys()

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
				hausD = hausdorff(variant, medoid)
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
					key = p+':'+o
					dist = distances[key]
					curCost += dist

				# if the total cost increased revert back, else do nothing

				if float(curCost) < float(prevCost):
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
	non_medoids = list(mutations)
	for key in clusters:
		non_medoids.remove(key)

	# assign each data point to a medoid
	for mutation in non_medoids:
		closest = 0					# the closest medoid
		dist = 100000000		# the shortest distance to a medoid
		for medoid in clusters:
			key = mutation+':'+medoid
			if closest == 0 or distances[key] < dist:
				closest = medoid
				dist = distances[key]
		clusters[closest][0].append(mutation)
		clusters[closest][1] += dist


	# print results to output
	out = open(out_path, 'w')
	for key in clusters:
		out.write("M = "+key+":\n")
		for p in clusters[key][0]:
			out.write(p+" ")
		out.write("\n\n")
	out.close()


	# return medoid set
	print sorted(clusters.keys())
	return clusters



# 2.
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

	for medoid in medoids:
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




# 3.
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



# 4.
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
	print hausD
	return hausD



distances = rmsdFromPDB("/Users/Chris/GitHub/thesis/mutagenesis/merged_isdb.pdb")
hausdorff(("440mGLU", "438mALA", "434mCYS"), ("293mGLY", "297mGLU", "281mASP"), distances)










# ****** TEST CASES ***********

# clusters = kMedCluster("/Users/Chris/GitHub/thesis/mutagenesis/ranked_mutations/D110m2.txt", 
# 		"/Users/Chris/GitHub/thesis/mutagenesis/merged_isdb.pdb", 8, 
# 		"/Users/Chris/GitHub/thesis/mutagenesis/results.txt")

# print checkCoverage(clusters, "/Users/Chris/GitHub/thesis/mutagenesis/ranked_mutations/D110m2.txt")










