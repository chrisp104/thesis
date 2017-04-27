import os
import random

# scripts to take output of analyze.py (ranked mutations that disrupt binding in all
# affected docking models) and run K-MEDOIDS (1-dimensional) clustering on them
#
# FUNCTIONS
# 1. kMedCluster() - cluster mutations using k-medoids algorithm
# 2. rmsdFromPDB()



# 1.
# kMedCluster()
#
# take a ranked mutation file that covers mutations aggregated for all 30 docking models
# and cluster them into k group based on k-medoids
# 
# ARGUMENTS
# 1. rankedFile: str - path of the ranked mutations file
# 2. isdb: str - path of isdb Ag pdb file
#
# RETURNS 
# 	clusters - returns set of k mutations representing the cluster medoids
def kMedCluster(rankedFile, isdb, k):

	# read in rankedFile and store mutations into array
	rFile = open(rankedFile, 'r')
	rlines = rFile.readlines()
	rFile.close()

	mutations = []
	for line in rlines:
		mutPosition = line[:3]
		# mutNumber = int(line.split('|')[1])		# number of docking models some sort of mutation at
		# 																			# this point will disrupt
		# for i in range(mutNumber):
		# 	mutations.append(mutPosition)
		mutations.append(mutPosition)

	# create distance matrix (dictionary)
	distances = rmsdFromPDB(isdb)



	# *** K-MEDOIDS ALGORITHM ***

	# 1. Choose k entities to become the medoids
	# select k random medoids from the data
	clusters = {}
	while len(clusters) != k:
		randomMut = random.choice(mutations)
		if not randomMut in clusters:
			clusters[randomMut] = [[], 0]

	# initialize some variables for the loop
	firstTime = True
	new_medoid_found = True
	

	# THE LOOP UNTIL MEDOIDS DO NOT CHANGE (CONVERGENCE)
	while new_medoid_found:
		print clusters.keys()

		# 2. Assign every entity to its closest medoid
		# remove the medoids from the list
		non_medoids = list(mutations)
		for key in clusters:
			non_medoids.remove(key)

		# assign each data point to a medoid
		for mutation in non_medoids:
			closest = 0				# the closest medoid
			dist = 100000000	# the shortest distance to a medoid
			for medoid in clusters:
				key = mutation+':'+medoid
				if closest == 0 or distances[key] < dist:
					closest = medoid
					dist = distances[key]
			clusters[closest][0].append(mutation)
			clusters[closest][1] += dist

		# # IF FIRST ITERATION calculate average distance total (per cluster)
		# if firstTime:
		# 	prevCost = 0
		# 	for medoid in cluster:
		# 		points = cluster[medoid]
		# 		for point in points:
		# 			key = medoid+':'+point
		# 			dist = distances[key]
		# 			prevCost += dist
		# 	prevCost /= float(k)
		# 	firstTime = False
		
		
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
		closest = 0				# the closest medoid
		dist = 100000000	# the shortest distance to a medoid
		for medoid in clusters:
			key = mutation+':'+medoid
			if closest == 0 or distances[key] < dist:
				closest = medoid
				dist = distances[key]
		clusters[closest][0].append(mutation)
		clusters[closest][1] += dist

	# return medoid set
	print clusters.keys()
	return clusters



# 2.
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



kMedCluster("/Users/Chris/GitHub/thesis/mutagenesis/ranked_mutations/D102m0.txt", 
	"/Users/Chris/GitHub/thesis/mutagenesis/merged_isdb.pdb", 9)





