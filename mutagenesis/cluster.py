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
		if not mutPosition in mutations:
			mutations.append(mutPosition)


	# create distance matrix (dictionary)
	distances = rmsdFromPDB(isdb)

	# *** K-MEDOIDS ALGORITHM ***

	# select k random medoids from the data
	medoids = []
	while len(medoids) != k:
		randomMut = random.choice(mutations)
		if not randomMut in medoids:
			medoids.append(randomMut)

	prevCost = 1000000000
	curCost = 100000000

	# while cost of configuration decreases
	while curCost < prevCost:
		print medoids
		print curCost
		prevCost = curCost

		# for each medoid m
		for i in range(len(medoids)):
			print i

			# for each non-medoid data point o
			for o in mutations:
				curCost = 0

				# swap m and o, using placeholder in case we need to revert
				placeholder = medoids[i]
				medoids[i] = o

				# *** recompute the cost (sum of dist of points to their medoid)
				# for each non-medoid, compare with each medoid, and add lowest score to curCost
				for mutation in mutations:
					low_dist = 0
					for medoid in medoids:
						key = mutation+':'+medoid
						if low_dist == 0:
							low_dist = distances[key]
						elif distances[key] < low_dist:
							low_dist = distances[key]

					curCost += low_dist

				# if the total cost increased revert back, else do nothing
				if curCost > prevCost:
					medoids[i] = placeholder

	# return medoid set
	return medoids



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
	"/Users/Chris/GitHub/thesis/mutagenesis/merged_isdb.pdb", 1)





