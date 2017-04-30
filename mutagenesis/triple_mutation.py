from cluster import *
# to use rmsdFromPDB for isdb

import itertools

# scripts to take ranked mutation output files and create triple mutation variants in which
# the mutations are within 12 angstroms of each other
#
# FUNCTIONS
# 1. createVariant



# 1.
# createVariant()
#
# see above
# 
# ARGUMENTS
# 1. ranked_file: str - the path of the ranked mutation file from which to create variants
# 2. isdb: str - path of isdb protein pdb file
# 3. k: int - the number of mutations per variant
# 4. out_file: str - path to write variants to
#
# RETURNS 
# 	array of arrays, each array representing a variant containing k mutations in the form 114mGLU
def createVariant(ranked_file, isdb, k, out_file):
	# read in rankedFile and store mutations into array
	rFile = open(ranked_file, 'r')
	rlines = rFile.readlines()
	rFile.close()

	# create distance matrix (dictionary)
	distances = rmsdFromPDB(isdb)

	# create all possible sets of k mutants
	mutation_list = []
	for line in rlines:
		if line == '\n': continue
		mutation = line[:7]
		score = float(line.split('|')[2])
		mutation_list.append((mutation, score))

	# create the variants
	variants = []
	for comb in itertools.combinations(mutation_list, k):
		variants.append(comb)

	# now for each k-mutation variant, check pairwise distances and make sure they
	# are all within 12 angstroms. Or else, remove.
	final_variants = []
	for variant in variants:
		#print variant
		remove = False
		for i in range(len(variant)):
			if remove:
				break
			for j in range(i+1, len(variant)):
				m1 = variant[i][0]
				res1 = m1[:3]
				m2 = variant[j][0]
				res2 = m2[:3]

				# if the mutations are at same resiude, remove
				if res1 == res2:
					remove = True

				# look up distance
				key = res1+":"+res2
				distance = distances[key]
				if float(distance) > 12:
					remove = True
					break

		if not remove:
			final_variants.append(variant)

	# write output
	out = open(out_file, 'w')
	for variant in final_variants:
		for i in range(k):
			out.write(variant[i][0]+' ')

		# calculate total disruption score and also write
		score = 0
		for mutation in variant:
			s = mutation[1]
			score += s
		out.write("| "+str(score))

		out.write('\n')

	out.close()


	#print final_variants
	#print len(final_variants)
	return final_variants

# createVariant("/Users/Chris/GitHub/thesis/mutagenesis/ranked_mutations/D102m0.txt", 
# 	"/Users/Chris/GitHub/thesis/mutagenesis/merged_isdb.pdb", 3, "/Users/Chris/GitHub/thesis/mutagenesis/variants/D102m2.txt")









