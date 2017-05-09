from cluster_triple import rmsdFromPDB
import os

# scripts to take the aggregated ranked mutations, the bins of the Abs,
# the final variant medoids / clusters and "perform" (cross-check with confirmed
# mutagenesis data) mutagenesis analysis to see which docking models to exclude
#
# FUNCTIONS
# 1. makeExclusions()
# 2. findExclusions() - loop through contact files and find which to exclude


# 1.
# makeExclusions()
#
# take medoids of mutations across all Abs, mutate them to Alanine,
# find which models are affected, find affected bins, for all bins eliminate
# docking models based on following:
#		if in affected bin: eliminate models whose contacts are not within N angstroms of correct mutation
#		if in not affected bin: eliminate models whose contacts are within N angstroms of correct mutation
# 
# ARGUMENTS
# 1. clustersFile: str - path of final clusters file
# 2. binsFile: str - path of file with known bins
# 3. rankedFile: str - path of aggregated ranked mutations file
# 4. experimentalFile: str - path of file with confirmed mutations (for validation)
# 5. isdb: str - path of isdb Ag pdb file
#
# RETURNS 
#		dictionary: key - Ab model (e.g. D102), value - list of eliminated contact models
def makeExclusions(clustersFile, binsFile, rankedFile, experimentalFile, isdb):
	# read in all files and convert to data structures
	variants = []
	bins = {}
	mutations = {}

	# clusters and medoids
	cFile = open(clustersFile, 'r')
	for line in cFile.readlines():
		if line[0] == 'M':
			variant = line[4:].strip().split(' ')
			variants.append(variant)
	cFile.close()

	# bins
	bFile = open(binsFile, 'r')
	for line in bFile.readlines():
		key = line[0]
		value = line[2:].strip().split(',')
		bins[key] = value
	bFile.close()

	# experimental mutation data
	mFile = open(experimentalFile, 'r')
	for line in mFile.readlines():
		splitLine = line.split(':')
		mutation = splitLine[0].split('/')
		tuple_mutation = tuple(mutation)
		affected = []
		models = splitLine[1].split(',')
		for m in models:
			m_split = m.split('-')
			tuple_model = (m_split[0], int(m_split[1].strip()))
			affected.append(tuple_model)
		mutations[tuple_mutation] = affected


	# 1. determine which medoid mutations are experimentally disruptive
	disrupted = {}		# hold Ab models that are disrupted by each mutation
	mutagens = []			# hold the mutagens that medoids found to be experimentally disruptive
	
	# loop through each mutation in experimental mutagen variant
	for mutagen in mutations:
		for mutation in mutagen:	# when there are multiple mutations
			
			# loop through each mutation in our generated variant
			for variant in variants:
				for medoid in variant:
					
					if medoid == mutation[:3] and not variant in mutagens:
						mutagens.append(variant)


	# for each mutagen that produced valid experimental results:
	# for mutagen in mutagens:
	# 	for m in mutagen:

	# 		# for each mutation within the mutagen
	# 		m_invalid = False	
	# 		for aa in m:
	# 			if aa[-3:] != "ALA":		# SET TO WHAT AA WE WANT TO MUTATE TO
	# 				m_invalid = True
	# 		if m_invalid: continue
	# for each mutagen that produced valid experimental results LET'S JUST ALA MUTATE FOR NOW:
	for variant in mutagens:
		variant = tuple(variant)
		for mutation in variant:
			m = mutation+"mALA"
			m = [m]
			m = tuple(m)
			if m in mutations.keys():
				# add entry to dictionary if valid mutation
				if not variant in disrupted:
					disrupted[variant] = []		# array to hold models

				data = mutations[m]
				for d in data:
					# if disruption score is below 50 (on 100 scale), add
					if int(d[1]) < 50:
						disrupted[variant].append(d[0])

	print disrupted

	
	# 2. add Ab models to disruption dictionary if in same bin as others
	print disrupted
	new_disrupted = {}
	for mutagen in disrupted:
		new_disrupted[mutagen] = []
		affected = disrupted[mutagen]
		for ab in affected:
			new_disrupted[mutagen].append(ab)
			# find it's bin
			for b in bins:
				ab_models = bins[b]
				# if this ab is in this bin, add all other models to disrupted
				if ab in ab_models:
					for ab_model in ab_models:
						if ab_model in new_disrupted[mutagen]: continue
						new_disrupted[mutagen].append(ab_model)
	disrupted = new_disrupted
	
	exclusions = findExclusions(disrupted, isdb)
	return exclusions




# 2.
# findExclusions()
#
# loop through docking models in all Ab directories
#		if in affected bin: eliminate models whose contacts are not within N angstroms of correct mutation
#		if in not affected bin: eliminate models whose contacts are within N angstroms of correct mutation
# 
# ARGUMENTS
# 1. disrupted: dict - the mutation and disrupted Abs from function 1
# 2. isdb: str - path of isdb Ag pdb file
#
# RETURNS 
#		array - eliminated contact models
def findExclusions(disrupted, isdb):
	exclusions = []
	distances = rmsdFromPDB(isdb)

	# loop through all files once per mutation and eliminate
	for variant in disrupted:
		mutations = []
		affected = disrupted[variant]
		
		for mutation in variant:
			mutations.append(mutation[:3])

		for ab in os.listdir("/Users/Chris/GitHub/thesis/mutagenesis/one_models"):
			if ab[0] != 'D': continue

			# set whether this ab should be disrupted for this mutation or not
			ab_disrupted = False
			if ab in affected:
				ab_disrupted = True

			for mod in os.listdir("/Users/Chris/GitHub/thesis/mutagenesis/one_models/"+ab):
				if mod[0] != 'm': continue
				for d in os.listdir("/Users/Chris/GitHub/thesis/mutagenesis/one_models/"+ab+'/'+mod):
					file = open("/Users/Chris/GitHub/thesis/mutagenesis/one_models/"+ab+'/'+mod+'/'+d)
					lines = file.readlines()

					# if disrupted model, then the file SHOULD contain the current mutation residue
					if ab_disrupted:
						# loop through mutations, since at least one of them SHOULD be in here
						variant_affects = False
						for m in mutations:
							for line in lines:
								if line[0] == '\n': break
								# variant does affect this model if it has contacts at a mutated residue
								if line[:3] == m and len(line) > 9:
									variant_affects = True
						if not variant_affects and not d[:9] in exclusions:
							exclusions.append(d[:9])

					# if NOT disrupted model, then the file should NOT contain any of the current mutation residue
					# so eliminate any that have all three as contacts
					if not ab_disrupted:

						# loop through mutations, and NONE of them should be in here
						all_present = True
						for m in mutations:
							for line in lines:
								if line[0] == '\n': break
								# variant does affect this model if it has contacts at all the mutated residues for the variant
								if line[:3] == m and len(line) == 9:
									if d[:4] == "D430":
										all_present = False
						if all_present and not d[:9] in exclusions:
							exclusions.append(d[:9])

					file.close()
					
	for e in exclusions:
		print e
	print "Ran EXCLUSIONS"
	return exclusions








# makeExclusions("/Users/Chris/GitHub/thesis/mutagenesis/final_clusters.txt",
# 	"/Users/Chris/GitHub/thesis/mutagenesis/bins.txt", 
# 	"/Users/Chris/GitHub/thesis/mutagenesis/all.txt",
# 	"/Users/Chris/GitHub/thesis/mutagenesis/confirmed_mutations.txt",
# 	"/Users/Chris/GitHub/thesis/mutagenesis/merged_isdb.pdb")
