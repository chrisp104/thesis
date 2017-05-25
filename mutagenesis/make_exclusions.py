from cluster_variants import rmsdFromPDB
import os
import time

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
# 1. variants: array - array of arrays where internal arrays are the variants' mutations
# 2. binsFile: str - path of file with known bins
# 3. rankedFile: str - path of aggregated ranked mutations file
# 4. experimentalFile: str - path of file with confirmed mutations (for validation)
# 5. isdb: str - path of isdb Ag pdb file
# 6. ab_list: array - representative Abs that are being "experimentally tested"
# 7. modelType: str - which directory of models to look at: separate or one
#
#
# RETURNS 
#		list: eliminated contact models
def makeExclusions(variants, binsFile, rankedFile, experimentalFile, isdb, ab_list, modelType='m'):
	# read in all files and convert to data structures
	#variants = []
	bins = {}
	mutations = {}

	# clusters and medoids
	# cFile = open(clustersFile, 'r')
	# for line in cFile.readlines():
	# 	if line[0] == 'M':
	# 		variant = line[4:].strip().split(' ')
	# 		variants.append(variant)
	# cFile.close()

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

			# skip if it's not one we "experimentally tested"
			if not m_split[0] in ab_list: continue

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
	# we will have a list of experimentally disrupted and not disrupted
	# dictionary disrupted: key = mutation, value = [disrupted], [non-disrupted]
	for variant in mutagens:
		variant = tuple(variant)
		for mutation in variant:
			print mutation
			m = mutation+"mALA"
			m = [m]
			m = tuple(m)
			if m in mutations.keys():
				# add entry to dictionary if valid mutation
				if not variant in disrupted:
					disrupted[variant] = [[], []]		# array to hold models

				data = mutations[m]
				for d in data:
					# if disruption score is below 50 (on 100 scale), add
					if int(d[1]) < 50 and not d[0] in disrupted[variant][0]:
						disrupted[variant][0].append(d[0])
					elif not d[0] in disrupted[variant][1]:
						disrupted[variant][1].append(d[0])

	
	# 2. add Ab models to disrupted domains dictionary if in same bin as others to
	# eliminate first based on if it's not in the right domain
	new_disrupted = {}			# holds the actual disruption data and extrapolated disruptions
	domain_disrupted = {}		# holds which ab bins were disrupted as a whole

	for mutagen in disrupted:
		print "mutagen"
		print mutagen
		new_disrupted[mutagen] = []
		domain_disrupted[mutagen] = []
		
		# determine if this mutation spot was more disruptive to more tested models 
		disruptive = False
		disrupted_models = disrupted[mutagen][0]
		non_disr_models = disrupted[mutagen][1]

		if len(disrupted_models) >= len(non_disr_models):
			disruptive = True

		affected = disrupted[mutagen][0]

		for ab in affected:

			new_disrupted[mutagen].append(ab)
			domain_disrupted[mutagen].append(ab)
			
			# find it's bin
			for b in bins:
				ab_models = bins[b]
				print ab_models
				# if this ab is in this bin, add all other models to disrupted that we 
				# don't have experimental data for - this is the extrapolation step  
				if ab in ab_models:
					for ab_model in ab_models:
						if not ab_model in domain_disrupted[mutagen]: 
							domain_disrupted[mutagen].append(ab_model)

						# first add the tested disruptives to the dictionary no matter what
						if not ab_model in new_disrupted[mutagen] and ab_model in disrupted_models:
							new_disrupted[mutagen].append(ab_model)

						# if we have deemed the bin to more likely be disruptive, then
						# add the rest of the non tested models into the disrupted list
						# as long as it wasn't tested to be non disruptive.
						if not ab_model in new_disrupted[mutagen] and disruptive and not ab_model in non_disr_models:
							new_disrupted[mutagen].append(ab_model)


	disrupted = new_disrupted
	print "Disrupted"
	print disrupted
	print "Domains disrupted"
	print domain_disrupted

	exclusions, remaining = findExclusions(disrupted, isdb, modelType, domain_disrupted)
	return exclusions, remaining





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
# 3. modelType: str - which directory of models to look at: separate or one
# 4. domain_disrupted: dict - the domains disrupted using tested positives regardless of what
#				the majority disrutive/non-disruptive consensus was so we can eliminate based off of domain first
#
# RETURNS 
#		array - eliminated contact models
#		remaining - dict: key = Ab, value = list of remaining models
def findExclusions(disrupted, isdb, modelType, domain_disrupted):
	if modelType == 'n':
		mutation_directory = 'sep_mutations'
	else:
		mutation_directory = 'one_mutations'


	exclusions = []
	remaining = {}
	distances = rmsdFromPDB(isdb)

	# loop through all files once per variant and eliminate
	for variant in disrupted:
		mutations = []
		affected = disrupted[variant]
		print domain_disrupted
		domain_affected = domain_disrupted[variant]
		
		for mutation in variant:
			mutations.append(mutation[:3])

		for ab in os.listdir("/Users/Chris/GitHub/thesis/mutagenesis/"+mutation_directory):
			if ab[0] != 'D': continue
			remaining[ab] = []

			# set whether this ab should be disrupted for this mutation or not
			ab_disrupted = False
			cur_domain_disrupted = False
			if ab in affected:
				ab_disrupted = True

			if ab in domain_affected:
				cur_domain_disrupted = True

			for mod in os.listdir("/Users/Chris/GitHub/thesis/mutagenesis/"+mutation_directory+'/'+ab):
				if mod[0] != modelType: continue
				for d in os.listdir("/Users/Chris/GitHub/thesis/mutagenesis/"+mutation_directory+'/'+ab+'/'+mod):
					file = open("/Users/Chris/GitHub/thesis/mutagenesis/"+mutation_directory+'/'+ab+'/'+mod+'/'+d)
					lines = file.readlines()
					# store the mutation contact information in dictionary
					incorrect_domain = True		# to use if this model SHOULD have the mutation spots possible in variant's positions
					model = {}
					for line in lines:
						model[line[:3]] = []

						# set incorrect_domain to true if this has residue spots in variant
						if cur_domain_disrupted and line[:3] in mutations:
							incorrect_domain = False

						if len(line) <= 7: continue
						model_muts = line.strip().split(':')[1].strip().split(',')
						for mut in model_muts:
							split_mut = mut.split('|')
							model[line[:3]].append(split_mut)

					if incorrect_domain and cur_domain_disrupted:
						exclusions.append(d[:9])
						continue

					# if disrupted model, then the file SHOULD contain the current mutation residue
					# excluding the ones that contain none of the residues
					if ab_disrupted:
						# loop through mutations, since at least one of them SHOULD be in here
						variant_affects = False
						for m in mutations:
							if not m in model.keys(): continue
							# variant does affect this model if it has contacts at a mutated residue
							if len(model[m]) != 0:
								variant_affects = True
						if not variant_affects and not d[:9] in exclusions:
							exclusions.append(d[:9])
						else: 
							remaining[ab].append(d[:9])

					# if NOT disrupted model, then the file should NOT contain any of the current mutation residue
					# for now, eliminate any that have A mutation of the variant that is significantly disruptive, >0.5
					if not ab_disrupted:
						# loop through mutations, and NONE of them should be disruptive.
						variant_affects = False
						for m in mutations:
							if not m in model.keys(): continue
							# variant does affect this model if it has contacts at a mutated residue
							cur_model_mutations = model[m]
							for cur_mut in cur_model_mutations:
								if cur_mut[0] == "ALA" and float(cur_mut[1]) > 0.25:	# *************************** CHANGE HERE
									variant_affects = True
						if variant_affects and not d[:9] in exclusions:
							exclusions.append(d[:9])
						else: 
							remaining[ab].append(d[:9])

					file.close()
					
	# for e in exclusions:
	# 	print e
	print len(exclusions)
	print "Ran EXCLUSIONS"
	return exclusions, remaining








# makeExclusions([['339', '407'], ['165', '166']],
# 	"/Users/Chris/GitHub/thesis/mutagenesis/bins.txt", 
# 	"/Users/Chris/GitHub/thesis/mutagenesis/all.txt",
# 	"/Users/Chris/GitHub/thesis/mutagenesis/confirmed_mutations.txt",
# 	"/Users/Chris/GitHub/thesis/mutagenesis/merged_isdb.pdb",
# 	["D110", "D331", "D410", "D229", "D324", "D214", "D431", "D302", "D106", "D305", "D204", "D318"],
# 	'n')
