import os
from collections import OrderedDict

# scripts to take output of mutate.py and identify mutations that
# disrupt binding for all docking models
#
# FUNCTIONS
# 1. rankForAb - orders Ag mutations by prevalence among docking models for a given Ab
# 2. rankForAll - orders Ag mutations by prevalence among docking models for all Abs
# 3. analyzePairs - take output from function 1 and count for each docking model pair how many similar
# disruptive mutations they have in common





# 1.
# rankForAb()
#
# take _mutations.txt output files for all 30 docking models, loop through them
# and rank the most prevalent mutations
# 
# ARGUMENTS
# 1. directory: str - the directory containing the 30 _mutations.txt files
# 2. out_path: str - directory
# 3. num_affected: int - the cutoff for number of affected docking models below which 
#		 a mutation is not counted
# 4. cutoff_score: float - the cumulative disruption score below which a mutation is not considered
# 5. exclusions: array - containing file names to exclude when ranking since they are invalid docking models
#
# RETURNS 
# 	dictionary:
#			key: tuple = (Ag res #, mutation)
#			value: array = [[models], # models, total score)]
def rankForAb(directory, out_path, num_affected, cutoff_score, exclusions=[]):
	out_dir = out_path[:-10]
	if not os.path.exists(out_dir):
		os.makedirs(out_dir)

	out = open(out_path, 'w')
	ranked_mutations = {}

	for fn in os.listdir(directory):
		if fn[-13:] != "mutations.txt" and fn[0] != 'D': continue
		
		# EXCLUSIONS
		if fn[:9] in exclusions:
			continue

		print fn

		file = open(directory+fn, 'r')
		lines = file.readlines()
		file.close()

		# for each mutations.txt file
		for line in lines:
			if len(line) == 6: 
				continue

			line_split = line.split(':')
			ag_num = line_split[0]
			mutations = line_split[1].strip().split(',')

			# loop through mutations
			for mutation in mutations:
				split_mutation = mutation.split('|')
				mut_aa = split_mutation[0]
				mut_score = str(split_mutation[1])

				key = (ag_num, mut_aa)

				# dictionary stuff

				# if the mutation has negative score for this model, don't add and make note
				if float(mut_score) < 0:
					ranked_mutations[key] = [0, 0, 0]
					continue
				# if this mutation was not disruptive for another model, then don't bother
				elif key in ranked_mutations and ranked_mutations[key][1] == 0:
					continue
				elif key in ranked_mutations:
					ranked_mutations[key][0].append(fn[-17:-15])
					ranked_mutations[key][1] += 1
					ranked_mutations[key][2] += float(mut_score)
				else:
					ranked_mutations[key] = [[fn[-17:-15]], 1, float(mut_score)]

	# sort the dictionary by # of models each mutation is disruptive in
	ranked_mutations = sorted(ranked_mutations.items(), key=lambda e: e[1][1], reverse=True)


	# write out results to a file
	for mutation in ranked_mutations:
		ag = mutation[0]
		ab = mutation[1]

		# if this mutation was not disruptive then don't print
		if ab[1] == 0:
			continue

		models = ab[0]
		num_models = ab[1]
		score = ab[2]

		# *** SKIP SWITCHES

		# skip if mutation only affects N models
		if num_models < num_affected:
			continue

		# skip if mutation score is not above N
		if score < cutoff_score:
			continue

		out.write(ag[0] + "m" + ag[1]+': ')
		
		# write the models
		first = True
		for model in models:
			if first:
				out.write(model)
				first = False
				continue
			out.write(','+model)
		
		# write num_models
		out.write('|'+str(num_models))

		# write the total score over all models for this mutation
		out.write('|'+str(score)+"\n")


	# return and terminate
	out.close()
	return ranked_mutations







# 2.
# rankForAll()
#
# take _mutations.txt output files for all 30 docking models for all Abs, loop through them
# and rank the most prevalent mutations
# 
# ARGUMENTS
# 1. directory: str - the directory containing the Ab directories containing the 30 _mutations.txt files
# 2. out_path: str - directory
# 3. num_affected: int - the cutoff for number of affected docking models below which 
#		 a mutation is not counted
# 4. cutoff_score: float - the cumulative disruption score below which a mutation is not considered
# 5. exclusions: array - containing file names to exclude when ranking since they are invalid docking models
# 6. type: char - 'm' or 'n'
#
# RETURNS 
# 	dictionary:
#			key: tuple = (Ag res #, mutation)
#			value: array = [[models], # models, total score)]
def rankForAll(directory, out_path, num_affected, cutoff_score, exclusions=[], modelType='m'):
	out = open(out_path, 'w')
	ranked_mutations = {}
	for antibody in os.listdir(directory):
		if antibody[0] != 'D': continue
		for fn in os.listdir(directory+antibody):
			if fn[0] != modelType: continue
			for fn2 in os.listdir(directory+antibody+'/'+fn):
				if fn2[:9] in exclusions: continue
				#print fn2
				if fn2[-13:] != "mutations.txt" and fn2[0] != 'D': continue
			
				file = open(directory+antibody+'/'+fn+'/'+fn2, 'r')
				lines = file.readlines()
				file.close()

				# for each mutations.txt file
				for line in lines:
					if len(line) == 6: 
						continue

					line_split = line.split(':')
					ag_num = line_split[0]
					mutations = line_split[1].strip().split(',')

					# loop through mutations
					for mutation in mutations:
						split_mutation = mutation.split('|')
						mut_aa = split_mutation[0]
						mut_score = str(split_mutation[1])

						key = (ag_num, mut_aa)

						# add to the dictionary
						if key in ranked_mutations:
							ranked_mutations[key][0].append(fn[-17:-15])
							ranked_mutations[key][1] += 1
							ranked_mutations[key][2] += float(mut_score)
						else:
							ranked_mutations[key] = [[fn[-17:-15]], 1, float(mut_score)]

	# sort the dictionary by # of models each mutation is disruptive in
	ranked_mutations = sorted(ranked_mutations.items(), key=lambda e: e[1][1], reverse=True)


	# write out results to a file
	for mutation in ranked_mutations:
		ag = mutation[0]
		ab = mutation[1]

		# if this mutation was not disruptive then don't print
		if ab[1] == 0:
			continue

		models = ab[0]
		num_models = ab[1]
		score = ab[2]

		# *** SKIP SWITCHES

		# skip if mutation only affects N models
		if num_models <= num_affected:
			continue

		# skip if mutation score is not above N
		if score < cutoff_score:
			continue

		out.write(ag[0] + "m" + ag[1]+': ')
		
		# write num_models
		out.write(str(num_models))

		# write the total score over all models for this mutation
		out.write('|'+str(score)+"\n")


	# return and terminate
	out.close()
	return ranked_mutations








# 2.
# analyzePairs()
#
# take output from function 1 and count for each docking model pair how many similar
# disruptive mutations they have in common
# 
# ARGUMENTS
# 1. data: str - path of output file from function1
# 2. out_path: str - path of output file 
#
# RETURNS 
# 	dictionary:
#			key: str - "Docking model #: docking model #"
#			value: int - number of mutations in common
def analyzePairs(data, out_path):
	out = open(out_path, 'w')
	file = open(data, 'r')
	lines = file.readlines()
	file.close()

	#initialize dictionary
	pairs = {}
	for i in range(30):
		for j in range(i+1, 30):
			key = str(i).zfill(2)+':'+str(j).zfill(2)
			pairs[key] = 0

	# go through lines and count
	for line in lines:
		models = line.split(':')[1].strip().split('|')[0]
		model_list = models.split(',')
		
		for i in range(len(model_list)):
			for j in range(i+1, len(model_list)):
				m1 = model_list[i]
				m2 = model_list[j]
				if m1 < m2:
					key = m1+':'+m2
				else:
					key = m2+':'+m1

				pairs[key] += 1

	keys = sorted(pairs)

	for key in keys:
		# if pairs[key] == 0: continue
		out.write(key+": "+str(pairs[key])+"\n")

	return pairs

# data = rankMutations("/Users/Chris/GitHub/thesis/mutagenesis/mutations/D102/m0/", "/Users/Chris/GitHub/thesis/mutagenesis/mutations/D102/m0/ranked.txt")
# analyzePairs("/Users/Chris/GitHub/thesis/mutagenesis/mutations/D102/m0/ranked.txt", "/Users/Chris/GitHub/thesis/mutagenesis/mutations/D102/m0/pairs.txt")









