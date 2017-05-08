# to count the number of mutations needed to cover all models
import os

high = 0
for fn in os.listdir("/Users/Chris/GitHub/thesis/mutagenesis/ranked_mutations/"):
	if fn[0] != 'D': continue
	file = open("/Users/Chris/GitHub/thesis/mutagenesis/ranked_mutations/"+fn, 'r')
	lines = file.readlines()
	file.close()

	counting = []
	i = 0
	for line in lines:
		i += 1
		models = line.split(':')[1].strip().split('|')[0].split(',')
		for model in models:
			if not model in counting:
				counting.append(model)

		if len(counting) == 30:
			if i > high:
				high = i
			break

	if len(counting) != 30:
		print fn

print "HIGH: "+str(high)