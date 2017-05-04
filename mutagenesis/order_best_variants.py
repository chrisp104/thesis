# take the best variant set output file and order them according to expected bins
import os

bestVariants = open("/Users/Chris/GitHub/thesis/mutagenesis/bestVariants_3.txt", 'r')
order = open("/Users/Chris/GitHub/thesis/mutagenesis/ab_order.txt", 'r')

out = open("/Users/Chris/GitHub/thesis/mutagenesis/bestVariants_3_ordered.txt", 'w')

data = []
for line in bestVariants.readlines():
	if line == '\n': continue
	data.append(line)

for ab in order.readlines():
	ab = ab[:4]
	print ab
	for datum in data:
		if ab == datum[:4]:
			out.write(datum)
			

bestVariants.close()
order.close()
out.close()