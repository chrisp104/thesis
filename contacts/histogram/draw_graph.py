import numpy as np
import matplotlib.pyplot as plt
import os

# Draws a histogram for a given total contacts file, returned by all_contacts functions.
#
# 1. 
# drawHistogram()
#
# ARGUMENTS
# 1. file: string - total_contacts file name
# 2. lower: int - lower bound
# 3. upper: int - upper bound
# 4. width: int - with of bin

def drawHistogram(file, lower, upper, width):
	print file

	# BINS
	bins = []
	for i in range(lower, upper, width):
		bins.append(i)

	# FREQUENCIES
	frequencies = []
	f = open(file, 'r')
	lines = f.readlines()

	bin_total = 0
	counter = 0
	for line in lines:
		counter += 1
		
		# add to frequencies if bin width exceeded
		if counter > width:
			frequencies.append(bin_total)
			bin_total = 0
			counter = 1
			
		value = int(line[10:])
		bin_total += value

	frequencies.append(bin_total)
	print frequencies

	# plotting
	pos = np.arange(len(bins))
	width = 1.0

	ax = plt.axes()
	ax.set_xticks(pos+(width / 2))
	ax.set_ylim([0,150])					# MAX Y VALUE HERE
	ax.set_xticklabels(bins, rotation=90, size=5)

	plt.bar(pos, frequencies, width, color='b')

	# save out the file with appropriate name
	if (file[:3] == "mut"):
		plt.savefig(file[15:25]+"_plot.png")
	if (file[:3] == "cry"):
		plt.savefig(file[16:22]+"_plot.png")

	plt.show()
	f.close()



# EXECUTE - NEED TO ONLY DO ONE FOLDER AT A TIME
for fn in os.listdir("/Users/Chris/GitHub/thesis/contacts/histogram/mutants/"):
	if (fn[:3] == "tot"):
		drawHistogram("mutants/"+fn, 119, 457, 10)

for fn in os.listdir("/Users/Chris/GitHub/thesis/contacts/histogram/crystals/"):
	if (fn[:3] == "tot"):
		drawHistogram("crystals/"+fn, 119, 457, 10)


