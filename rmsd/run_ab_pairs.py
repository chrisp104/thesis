# to find pairwise rmsds between all 10 models in a given Ab model directory

from rmsd_pairs import *
import os

for d in os.listdir("/Users/Chris/GitHub/thesis/antibodies/"):
	if d[0:1] != 'D': continue

	rmsdPairMultiple("/Users/Chris/GitHub/thesis/antibodies/"+d+"/decoys/", ['H', 'L'],
		"/Users/Chris/GitHub/thesis/antibodies/"+d+"/decoys/rmsd_pairs.txt")