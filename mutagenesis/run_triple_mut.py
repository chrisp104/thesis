from triple_mutation import *
import os


for fn in os.listdir("/Users/Chris/GitHub/thesis/mutagenesis/ranked_mutations"):
	if fn[0] != 'D': continue
	print fn
	createVariant("/Users/Chris/GitHub/thesis/mutagenesis/ranked_mutations/"+fn, 
	"/Users/Chris/GitHub/thesis/mutagenesis/merged_isdb.pdb", 3, "/Users/Chris/GitHub/thesis/mutagenesis/variants/"+fn)