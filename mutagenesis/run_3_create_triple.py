from create_triple import *
import os


for fn in os.listdir("/Users/Chris/GitHub/thesis/mutagenesis/one_ranked_mutations"):
	if fn[0] != 'D': continue
	print fn
	createVariant("/Users/Chris/GitHub/thesis/mutagenesis/one_ranked_mutations/"+fn, 
	"/Users/Chris/GitHub/thesis/mutagenesis/merged_isdb.pdb", 3, "/Users/Chris/GitHub/thesis/mutagenesis/one_variants/"+fn)