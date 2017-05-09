from create_triple import *
import os


for fn in os.listdir("/Users/Chris/GitHub/thesis/mutagenesis/one_ranked_mutations_2"):
	if fn[0] != 'D': continue
	print fn
	createVariant("/Users/Chris/GitHub/thesis/mutagenesis/one_ranked_mutations_2/"+fn, 
	"/Users/Chris/GitHub/thesis/mutagenesis/merged_isdb.pdb", 3, "/Users/Chris/GitHub/thesis/mutagenesis/one_variants_2/"+fn)