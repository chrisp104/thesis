from create_variants import *
import os


for fn in os.listdir("/Users/Chris/GitHub/thesis/mutagenesis/run_1/ranked_mutations"):
	if fn[0] != 'D': continue
	print fn
	createVariant("/Users/Chris/GitHub/thesis/mutagenesis/run_1/ranked_mutations/"+fn, 
	"/Users/Chris/GitHub/thesis/mutagenesis/merged_isdb.pdb", 2, "/Users/Chris/GitHub/thesis/mutagenesis/run_1/variants_2/"+fn)