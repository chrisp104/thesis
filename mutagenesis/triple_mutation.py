# scripts to take ranked mutation output files and create triple mutation variants in which
# the mutations are within 12 angstroms of each other
#
# FUNCTIONS
# 1. createVariant



# 1.
# createVariant()
#
# see above
# 
# ARGUMENTS
# 1. ranked_file: str - the path of the ranked mutation file from which to create variants
# 2. isdb: str - path of isdb protein pdb file
# 3. k: int - the number of mutations per variant
# 4. out_file: str - path to write variants to
#
# RETURNS 
# 	array of arrays, each array representing a variant containing k mutations in the form 114mGLU
def createVariant(ranked_file, isdb, k, out_file):
	