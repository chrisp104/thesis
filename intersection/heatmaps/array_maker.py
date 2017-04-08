# array_maker.py
#
#	for taking in symmetrical results and making an A x B array of arrays
#
# ARGUMENTS
# 1. file: string - name of file with results
# 2. x: number of items in each row 

def makeArray(file, x):
	matrix = []
	
	# file
	data = open(file, 'r')
	lines = data.readlines()
	data.close()

	arr = []	# for each row
		
	for i in range(len(lines)):

		line = lines[i]
		num = int(line[17:19])
		arr.append(num)

		# if filled row, add to matrix and then clear arr
		if (i+1)%x == 0 and i != 0:
			matrix.append(arr)
			arr = []

	#matrix.append(arr)

	print matrix

makeArray("output.txt", 9)

