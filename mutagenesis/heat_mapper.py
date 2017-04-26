import plotly as py
import plotly.graph_objs as go
# heat_mapper.py
#
#	for taking in symmetrical results and making an A x B array of arrays
# turning it into a 2d array and then making heat map from it
#
# ARGUMENTS
# 1. file: string - name of file with results
# 2. x: number of items in each row 
#
# RETURNS
# 	matrix - the matrix of values for heat map
#		labels - the labels corresponding to the order of the matrix entries

def makeArray(file, x):
	HIGH = 200
	
	# MAKING THE MATRIX
	matrix = []
	
	# file
	data = open(file, 'r')
	lines = data.readlines()
	data.close()

	arr = []	# for each row
		
	for i in range(len(lines)):
		if lines[i] == "\n": break

		line = lines[i]
		print line
		num = float(line[7:])
		if num > HIGH:
			arr.append(HIGH)
		else:
			arr.append(num)


		# if filled row, add to matrix and then clear arr
		if (i+1)%x == 0 and i != 0:
			matrix.append(arr)
			arr = []

	print matrix


	# MAKING THE LABEL ARRAY
	labels = []
	for i in range(30):
		labels.append(str(i).zfill(2))

	print labels
	return matrix, labels



matrix, labels = makeArray("/Users/Chris/GitHub/thesis/mutagenesis/results/pairs.txt", 30)

trace = go.Heatmap(
	z=matrix,
	x=labels,
	y=labels)
data=[trace]
py.offline.plot(data, filename='pairs.html')


