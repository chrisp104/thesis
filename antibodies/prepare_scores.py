# script to pull out total scores for each of the 10 Ab models

import os


for directory in os.listdir("./"):
	if directory[0] != 'D': continue
	print directory

	out = open("./"+directory+"/decoys/scores.txt", 'w')

	try:
		file = open("./"+directory+"/decoys/score.fasc", 'r')
	except IOError:
		continue
	scores = file.readlines()
	file.close()

	score_dict = {}

	for score in scores:
		info = score.split()
		# put into dictionary so we can sort in output
		if info[0] == "SCORE:" and info[1] != "total_score":
			score_dict[float(info[1])] = info[43]

	for key in sorted(score_dict):		
		out.write(score_dict[key]+": "+str(key)+"\n")

	out.close()