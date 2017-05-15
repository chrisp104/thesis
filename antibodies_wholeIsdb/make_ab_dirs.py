import os

file = open("abNames.txt", 'r')
antibodies = file.readlines()
for antibody in antibodies:
	if antibody[0] != "#":
		name = antibody[0:2]+antibody[3:5]
		os.makedirs(name)

file.close()