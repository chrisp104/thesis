# script to take the mutation profiles from the three crystal structures and write the 
# "experimentally disruptive" mutations to a file
import os

out = open("/Users/Chris/GitHub/thesis/mutagenesis/confirmed_mutations.txt", 'a')
os.chdir("/Users/Chris/GitHub/thesis/contacts/crystals/")
for fn in os.listdir("/Users/Chris/GitHub/thesis/contacts/crystals/"):
	if fn[-13:] != 'mutations.txt': continue

	ab = ''
	if fn == "D5d1q_mutations.txt": ab = "D206"
	if fn == "D5d1x_mutations.txt": ab = "D430"
	if fn == "D5d1z_mutations.txt": ab = "D410"
	print fn, ab


	file = open(fn, 'r')
	for line in file.readlines():
		if len(line) <= 6: continue

		ag_res = line[:3]
		data = line.strip().split(' ')[1]
		mutations = data.split(',')

		for mutation in mutations:
			msplit = mutation.split('|')
			res = msplit[0]
			score = msplit[1]

			if float(score) > 0.5:	#********************* CHANGE HERE
				out.write(ag_res+'m'+res+':'+ab+'-0\n')

	file.close()

	