# Script to modify downloaded data in directories for Ab docking

# loop through each Ab directory
for ab in */ ; do
	cd "$ab"

	# IN Ab DIRECTORY
	
	# loop through files
	for f1 in * ; do
		# start at position 1 for length 1
		if [[ ${f1:1:1} == "M" ]]; then

			# remove HEADER from files
			find "$f1" -type f -exec sed -i '' -e "/HEADER/d" {} \;
			# remove END from files
			find "$f1" -type f -exec sed -i '' -e "/END/d" {} \;
			# remove the second to last column from PDB files
			find "$f1" -type f -name '*.pdb' -exec sed -i '' -e 's/^\(.\{72\}\).\{3\}/\1   /' {} \;

			cd "$f1"

			# IN MODEL DIRECTORY
			for f2 in * ; do
				#echo "$f2"
				#model="${f1:2:3}"

				rename "s/model.000./${ab:0:4}m${f1:2:1}d/" *
			done

			cd ..
		fi
	done

	cd ..
done
