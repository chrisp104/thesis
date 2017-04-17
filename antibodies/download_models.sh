# Script to SCP decoys from remote to local directory

# loop through each Ab directory
for ab in */ ; do
	echo "$ab"
	# scp -rp cpark@anthill.cs.dartmouth.edu:/home/anthill/cpark/yeung/${ab}/decoys ./${ab}
	scp -rp cpark@anthill.cs.dartmouth.edu:/home/anthill/cpark/yeung/${ab}/score.fasc ./${ab}/decoys
done