# Script to run rmsd.py for all Ab directories and their models

cd ./D206/1M2D
python -c 'import rmsd; rmsd.rmsdMultiple("5d1q_aligned.pdb", "D206m2d", "E", "output.txt")'
cd ../1M9D
python -c 'import rmsd; rmsd.rmsdMultiple("5d1q_aligned.pdb", "D206m9d", "E", "output.txt")'
