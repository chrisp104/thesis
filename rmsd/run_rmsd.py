from rmsd_chain import *
import os

# print os.getcwd()
os.chdir("/Users/Chris/Desktop/Data/D206/2M2D")
rmsdMultiple("5d1q_aligned.pdb", "D206m2d", "E", "output_rmsd.txt")
os.chdir("/Users/Chris/Desktop/Data/D206/2M9D")
rmsdMultiple("5d1q_aligned.pdb", "D206m9d", "E", "output_rmsd.txt")
os.chdir("/Users/Chris/Desktop/Data/D206/2MRD")
rmsdMultiple("5d1q_aligned.pdb", "D206mRd", "E", "output_rmsd.txt")

os.chdir("/Users/Chris/Desktop/Data/D430/2M7D")
rmsdMultiple("5d1x_aligned.pdb", "D430m7d", "E", "output_rmsd.txt")
os.chdir("/Users/Chris/Desktop/Data/D430/2M8D")
rmsdMultiple("5d1x_aligned.pdb", "D430m8d", "E", "output_rmsd.txt")
os.chdir("/Users/Chris/Desktop/Data/D430/2MRD")
rmsdMultiple("5d1x_aligned.pdb", "D430mRd", "E", "output_rmsd.txt")

os.chdir("/Users/Chris/Desktop/Data/D410/2M1D")
rmsdMultiple("5d1z_aligned.pdb", "D410m1d", "I", "output_rmsd.txt")
os.chdir("/Users/Chris/Desktop/Data/D410/2M5D")
rmsdMultiple("5d1z_aligned.pdb", "D410m5d", "I", "output_rmsd.txt")
os.chdir("/Users/Chris/Desktop/Data/D410/2MRD")
rmsdMultiple("5d1z_aligned.pdb", "D410mRd", "I", "output_rmsd.txt")