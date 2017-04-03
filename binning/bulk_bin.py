from ag_binner import *
import os

# print os.getcwd()
os.chdir("/Users/Chris/Desktop/Data/D206/2M2D")
find_bin("D206m2d", ['O', 'R', 'T'], 8, "output_bin.txt")
os.chdir("/Users/Chris/Desktop/Data/D206/2M9D")
find_bin("D206m9d", ['O', 'R', 'T'], 8, "output_bin.txt")
os.chdir("/Users/Chris/Desktop/Data/D206/2MRD")
find_bin("D206mRd", ['O', 'R', 'T'], 8, "output_bin.txt")

os.chdir("/Users/Chris/Desktop/Data/D430/2M7D")
find_bin("D430m7d", ['O', 'R', 'T'], 8, "output_bin.txt")
os.chdir("/Users/Chris/Desktop/Data/D430/2M8D")
find_bin("D430m8d", ['O', 'R', 'T'], 8, "output_bin.txt")
os.chdir("/Users/Chris/Desktop/Data/D430/2MRD")
find_bin("D430mRd", ['O', 'R', 'T'], 8, "output_bin.txt")

os.chdir("/Users/Chris/Desktop/Data/D410/2M1D")
find_bin("D410m1d", ['O', 'R', 'T'], 8, "output_bin.txt")
os.chdir("/Users/Chris/Desktop/Data/D410/2M5D")
find_bin("D410m5d", ['O', 'R', 'T'], 8, "output_bin.txt")
os.chdir("/Users/Chris/Desktop/Data/D410/2MRD")
find_bin("D410mRd", ['O', 'R', 'T'], 8, "output_bin.txt")