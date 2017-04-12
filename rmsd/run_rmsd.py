from rmsd import *
import os

print os.getcwd()
os.chdir("/Users/Chris/GitHub/thesis/antibodies/D206/2M2D")
rmsdMultiple("5d1q_aligned.pdb", "D206m2d", "E", ['O', 'R', 'T'], "/Users/Chris/GitHub/thesis/rmsd/outputs/rmsd_D206m2.txt")
os.chdir("/Users/Chris/GitHub/thesis/antibodies/D206/2M9D")
rmsdMultiple("5d1q_aligned.pdb", "D206m9d", "E", ['O', 'R', 'T'], "/Users/Chris/GitHub/thesis/rmsd/outputs/rmsd_D206m9.txt")
os.chdir("/Users/Chris/GitHub/thesis/antibodies/D206/2MRD")
rmsdMultiple("5d1q_aligned.pdb", "D206mRd", "E", ['O', 'R', 'T'], "/Users/Chris/GitHub/thesis/rmsd/outputs/rmsd_D206mR.txt")
os.chdir("/Users/Chris/GitHub/thesis/antibodies/D206/MR_m165")
rmsdMultiple("5d1q_aligned.pdb", "D206mut165d", "E", ['O', 'R', 'T'], "/Users/Chris/GitHub/thesis/rmsd/outputs/rmsd_D206mut165.txt")
os.chdir("/Users/Chris/GitHub/thesis/antibodies/D206/MR_m390")
rmsdMultiple("5d1q_aligned.pdb", "D206mut390d", "E", ['O', 'R', 'T'], "/Users/Chris/GitHub/thesis/rmsd/outputs/rmsd_D206mut390.txt")
os.chdir("/Users/Chris/GitHub/thesis/antibodies/D206/MR_m435")
rmsdMultiple("5d1q_aligned.pdb", "D206mut435d", "E", ['O', 'R', 'T'], "/Users/Chris/GitHub/thesis/rmsd/outputs/rmsd_D206mut435.txt")


os.chdir("/Users/Chris/GitHub/thesis/antibodies/D430/2M7D")
rmsdMultiple("5d1x_aligned.pdb", "D430m7d", "E", ['O', 'R', 'T'], "/Users/Chris/GitHub/thesis/rmsd/outputs/rmsd_D430m7.txt")
os.chdir("/Users/Chris/GitHub/thesis/antibodies/D430/2M8D")
rmsdMultiple("5d1x_aligned.pdb", "D430m8d", "E", ['O', 'R', 'T'], "/Users/Chris/GitHub/thesis/rmsd/outputs/rmsd_D430m8.txt")
os.chdir("/Users/Chris/GitHub/thesis/antibodies/D430/2MRD")
rmsdMultiple("5d1x_aligned.pdb", "D430mRd", "E", ['O', 'R', 'T'], "/Users/Chris/GitHub/thesis/rmsd/outputs/rmsd_D430mR.txt")
os.chdir("/Users/Chris/GitHub/thesis/antibodies/D430/MR_m165")
rmsdMultiple("5d1x_aligned.pdb", "D430mut165d", "E", ['O', 'R', 'T'], "/Users/Chris/GitHub/thesis/rmsd/outputs/rmsd_D430mut165.txt")
os.chdir("/Users/Chris/GitHub/thesis/antibodies/D430/MR_m390")
rmsdMultiple("5d1x_aligned.pdb", "D430mut390d", "E", ['O', 'R', 'T'], "/Users/Chris/GitHub/thesis/rmsd/outputs/rmsd_D430mut390.txt")
os.chdir("/Users/Chris/GitHub/thesis/antibodies/D430/MR_m435")
rmsdMultiple("5d1x_aligned.pdb", "D430mut435d", "E", ['O', 'R', 'T'], "/Users/Chris/GitHub/thesis/rmsd/outputs/rmsd_D430mut435.txt")


os.chdir("/Users/Chris/GitHub/thesis/antibodies/D410/2M1D")
rmsdMultiple("5d1z_aligned.pdb", "D410m1d", "I", ['O', 'R', 'T'], "/Users/Chris/GitHub/thesis/rmsd/outputs/rmsd_D410m1.txt")
os.chdir("/Users/Chris/GitHub/thesis/antibodies/D410/2M5D")
rmsdMultiple("5d1z_aligned.pdb", "D410m5d", "I", ['O', 'R', 'T'], "/Users/Chris/GitHub/thesis/rmsd/outputs/rmsd_D410m5.txt")
os.chdir("/Users/Chris/GitHub/thesis/antibodies/D410/2MRD")
rmsdMultiple("5d1z_aligned.pdb", "D410mRd", "I", ['O', 'R', 'T'], "/Users/Chris/GitHub/thesis/rmsd/outputs/rmsd_D410mR.txt")
os.chdir("/Users/Chris/GitHub/thesis/antibodies/D410/MR_m165")
rmsdMultiple("5d1z_aligned.pdb", "D410mut165d", "I", ['O', 'R', 'T'], "/Users/Chris/GitHub/thesis/rmsd/outputs/rmsd_D410mut165.txt")
os.chdir("/Users/Chris/GitHub/thesis/antibodies/D410/MR_m390")
rmsdMultiple("5d1z_aligned.pdb", "D410mut390d", "I", ['O', 'R', 'T'], "/Users/Chris/GitHub/thesis/rmsd/outputs/rmsd_D410mut390.txt")
os.chdir("/Users/Chris/GitHub/thesis/antibodies/D410/MR_m435")
rmsdMultiple("5d1z_aligned.pdb", "D410mut435d", "I", ['O', 'R', 'T'], "/Users/Chris/GitHub/thesis/rmsd/outputs/rmsd_D410mut435.txt")


# THESE WERE THE CRYSTAL STRUCTURES DOCKED ONCE MORE TO SEE IF BETTER PREDICTIONS WOULD HAPPEN BUT THEY ARE THE EXACT SAME
# os.chdir("/Users/Chris/GitHub/thesis/antibodies/D206/3MRD")
# rmsdMultiple("5d1q_aligned.pdb", "D206mRd", "E", ['O', 'R', 'T'], "/Users/Chris/GitHub/thesis/rmsd/outputs/rmsd_D206mR_2.txt")
# os.chdir("/Users/Chris/GitHub/thesis/antibodies/D430/3MRD")
# rmsdMultiple("5d1x_aligned.pdb", "D430mRd", "E", ['O', 'R', 'T'], "/Users/Chris/GitHub/thesis/rmsd/outputs/rmsd_D430mR_2.txt")
# os.chdir("/Users/Chris/GitHub/thesis/antibodies/D410/3MRD")
# rmsdMultiple("5d1z_aligned.pdb", "D410mRd", "I", ['O', 'R', 'T'], "/Users/Chris/GitHub/thesis/rmsd/outputs/rmsd_D410mR_2.txt")