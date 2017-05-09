from res_contacts import *
import os

# ***************** FOR EACH AB:AG DIRECTORY ******************

# print os.getcwd()
# os.chdir("/Users/cpark/thesis/antibodies/D206/2M2D")
# bulkContacts("D206m2d", ['O', 'R', 'T'], 10, "/Users/cpark/thesis/contacts/ab_contact_info/")
# os.chdir("/Users/cpark/thesis/antibodies/D206/2M9D")
# bulkContacts("D206m9d", ['O', 'R', 'T'], 10, "/Users/cpark/thesis/contacts/ab_contact_info/")
# os.chdir("/Users/cpark/thesis/antibodies/D206/2MRD")
# bulkContacts("D206mRd", ['O', 'R', 'T'], 10, "/Users/cpark/thesis/contacts/ab_contact_info/")

# os.chdir("/Users/cpark/thesis/antibodies/D430/2M7D")
# bulkContacts("D430m7d", ['O', 'R', 'T'], 10, "/Users/cpark/thesis/contacts/ab_contact_info/")
# os.chdir("/Users/cpark/thesis/antibodies/D430/2M8D")
# bulkContacts("D430m8d", ['O', 'R', 'T'], 10, "/Users/cpark/thesis/contacts/ab_contact_info/")
# os.chdir("/Users/cpark/thesis/antibodies/D430/2MRD")
# bulkContacts("D430mRd", ['O', 'R', 'T'], 10, "/Users/cpark/thesis/contacts/ab_contact_info/")

# os.chdir("/Users/cpark/thesis/antibodies/D410/2M1D")
# bulkContacts("D410m1d", ['O', 'R', 'T'], 10, "/Users/cpark/thesis/contacts/ab_contact_info/")
# os.chdir("/Users/cpark/thesis/antibodies/D410/2M5D")
# bulkContacts("D410m5d", ['O', 'R', 'T'], 10, "/Users/cpark/thesis/contacts/ab_contact_info/")
# os.chdir("/Users/cpark/thesis/antibodies/D410/2MRD")
# bulkContacts("D410mRd", ['O', 'R', 'T'], 10, "/Users/cpark/thesis/contacts/ab_contact_info/")



# ************ DIRECTORY BASED CONTACT COUNTING **************

# os.chdir("/Users/cpark/thesis/contacts/crystal/")
# bulkDirectory(['E', 'I'], 8, "/Users/cpark/thesis/contacts/out_contacts/")

os.chdir("/Users/Chris/GitHub/thesis/antibodies/D206/mR")
bulkDirectory(['O', 'R', 'T'], 8, "/Users/Chris/GitHub/thesis/contacts/out_contacts/")

os.chdir("/Users/Chris/GitHub/thesis/antibodies/D410/mR")
bulkDirectory(['O', 'R', 'T'], 8, "/Users/Chris/GitHub/thesis/contacts/out_contacts/")

os.chdir("/Users/Chris/GitHub/thesis/antibodies/D430/mR")
bulkDirectory(['O', 'R', 'T'], 8, "/Users/Chris/GitHub/thesis/contacts/out_contacts/")


# ********************* PERCENT CORRECT ***********************

# ****** Best RMSD models moved to contacts dir ******

# os.chdir("/Users/cpark/thesis/contacts/D206/")
# percentContact("5d1q.pdb", ['I', 'E', 'O', 'R', 'T'], 8, "/Users/cpark/thesis/contacts/D206/percent_to_crystal.txt")

# os.chdir("/Users/cpark/thesis/contacts/D410/")
# percentContact("5d1z.pdb", ['I', 'E', 'O', 'R', 'T'], 8, "/Users/cpark/thesis/contacts/D410/percent_to_crystal.txt")

# os.chdir("/Users/cpark/thesis/contacts/D430/")
# percentContact("5d1x.pdb", ['I', 'E', 'O', 'R', 'T'], 8, "/Users/cpark/thesis/contacts/D430/percent_to_crystal.txt")


# ****** All models in antibodies' directories ******

# os.chdir("/Users/cpark/thesis/antibodies/D430/3MRD")
# percentContact("5d1x.pdb", ['I', 'E', 'O', 'R', 'T'], 8, "/Users/cpark/thesis/antibodies/D430/3MRD/percent_to_crystal.txt")
# os.chdir("/Users/cpark/thesis/antibodies/D430/MR_m54")
# percentContact("5d1x.pdb", ['I', 'E', 'O', 'R', 'T'], 8, "/Users/cpark/thesis/antibodies/D430/MR_m54/percent_to_crystal.txt")
# os.chdir("/Users/cpark/thesis/antibodies/D430/MR_m165")
# percentContact("5d1x.pdb", ['I', 'E', 'O', 'R', 'T'], 8, "/Users/cpark/thesis/antibodies/D430/MR_m165/percent_to_crystal.txt")
# os.chdir("/Users/cpark/thesis/antibodies/D430/MR_m390")
# percentContact("5d1x.pdb", ['I', 'E', 'O', 'R', 'T'], 8, "/Users/cpark/thesis/antibodies/D430/MR_m390/percent_to_crystal.txt")
# os.chdir("/Users/cpark/thesis/antibodies/D430/MR_m435")
# percentContact("5d1x.pdb", ['I', 'E', 'O', 'R', 'T'], 8, "/Users/cpark/thesis/antibodies/D430/MR_m435/percent_to_crystal.txt")

# os.chdir("/Users/cpark/thesis/antibodies/D410/3MRD")
# percentContact("5d1z.pdb", ['I', 'E', 'O', 'R', 'T'], 8, "/Users/cpark/thesis/antibodies/D410/3MRD/percent_to_crystal.txt")
# os.chdir("/Users/cpark/thesis/antibodies/D410/MR_m52")
# percentContact("5d1z.pdb", ['I', 'E', 'O', 'R', 'T'], 8, "/Users/cpark/thesis/antibodies/D410/MR_m52/percent_to_crystal.txt")
# os.chdir("/Users/cpark/thesis/antibodies/D410/MR_m165")
# percentContact("5d1z.pdb", ['I', 'E', 'O', 'R', 'T'], 8, "/Users/cpark/thesis/antibodies/D410/MR_m165/percent_to_crystal.txt")
# os.chdir("/Users/cpark/thesis/antibodies/D410/MR_m390")
# percentContact("5d1z.pdb", ['I', 'E', 'O', 'R', 'T'], 8, "/Users/cpark/thesis/antibodies/D410/MR_m390/percent_to_crystal.txt")
# os.chdir("/Users/cpark/thesis/antibodies/D410/MR_m435")
# percentContact("5d1z.pdb", ['I', 'E', 'O', 'R', 'T'], 8, "/Users/cpark/thesis/antibodies/D410/MR_m435/percent_to_crystal.txt")

# os.chdir("/Users/cpark/thesis/antibodies/D206/3MRD")
# percentContact("5d1q.pdb", ['I', 'E', 'O', 'R', 'T'], 8, "/Users/cpark/thesis/antibodies/D206/3MRD/percent_to_crystal.txt")
# os.chdir("/Users/cpark/thesis/antibodies/D206/MR_m54")
# percentContact("5d1q.pdb", ['I', 'E', 'O', 'R', 'T'], 8, "/Users/cpark/thesis/antibodies/D206/MR_m54/percent_to_crystal.txt")
# os.chdir("/Users/cpark/thesis/antibodies/D206/MR_m165")
# percentContact("5d1q.pdb", ['I', 'E', 'O', 'R', 'T'], 8, "/Users/cpark/thesis/antibodies/D206/MR_m165/percent_to_crystal.txt")
# os.chdir("/Users/cpark/thesis/antibodies/D206/MR_m390")
# percentContact("5d1q.pdb", ['I', 'E', 'O', 'R', 'T'], 8, "/Users/cpark/thesis/antibodies/D206/MR_m390/percent_to_crystal.txt")
# os.chdir("/Users/cpark/thesis/antibodies/D206/MR_m435")
# percentContact("5d1q.pdb", ['I', 'E', 'O', 'R', 'T'], 8, "/Users/cpark/thesis/antibodies/D206/MR_m435/percent_to_crystal.txt")






