from all_contacts import *
import os

#print os.getcwd()

# os.chdir("/Users/Chris/GitHub/thesis/antibodies/D206/2M2D")
# bulkAllContacts("/Users/Chris/GitHub/thesis/contacts/all_bulk_test/", "D206m2d", ['O', 'R', 'T'], 8, normalized=False)
# os.chdir("/Users/Chris/GitHub/thesis/antibodies/D206/2M9D")
# bulkAllContacts("/Users/Chris/GitHub/thesis/contacts/all_bulk_test/", "D206m9d", ['O', 'R', 'T'], 8, normalized=False)
# os.chdir("/Users/Chris/GitHub/thesis/antibodies/D206/2MRD")
# bulkAllContacts("/Users/Chris/GitHub/thesis/contacts/all_bulk_test/", "D206mRd", ['O', 'R', 'T'], 8, normalized=False)
# os.chdir("/Users/Chris/GitHub/thesis/antibodies/D206/MR_m165")
# bulkAllContacts("/Users/Chris/GitHub/thesis/contacts/all_bulk_test/", "D206mut165d", ['O', 'R', 'T'], 8, normalized=False)
# os.chdir("/Users/Chris/GitHub/thesis/antibodies/D206/MR_m390")
# bulkAllContacts("/Users/Chris/GitHub/thesis/contacts/all_bulk_test/", "D206mut390d", ['O', 'R', 'T'], 8, normalized=False)
# os.chdir("/Users/Chris/GitHub/thesis/antibodies/D206/MR_m435")
# bulkAllContacts("/Users/Chris/GitHub/thesis/contacts/all_bulk_test/", "D206mut435d", ['O', 'R', 'T'], 8, normalized=False)

# os.chdir("/Users/Chris/GitHub/thesis/antibodies/D430/2M7D")
# bulkAllContacts("/Users/Chris/GitHub/thesis/contacts/all_bulk_test/", "D430m7d", ['O', 'R', 'T'], 8, normalized=False)
# os.chdir("/Users/Chris/GitHub/thesis/antibodies/D430/2M8D")
# bulkAllContacts("/Users/Chris/GitHub/thesis/contacts/all_bulk_test/", "D430m8d", ['O', 'R', 'T'], 8, normalized=False)
# os.chdir("/Users/Chris/GitHub/thesis/antibodies/D430/2MRD")
# bulkAllContacts("/Users/Chris/GitHub/thesis/contacts/all_bulk_test/", "D430mRd", ['O', 'R', 'T'], 8, normalized=False)
# os.chdir("/Users/Chris/GitHub/thesis/antibodies/D410/MR_m165")
# bulkAllContacts("/Users/Chris/GitHub/thesis/contacts/all_bulk_test/", "D410mut165d", ['O', 'R', 'T'], 8, normalized=False)
# os.chdir("/Users/Chris/GitHub/thesis/antibodies/D410/MR_m390")
# bulkAllContacts("/Users/Chris/GitHub/thesis/contacts/all_bulk_test/", "D410mut390d", ['O', 'R', 'T'], 8, normalized=False)
# os.chdir("/Users/Chris/GitHub/thesis/antibodies/D410/MR_m435")
# bulkAllContacts("/Users/Chris/GitHub/thesis/contacts/all_bulk_test/", "D410mut435d", ['O', 'R', 'T'], 8, normalized=False)

# os.chdir("/Users/Chris/GitHub/thesis/antibodies/D410/2M1D")
# bulkAllContacts("/Users/Chris/GitHub/thesis/contacts/all_bulk_test/", "D410m1d", ['O', 'R', 'T'], 8, normalized=False)
# os.chdir("/Users/Chris/GitHub/thesis/antibodies/D410/2M5D")
# bulkAllContacts("/Users/Chris/GitHub/thesis/contacts/all_bulk_test/", "D410m5d", ['O', 'R', 'T'], 8, normalized=False)
# os.chdir("/Users/Chris/GitHub/thesis/antibodies/D410/2MRD")
# bulkAllContacts("/Users/Chris/GitHub/thesis/contacts/all_bulk_test/", "D410mRd", ['O', 'R', 'T'], 8, normalized=False)
# os.chdir("/Users/Chris/GitHub/thesis/antibodies/D430/MR_m165")
# bulkAllContacts("/Users/Chris/GitHub/thesis/contacts/all_bulk_test/", "D430mut165d", ['O', 'R', 'T'], 8, normalized=False)
# os.chdir("/Users/Chris/GitHub/thesis/antibodies/D430/MR_m390")
# bulkAllContacts("/Users/Chris/GitHub/thesis/contacts/all_bulk_test/", "D430mut390d", ['O', 'R', 'T'], 8, normalized=False)
# os.chdir("/Users/Chris/GitHub/thesis/antibodies/D430/MR_m435")
# bulkAllContacts("/Users/Chris/GitHub/thesis/contacts/all_bulk_test/", "D430mut435d", ['O', 'R', 'T'], 8, normalized=False)


# *********** percentContacts() for all pairs in directory ***********
os.chdir("/Users/Chris/GitHub/thesis/contacts/pairs/")
percentContact(['I', 'E', 'O', 'R', 'T'], 8, "/Users/Chris/GitHub/thesis/contacts/pairs/percent_pairs.txt")