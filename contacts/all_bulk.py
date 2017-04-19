from all_contacts import *
import os

#print os.getcwd()

os.chdir("/Users/Chris/GitHub/thesis/antibodies/D206/m2")
bulkAllContacts("/Users/Chris/GitHub/thesis/contacts/out_contacts_2/", ['O', 'R', 'T'], 8, normalized=False)
# os.chdir("/Users/Chris/GitHub/thesis/antibodies/D206/m9")
# bulkAllContacts("/Users/Chris/GitHub/thesis/contacts/all_bulk_test/", "D206m9d", ['O', 'R', 'T'], 8, normalized=False)
# os.chdir("/Users/Chris/GitHub/thesis/antibodies/D206/mR")
# bulkAllContacts("/Users/Chris/GitHub/thesis/contacts/all_bulk_test/", "D206mRd", ['O', 'R', 'T'], 8, normalized=False)
# os.chdir("/Users/Chris/GitHub/thesis/antibodies/D206/mR_m165")
# bulkAllContacts("/Users/Chris/GitHub/thesis/contacts/all_bulk_test/", "D206mut165d", ['O', 'R', 'T'], 8, normalized=False)
# os.chdir("/Users/Chris/GitHub/thesis/antibodies/D206/mR_m390")
# bulkAllContacts("/Users/Chris/GitHub/thesis/contacts/all_bulk_test/", "D206mut390d", ['O', 'R', 'T'], 8, normalized=False)
# os.chdir("/Users/Chris/GitHub/thesis/antibodies/D206/mR_m435")
# bulkAllContacts("/Users/Chris/GitHub/thesis/contacts/all_bulk_test/", "D206mut435d", ['O', 'R', 'T'], 8, normalized=False)

# os.chdir("/Users/Chris/GitHub/thesis/antibodies/D430/m7")
# bulkAllContacts("/Users/Chris/GitHub/thesis/contacts/all_bulk_test/", "D430m7d", ['O', 'R', 'T'], 8, normalized=False)
# os.chdir("/Users/Chris/GitHub/thesis/antibodies/D430/m8")
# bulkAllContacts("/Users/Chris/GitHub/thesis/contacts/all_bulk_test/", "D430m8d", ['O', 'R', 'T'], 8, normalized=False)
# os.chdir("/Users/Chris/GitHub/thesis/antibodies/D430/mR")
# bulkAllContacts("/Users/Chris/GitHub/thesis/contacts/all_bulk_test/", "D430mRd", ['O', 'R', 'T'], 8, normalized=False)
# os.chdir("/Users/Chris/GitHub/thesis/antibodies/D410/mR_m165")
# bulkAllContacts("/Users/Chris/GitHub/thesis/contacts/all_bulk_test/", "D410mut165d", ['O', 'R', 'T'], 8, normalized=False)
# os.chdir("/Users/Chris/GitHub/thesis/antibodies/D410/mR_m390")
# bulkAllContacts("/Users/Chris/GitHub/thesis/contacts/all_bulk_test/", "D410mut390d", ['O', 'R', 'T'], 8, normalized=False)
# os.chdir("/Users/Chris/GitHub/thesis/antibodies/D410/mR_m435")
# bulkAllContacts("/Users/Chris/GitHub/thesis/contacts/all_bulk_test/", "D410mut435d", ['O', 'R', 'T'], 8, normalized=False)

# os.chdir("/Users/Chris/GitHub/thesis/antibodies/D410/m1")
# bulkAllContacts("/Users/Chris/GitHub/thesis/contacts/all_bulk_test/", "D410m1d", ['O', 'R', 'T'], 8, normalized=False)
# os.chdir("/Users/Chris/GitHub/thesis/antibodies/D410/m5")
# bulkAllContacts("/Users/Chris/GitHub/thesis/contacts/all_bulk_test/", "D410m5d", ['O', 'R', 'T'], 8, normalized=False)
# os.chdir("/Users/Chris/GitHub/thesis/antibodies/D410/mR")
# bulkAllContacts("/Users/Chris/GitHub/thesis/contacts/all_bulk_test/", "D410mRd", ['O', 'R', 'T'], 8, normalized=False)
# os.chdir("/Users/Chris/GitHub/thesis/antibodies/D430/mR_m165")
# bulkAllContacts("/Users/Chris/GitHub/thesis/contacts/all_bulk_test/", "D430mut165d", ['O', 'R', 'T'], 8, normalized=False)
# os.chdir("/Users/Chris/GitHub/thesis/antibodies/D430/mR_m390")
# bulkAllContacts("/Users/Chris/GitHub/thesis/contacts/all_bulk_test/", "D430mut390d", ['O', 'R', 'T'], 8, normalized=False)
# os.chdir("/Users/Chris/GitHub/thesis/antibodies/D430/mR_m435")
# bulkAllContacts("/Users/Chris/GitHub/thesis/contacts/all_bulk_test/", "D430mut435d", ['O', 'R', 'T'], 8, normalized=False)


# *********** percentContacts() for all pairs in directory ***********
# os.chdir("/Users/Chris/GitHub/thesis/contacts/pairs/")
# percentContact(['I', 'E', 'O', 'R', 'T'], 8, "/Users/Chris/GitHub/thesis/contacts/pairs/percent_pairs.txt")