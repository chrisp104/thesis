from res_contacts import *
import os

# ********* FOR EACH AB:AG DIRECTORY **************

# print os.getcwd()
# os.chdir("/Users/Chris/GitHub/thesis/antibodies/D206/2M2D")
# bulkContacts("D206m2d", ['O', 'R', 'T'], 10, "/Users/Chris/GitHub/thesis/contacts/ab_contact_info/")
# os.chdir("/Users/Chris/GitHub/thesis/antibodies/D206/2M9D")
# bulkContacts("D206m9d", ['O', 'R', 'T'], 10, "/Users/Chris/GitHub/thesis/contacts/ab_contact_info/")
# os.chdir("/Users/Chris/GitHub/thesis/antibodies/D206/2MRD")
# bulkContacts("D206mRd", ['O', 'R', 'T'], 10, "/Users/Chris/GitHub/thesis/contacts/ab_contact_info/")

# os.chdir("/Users/Chris/GitHub/thesis/antibodies/D430/2M7D")
# bulkContacts("D430m7d", ['O', 'R', 'T'], 10, "/Users/Chris/GitHub/thesis/contacts/ab_contact_info/")
# os.chdir("/Users/Chris/GitHub/thesis/antibodies/D430/2M8D")
# bulkContacts("D430m8d", ['O', 'R', 'T'], 10, "/Users/Chris/GitHub/thesis/contacts/ab_contact_info/")
# os.chdir("/Users/Chris/GitHub/thesis/antibodies/D430/2MRD")
# bulkContacts("D430mRd", ['O', 'R', 'T'], 10, "/Users/Chris/GitHub/thesis/contacts/ab_contact_info/")

# os.chdir("/Users/Chris/GitHub/thesis/antibodies/D410/2M1D")
# bulkContacts("D410m1d", ['O', 'R', 'T'], 10, "/Users/Chris/GitHub/thesis/contacts/ab_contact_info/")
# os.chdir("/Users/Chris/GitHub/thesis/antibodies/D410/2M5D")
# bulkContacts("D410m5d", ['O', 'R', 'T'], 10, "/Users/Chris/GitHub/thesis/contacts/ab_contact_info/")
# os.chdir("/Users/Chris/GitHub/thesis/antibodies/D410/2MRD")
# bulkContacts("D410mRd", ['O', 'R', 'T'], 10, "/Users/Chris/GitHub/thesis/contacts/ab_contact_info/")



# ************ DIRECTORY BASED CONTACT COUNTING **************

os.chdir("/Users/Chris/GitHub/thesis/contacts/contacts_crystal/")
bulkDirectory(['O', 'R', 'T'], 8, "/Users/Chris/GitHub/thesis/contacts/out_crystal/")

os.chdir("/Users/Chris/GitHub/thesis/contacts/contacts_mutant/")
bulkDirectory(['O', 'R', 'T'], 8, "/Users/Chris/GitHub/thesis/contacts/out_mutant/")