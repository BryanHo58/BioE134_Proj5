import os
from Bio import Entrez

"""
Various functions to collect new set of source data for project 5.

Remember to change filepaths for personal directory
I've listed every filepath with a "# CHANGE_ME" to make it easier to ctrl-F

@author Bryan Ho
"""



# Takes in text file of tab-seperated values and parses them to obtain the names of each chemical
# Returns a list of chemical names
def initiate():
    Ecoli_symless = 'C:/Users/bryan/Documents/Berkeley/4thYear/BioE134/BioE134_Proj5/E._coli_W_Chemicals_SymLess.txt'       # CHANGE_ME

    chemicals_list = []
    hundred = 1
    with open(Ecoli_symless , 'r', encoding="utf-8") as Ecoli_data:
        Ecoli_lines = Ecoli_data.readlines()
    for line in Ecoli_lines:
        data = line.split("\t")
        if data[2] == '\n':
            chemicals_list.append(data[1].strip("\""))
    print(len(chemicals_list))

    return chemicals_list

# Takes in text file of names and list of PubChem ids and parses them to a dictionary
# Returns a dictionary object with the key = name and values = list of PubChem ids
def initiate2():
    PubChem_ID = 'C:/Users/bryan/Documents/Berkeley/4thYear/BioE134/BioE134_Proj5/E.coli_ID_Pubchem.txt'     # CHANGE_ME

    ID_lib = {}
    with open(PubChem_ID , 'r', encoding="utf-8") as PCID:
        PCID_lines = PCID.readlines()
    for line in PCID_lines:
        name, IDs = line.split("\t")
        if IDs != "\n":
            ID_list = IDs.strip("\n").split(",")
            ID_lib[name] = ID_list

    return ID_lib

Entrez.email = "bryanho58@berkeley.edu"  # Always tell NCBI who you are     # CHANGE_ME

# Takes in a dictionary object with the key = name and values = list of PubChem ids and generates a text file of only PubChem IDs
# Outputs a "best" text file consisting of only PubChem IDs to be processed by the online PubChem Download Service: https://pubchem.ncbi.nlm.nih.gov/pc_fetch/pc_fetch.cgi
# Outputs a "labeled" text file representing the dictionary object for future use
def obtain_SMILES(ID_lib):
    text = ""
    text2 = ""
    for name in ID_lib.keys():
        ID_list = ID_lib[name]
        text += name + "\t" + ",".join(ID_list) + "\n"
        text2 += "\n".join(ID_list) + "\n"

    best = 'C:/Users/bryan/Documents/Berkeley/4thYear/BioE134/BioE134_Proj5/E.coli_ID_Pubchem_best.txt'   # CHANGE_ME
    out_best = open(best, "w")
    out_best.write(text2)
    out_best.close()
    labeled = 'C:/Users/bryan/Documents/Berkeley/4thYear/BioE134/BioE134_Proj5/E.coli_ID_Pubchem_labeled.txt' # CHANGE_ME
    out_labeled = open(labeled, "w")
    out_labeled.write(text)
    out_labeled.close()

# Takes in three files:
# "labeled" contains name search query to list of representative PubChem IDs
# "SMILES" contains the ID to SMILES data obtained from the online PubChem Download Service
# "InCHI" contains the ID to SMILES data obtained from the online PubChem Download Service
# Organizes the three files into a master document
# Outputs a master document file with format: "Name/Search query    PubChem ID  SMILES  InCl"
def combine_master_doc():
    labeled = "C:/Users/bryan/Documents/Berkeley/4thYear/BioE134/BioE134_Proj5/E.coli_ID_Pubchem_labeled.txt" # CHANGE_ME
    best_SMILES = "C:/Users/bryan/Documents/Berkeley/4thYear/BioE134/BioE134_Proj5/E.coli_ID_Pubchem_best_SMILES.txt" # CHANGE_ME
    best_InCHI = "C:/Users/bryan/Documents/Berkeley/4thYear/BioE134/BioE134_Proj5/E.coli_ID_Pubchem_best_InCHI.txt"   # CHANGE_ME

    text = "Name/Search Query\tPubChem ID\tSMILES\tInChl\n"

    ID_lib = {}
    with open(labeled , 'r', encoding="utf-8") as labels:
        labels_lines = labels.readlines()
    for line in labels_lines:
        name, IDs = line.strip("\n").split("\t")
        if IDs != "\n":
            ID_list = IDs.split(",")
            for id in ID_list:
                if ID_lib.get(id, -1) == -1:
                    ID_lib[id] = [name]
                else:
                    ID_lib[id].append(name)

    with open(best_SMILES , 'r', encoding="utf-8") as b_smiles:
        SMILES_lines = b_smiles.readlines()
    with open(best_InCHI , 'r', encoding="utf-8") as b_inchi:
        InCHI_lines = b_inchi.readlines()
    for index in range(len(SMILES_lines)):
        SMILES_line = SMILES_lines[index].strip("\n")
        InCHI_line = InCHI_lines[index].strip("\n")

        id, smiles = SMILES_line.split("\t")
        id, inchi = InCHI_line.split("\t")

        for name in ID_lib[id]:
            text += name + "\t" + id + "\t" + smiles + "\t" + inchi + "\n"
            # smiles = "" #use to save space, removes repeat data
            # inchi = ""  #use to save space, removes repeat data

    PubChem_out = 'C:/Users/bryan/Documents/Berkeley/4thYear/BioE134/BioE134_Proj5/E.coli_W_Pubchem.txt' # CHANGE_ME
    out_PubChem = open(PubChem_out, "w")
    out_PubChem.write(text)
    out_PubChem.close()

# Takes in list of names and searches them in the PubChem database
# Outputs a text file with the name and all PubChem IDs that returned from search
# Text file format = name \t comma-separated list of ids \n
def obtain_uids(master_list):
    count = 0
    text = ""
    for query in master_list:
        # Split into partition files to save memory space
        if count%1000 == 0:
            PubChem_out = 'C:/Users/bryan/Documents/Berkeley/4thYear/BioE134/BioE134_Proj5/E.coli_ID_Pubchem_' + str(count//1000) + '.txt'    # CHANGE_ME
            out_PubChem = open(PubChem_out, "w")
            out_PubChem.write(text)
            out_PubChem.close()
            text = ""

        handle = Entrez.esearch(db="pccompound", term=query)
        record = Entrez.read(handle)

        if len(list(record.keys())) > 0:
            count += 1
            key = list(record.keys())[3]
            text += query + "\t"
            text += ",".join(record[key]) + "\n"
            handle.close()

    remainder_out = 'C:/Users/bryan/Documents/Berkeley/4thYear/BioE134/BioE134_Proj5/E.coli_ID_Pubchem_remainder.txt'  # CHANGE_ME
    out_remainder = open(remainder_out, "w")
    out_remainder.write(text)
    out_remainder.close()

#NOTE= Each Step must be done one at a time, check files after each step to ensure accuracy

#Step 1
# master = initiate()
# obtain_uids(master)

#Step 2
# ID_lib = initiate2()
# obtain_SMILES(ID_lib)

#Step 3
combine_master_doc()
