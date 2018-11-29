#!/usr/bin/env python
# coding: utf-8

# In[2]:


import Bio.PDB as bpdb
from os.path import isfile
from os import rmdir
from Bio.PDB.Polypeptide import three_to_one
import logging
logging.basicConfig(level=logging.DEBUG, format='%(levelname)s - %(message)s')


# In[1]:


# TODO: Create a list of warnings/errors for the user

def get_pdb_files(id_list, job = "J0"):
    """ Returns a dictionary containing {result[id]: 'filename'}.
    """
    file_list = dict()
    pdbl = bpdb.PDBList()
    for i, pid in enumerate(id_list, 1):
        filename = pdbl.retrieve_pdb_file(pid, pdir="P%05d"%i, file_format="pdb")
        if isfile(filename):
            file_list[pid] = filename
        else:
            logging.warning("No structure found in PDB for the ID '"+ pid+ "'")
            try:
                rmdir(pid)
            except:
                logging.error("Couldn't delete directory for "+ str(pid))
    return file_list

def get_pdb_files2(id_list, job = "J0"):
    """ Returns a dictionary containing {result[id]: 'filename'}.
    """
    result = []
    pdbl = bpdb.PDBList()
    for i, pid in enumerate(id_list, 1):
        filename = pdbl.retrieve_pdb_file(pid, pdir="P%05d"%i, file_format="pdb")
        if isfile(filename):
            result += [{"name": pid, 
                        "folder": "P%05d"%i}]
        else:
            logging.warning("No structure found in PDB for the ID '"+ pid+ "'")
            try:
                rmdir(pid)
            except:
                logging.error("Couldn't delete directory for "+ str(pid))
    return file_list

def text_to_list(userinput):
    """Takes the input from the user as a free text containing protein IDs separated
    by commas or semi-colons.
        
    Keyword arguments:
    userinput -- Text provided by the user.
    
    Returns:
    A list of protein IDs"""
    result = "".join(userinput.split(" ")) # remove spaces
    result = ",".join(result.split(";")) # converts ';' in ','
    # if there are other symbols to consider, follow the pattern of the commands above
    return result.split(",") 

if __name__ == "__main__" and False:
    t = text_to_list("1S0L")
    get_pdb_files(t)


# In[3]:


# I COULDNT FIGURE OUT HOW TO WORK WITH mmCIF FILES.
# MORE RESEARCH NEEDS TO BE DONE IN ORDER TO OBTAIN THE PRIMARY SEQ FROM THOSE FILES.

def cif_to_fasta(filename):
    """Parses a cif(pdb) file to recover the primary sequence.
    
    Expects to find the tag '_entity_poly.pdbx_seq_one_letter_code'\
    in the file followed by the sequence (delimited by ';') to mark this field."""
    with open(filename) as f: 
        lines = "".join(f.read().splitlines())
    pos = lines.find("_entity_poly.pdbx_seq_one_letter_code")
    if pos==-1:
        logging.warning("File "+ filename+ " does not contain information about the primary sequence. \nSkipping.") # TODO: error list
        return ""
    else:
        seq = lines[pos+45:]
        return seq[:seq.find(";")]


# In[8]:



from Bio.PDB.Polypeptide import three_to_one

def pdb_to_fasta(filename):
    raise Exception("addapt this function to work with the new format of the dictionary")
    """Parses a pdb file to recover the primary sequence.
    
    Expects to find the tag 'SEQRES' in the file to mark this field."""
    with open(filename) as f: 
        lines = f.readlines()
    lines = [x[19:] for x in lines if x[:6]=="SEQRES"]
    #print(lines)
    if lines:
        lines = [" ".join(x.split()) for x in lines] # remove extra whitespaces
        seq = [x.split() for x in lines] # divide in 3 letter strings
        result = ""
        try:
            for l in seq:
                result += "".join([three_to_one(x) for x in l])
        except KeyError as e:
            logging.warning("File "+ filename+ " contains nonstandard amino acids. \nSkipping.") # TODO: error list
            return ""
        return (result)
    else:
        logging.warning("File "+ filename+ " does not contain information about the primary sequence. \nSkipping.") # TODO: error list

def structure_to_fasta(filename):
    if filename[-3:].upper() in ["PDB", "ENT"]:
        seq = pdb_to_fasta(filename)
    elif filename[-3:].upper() == "CIF":
        logging.warning("Warning: Not ready to handle this type of files on SeqD-HBM") 
        logging.warning(filename)
        # TODO: MORE RESEARCH NEEDS TO BE DONE IN ORDER TO OBTAIN THE PRIMARY SEQ FROM THOSE FILES.
        #seq = cif_to_fasta(filename)
        seq = ""
    else:
        raise Exception("Error: unexpected file type")
    """ # uncomment if we ever need the fasta file for something in the future.
    if seq:
        fastafile = filename[:-3]+"fasta"
        with open(fastafile, "w") as f:
            f.write(">conv|"+fastafile+"\n") #TODO: proper header
            f.write(seq+"\n")"""
    return seq


# In[9]:


if __name__ == "__main__":
    test = ["1S0L/1s0l.cif", "1I1B/1i1b.cif", "1HIB/1hib.cif", "1I8H/pdb1i8h.ent"]
    for i in test:
        print("#"*50)
        structure_to_fasta(i)


# In[ ]:




