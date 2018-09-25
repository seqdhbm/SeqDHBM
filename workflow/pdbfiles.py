
# coding: utf-8

# In[1]:


import Bio.PDB as bpdb
from os.path import isfile
from os import rmdir
from Bio.PDB.Polypeptide import three_to_one


# In[2]:


# TODO: Create a list of warnings/errors for the user

def get_pdb_files(id_list):
    """ Returns a dictionary containing {result[id]: 'filename'}.
    """
    file_list = dict()
    for pid in id_list:
        pdbl = bpdb.PDBList()
        filename = pdbl.retrieve_pdb_file(pid, pdir=pid, file_format="pdb")
        if isfile(filename):
            file_list[pid] = filename
        else:
            print("Warning: No structure found in PDB for the ID '"+ pid+ "'")
            try:
                rmdir(pid)
            except:
                print("Couldn't delete directory for", pid)
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

if False:
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
        print("Warning: File", filename, "does not contain information about the primary sequence. Skipping.") # TODO: error list
        return ""
    else:
        seq = lines[pos+45:]
        print("dbg",seq[:seq.find(";")])
        return seq[:seq.find(";")]


# In[8]:



from Bio.PDB.Polypeptide import three_to_one

def pdb_to_fasta(filename):
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
            print("Warning: File", filename, "contains nonstandard amino acids. Skipping.") # TODO: error list
            return ""
        return (result)
    else:
        print("Warning: File", filename, "does not contain information about the primary sequence. Skipping.") # TODO: error list

def structure_to_fasta(filename):
    if filename[-3:].upper() in ["PDB", "ENT"]:
        seq = pdb_to_fasta(filename)
    elif filename[-3:].upper() == "CIF":
        print("Warning: Not ready to handle this type of files on SeqD-HBM") 
        print(filename)
        # TODO: MORE RESEARCH NEEDS TO BE DONE IN ORDER TO OBTAIN THE PRIMARY SEQ FROM THOSE FILES.
        #seq = cif_to_fasta(filename)
        seq = ""
    else:
        raise Exception("Error: unexpected file type")
    if seq:
        fastafile = filename[:-3]+"fasta"
        with open(fastafile, "w") as f:
            f.write(">conv|"+fastafile+"\n") #TODO: proper header
            f.write(seq+"\n")
        return seq
"""
test = ["1S0L/1s0l.cif", "1I1B/1i1b.cif", "1HIB/1hib.cif"]
for i in test:
    print(structure_to_fasta("/home/imhof_team/Public/mauricio/workflow/"+i))
# _entity_poly.pdbx_seq_one_letter_code
# _pdbx_poly_seq_scheme.hetero""";


# In[9]:


if __name__ == "__main__":
    test = ["1S0L/1s0l.cif", "1I1B/1i1b.cif", "1HIB/1hib.cif", "1I8H/pdb1i8h.ent"]
    for i in test:
        print("#"*50)
        structure_to_fasta(i)


# In[1]:


"a".split()

