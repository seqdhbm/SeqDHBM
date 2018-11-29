
# coding: utf-8

# # Module to work with fasta files
# 
# ## Working with a fasta file containing multiple sequences

# In[18]:


import sys, os

def break_fasta_sequence(seq):
    assert type(seq)==str
    return [seq[x:x+70] for x in range(0, len(seq), 70)]  

def fasta_to_seq2(in_file, jobfolder="J0"):
    fullpath = os.path.join(jobfolder, in_file)
    if not os.path.isfile(fullpath):
        print("file not found.")
        return []
    with open(fullpath) as f:
        # reads the file in a list
        lines = f.read().splitlines()
    seq = ''     # stores the sequence without linebreaks
    seq_dict={}  # stores the dictionary to be returned
    result = []
    head= ""
    seq_id = 1
    for i in lines:
        if i and i[0] == '>':
            # Header found
            if seq:   # if seq is empty, this might be the first line of the file
                # get the first section of the header, hopefully the ID
                result+= [{"name": head[1:], 
                           "seq":seq, 
                           "folder": os.path.join(jobfolder, "FF%05d"%seq_id),
                           "file": "%d.fasta"%(seq_id),
                           "submited_as": "Fasta file"}]
                seq_id+=1
                seq=''
            head=i
        else: # not a header - append to sequnce
            seq+=i
    if seq and head: 
        # after finishing reading the file, check if there is a sequence in the 
        # variable and add it to the dictionary accordingly.
        result+= [{"name": head[1:], 
                   "seq":seq, 
                   "folder": os.path.join(jobfolder, "FF%05d"%seq_id),
                   "file": "%d.fasta"%(seq_id),
                   "submited_as": "Fasta file"}]
    return result



def organize_sequences2(seq_list):
    """Creates separated fasta files organized in folders from a single fasta file with multiple sequences.
    
    Keyword arguments:
    in_file -- Original fasta file."""
    for k in seq_list:
        os.makedirs(k["folder"], exist_ok=True)
        file = os.path.join(k["folder"], k["file"])
        with open(file, 'w') as f:
            f.write(">%s\n"%k["name"])
            lines=break_fasta_sequence(k["seq"])
            f.write('\n'.join(lines))
            f.write('\n')
            


# In[23]:


if (__name__ == "__main__"):
    print("Unit testing")
    a = fasta_to_seq2("mock1.fasta", jobfolder="/home/mau/work")
    organize_sequences2(a)


# In[7]:


'''def fasta_to_seq(in_file):
    """ Parses a fasta file with multiple sequences and creates a dictionary out of the records.
    
    Returns a dictionary of sequences: The dictionary key
    will be the folder name and the values
    are tuples containing the full header and the sequence"""
    if not (in_file):
        return {}
    with open(in_file) as f:
        # reads the file in a list
        lines = f.read().splitlines()
    seq = ''     # stores the sequence without linebreaks
    seq_dict={}  # stores the dictionary to be returned
    for i in lines:
        if i and i[0] == '>':
            # Header found
            if seq:   # if seq is empty, this might be the first line of the file
                # get the first section of the header, hopefully the ID
                outfile=head[1:].split('|')[1]+'('+ str(len(seq))+')'
                cnt=1
                # check if there are sequences with same ID/length
                while outfile in seq_dict:
                    outfile=i[1:].split('|')[1]+'('+ str(len(seq))+')-'+str(cnt)
                    cnt+=1
                seq_dict[outfile]=(head, seq)
                seq=''
            head=i
        else: # not a header - append to sequnce
            seq+=i
    if seq: 
        # after finishing reading the file, check if there is a sequence in the 
        # variable and add it to the dictionary accordingly.
        outfile=head[1:].split('|')[1]+'('+ str(len(seq))+')'
        cnt=1
        while outfile in seq_dict:
            outfile=i[1:].split('|')[1]+'('+ str(len(seq))+')-'+str(cnt)
            cnt+=1
        seq_dict[outfile]=(head, seq)
    return seq_dict''';
   
'''
def organize_sequences(seq_dict):
    """Creates separated fasta files organized in folders from a single fasta file with multiple sequences.
    
    Keyword arguments:
    in_file -- Original fasta file."""
    for k, (h, s) in seq_dict.items():
        os.makedirs(k, exist_ok=True)
        with open(k+'/'+k+'.fasta', 'w') as f:
            f.write(h+'\n')
            lines=break_fasta_sequence(s)
            f.write('\n'.join(lines))
            ''';


# In[10]:




