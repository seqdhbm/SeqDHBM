
# coding: utf-8

# # Module to work with fasta files
# 
# ## Working with a fasta file containing multiple sequences

# In[2]:


import sys, os

def fasta_to_seq(in_file):
    """ Parses a fasta file with multiple sequences and creates a dictionary out of the records.
    Returns a dictionary of sequences: The dictionary key will be the folder name and the values
    are tuples containing the full header and the sequence"""
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
    return seq_dict

def organize_sequences(seq_dict):
    """Creates separated fasta files organized in folders from a single fasta file with multiple sequences.
    
    Keyword arguments:
    in_file -- Original fasta file."""
    for k, (h, s) in seq_dict.items():
        os.makedirs(k, exist_ok=True)
        with open(k+'/'+k+'.fasta', 'w') as f:
            f.write(h+'\n')
            lines=[s[x:x+60] for x in range(0, len(s), 60)]
            f.write('\n'.join(lines))


# In[ ]:


if (__name__ == "__main__"):
    seq_dict = fasta_to_seq("IL-1_formatted.fasta")
    organize_sequences(seq_dict)
    for folder in seq_dict.keys():
        print(folder)
        yasara_homology(folder+'/'+folder+'.fasta')

