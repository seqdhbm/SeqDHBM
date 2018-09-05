
# coding: utf-8

# In[58]:


import sys, os
in_file='IL-1_formatted.fasta'


def fasta_to_seq(in_file):
    """ Reads a fasta file with multiple sequences and creates a dictionary out of the records.
The folder name will be the key and the values are tuples containing the header and the sequence"""
    with open(in_file) as f:
        lines = f.read().splitlines()
    seq = ''
    seq_dict={}
    for i in lines:
        if i[0] == '>':
            if seq:
                outfile=head[1:].split('|')[0]+'('+ str(len(seq))+')'
                cnt=1
                while outfile in seq_dict:
                    outfile=i[1:].split('|')[0]+'('+ str(len(seq))+')-'+str(cnt)
                    cnt+=1
                seq_dict[outfile]=(head, seq)
                seq=''
            head=i
        else:
            seq+=i
    if seq:
        outfile=head[1:].split('|')[0]+'('+ str(len(seq))+')'
        cnt=1
        while outfile in seq_dict:
            outfile=i[1:].split('|')[0]+'('+ str(len(seq))+')-'+str(cnt)
            cnt+=1
        seq_dict[outfile]=(head, seq)
    return seq_dict

def organize_sequences(in_file):
    """ Create python doc """
    seqs = fasta_to_seq(in_file)
    for k, (h, s) in seqs.items():
        os.makedirs(k, exist_ok=True)
        with open(k+'/'+k+'.fasta', 'w') as f:
            f.write(h+'\n')
            lines=[s[x:x+60] for x in range(0, len(s), 60)]
            f.write('\n'.join(lines))
            
if __name__ == '__main__':
    assert len(sys.argv) == 2, 'Usage is "python organize_seqs <fastafile>"'
    organize_sequences(sys.argv[1])


# In[49]:




