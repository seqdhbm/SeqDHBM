# coding: utf-8

"""Functions to handle user submitted fasta files or sequences"""

import os


def break_fasta_sequence(seq):
    assert type(seq) == str
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
    result = []
    head = ""
    seq_id = 1
    for i in lines:
        if i and i[0] == '>':
            # Header found
            if seq:  # if seq is empty, this might be the 1st line of the file
                # get the first section of the header, hopefully the ID
                result += [{"name": head[1:],
                            "seq":seq,
                            "folder": os.path.join(jobfolder, "FF%05d" % seq_id),
                            "file": "%d.fasta" % seq_id,
                            "submitted_as": "Fasta File"}]
                seq_id += 1
                seq = ''
            head = i
        else:  # not a header - append to sequnce
            seq += i
    if seq and head:
        # after finishing reading the file, check if there is a sequence in the
        # variable and add it to the dictionary accordingly.
        result += [{"name": head[1:],
                    "seq":seq,
                    "folder": os.path.join(jobfolder, "FF%05d" % seq_id),
                    "file": "%d.fasta" % seq_id,
                    "submitted_as": "Fasta File"}]
    return result


def organize_sequences2(seq_list):
    """Creates separated fasta files organized in folders from a single fasta file with multiple sequences.

    Keyword arguments:
    in_file -- Original fasta file."""
    print("debugging")
    for k in seq_list:
        os.makedirs(k["folder"], exist_ok=True)
        file = os.path.join(k["folder"], k["file"])
        with open(file, 'w') as f:
            f.write(">%s\n" % k["name"])
            lines = break_fasta_sequence(k["seq"])
            f.write('\n'.join(lines))
            f.write('\n')
