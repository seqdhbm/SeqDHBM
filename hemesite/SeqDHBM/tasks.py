# Create your tasks here
from __future__ import absolute_import, unicode_literals
from celery import shared_task
from seqdhbm import fasta
import os

@shared_task
def assync_organize_seq(seq_list):
    for k in seq_list:
        os.makedirs(k["folder"], exist_ok=True)
        file = os.path.join(k["folder"], k["file"])
        with open(file, 'w') as f:
            f.write(">%s\n"%k["name"])
            lines=fasta.break_fasta_sequence(k["seq"])
            f.write('\n'.join(lines))
            f.write('\n')
    return "Files saved"
