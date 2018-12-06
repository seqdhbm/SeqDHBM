# Create your tasks here
from __future__ import absolute_import, unicode_literals
from celery import shared_task
from seqdhbm import fasta
from SeqDHBM import models
import os
import time
# About Celery
# # http://docs.celeryproject.org/en/latest/django/first-steps-with-django.html
# # http://docs.celeryproject.org/en/latest/getting-started/first-steps-with-celery.html#first-steps
# hemesitefolder$ celery -A hemesite worker -l debug --concurrency=15
# hemesitefolder$ celery -A hemesite worker -l debug --concurrency=1 -Q wesa,celery

@shared_task(queue="files")
def assync_organize_seq(seq_list):
    time.sleep(1)
    fasta.organize_sequences2(seq_list)
    return "Files saved"

@shared_task
def assync_save_results(seq_idx, results_dict):
    time.sleep(1)
    seq_obj = models.Sequence.objects.get(pk=seq_idx)
    for coord, res in results_dict.items():
        res_obj = models.Result_HBM(
            sequence = seq_obj,
            coord_atom = coord,
            ninemer = res["ninemer"],
            net_charge = res['netcharge'],
            disulfide_possible = bool(res['comment'].strip())
        )
        res_obj.save()
    return "Results saved"
