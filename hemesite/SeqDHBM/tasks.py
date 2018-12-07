# Create your tasks here
from __future__ import absolute_import, unicode_literals
from celery import shared_task
from celery.task import task
import re
from seqdhbm import fasta
from seqdhbm import SeqDHBM
from SeqDHBM import models
import os
import time
# About Celery
# # http://docs.celeryproject.org/en/latest/django/first-steps-with-django.html
# # http://docs.celeryproject.org/en/latest/getting-started/first-steps-with-celery.html#first-steps
# hemesitefolder$ celery -A hemesite worker -l debug --concurrency=15
# hemesitefolder$ celery -A hemesite worker -l debug --concurrency=1 -Q wesa,celery
# hemesitefolder$ celery -A hemesite beat -l info --scheduler django_celery_beat.schedulers:DatabaseScheduler

@shared_task
def assync_organize_seq(seq_list):
    """Save the files in the jobs folder for future use in yasara"""

    time.sleep(1)
    fasta.organize_sequences2(seq_list)
    return "Files saved"


@shared_task
def assync_save_results(seq_idx, results_dict):
    """Save the results in the db."""

    time.sleep(.5)
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
    if seq_obj.status_hbm != models.Sequence.STATUS_FAILED:
        seq_obj.status_hbm = models.Sequence.STATUS_PROCESSED
    seq_obj.save()
    return "Results saved"


@shared_task
def access_wesa(seq_idx):
    """ Get the results from the  WESA repository """

    wesa_result = SeqDHBM.GetResultsFromWESA(seq_idx)
    if wesa_result:
        print("tasks.py line 59ish: ", wesa_result)
        for result in models.Result_HBM.filter(sequence=seq):
            # check if the coordinating aa is not exposed
            if not wesa_result[int(result.coord_atom[1:])]:
                result.delete()
        seq.status_hbm = models.Sequence.STATUS_PROCESSED

# Scheduled in the file celery.py
@task(name='check_wesa_server')
def check_for_pending_sequences():
    """ Reads if there are any pending sequence. Try to read from the
    WESA server and update the records in the database."""

    queryset = models.Sequence.objects.filter(
            status_hbm=models.Sequence.STATUS_QUEUED,
            mode=models.Sequence.WESA_MODE)
    for seq in queryset:
        access_wesa.delay(seq.id)

    # TODO manage the status of the seq and jobs (wesa?)
    # DONE Load the sequence and all their results from db
    # DONE Check which of them are buried and delete
    # TODO If no coord aa is left, add the warning
    # TODO update the partial analysis
    # TODO Check if there are still pending sequences for WESA_tmp
    # TODO If all were processed, update the job to finished
    # TODO and compile the full analysis
    # TODO and send an email to the submitter
    # TODO only one instance of this function can be running at a time
    # TODO Continue from here
    return f"There were {len(queryset)} sequences pending"
