# Create your tasks here
from __future__ import absolute_import, unicode_literals
from celery import shared_task
from celery.task import task
from django.core.mail import EmailMessage
from prettytable import PrettyTable
from seqdhbm import fasta
from seqdhbm import SeqDHBM
from SeqDHBM import models
import time

# About Celery
# # http://docs.celeryproject.org/en/latest/django/first-steps-with-django.html
# # http://docs.celeryproject.org/en/latest/getting-started/first-steps-with-celery.html#first-steps
# hemesitefolder$ celery -A hemesite worker -l debug --concurrency=15
# hemesitefolder$ celery -A hemesite worker -l debug --concurrency=1 -Q wesa,celery
# hemesitefolder$ celery -A hemesite beat -l info --scheduler django_celery_beat.schedulers:DatabaseScheduler
# https://www.youtube.com/watch?v=VoTxTM6kBuU&list=RDEM7S2Y1zXF-I77fRXilXb6Ew&index=27

@shared_task
def assync_organize_seq(seq_list):
    """
    GIVEN: The website is waiting for submissions from the user
    WHEN: a new job is created and the user submits multiple sequences inside a fasta file.
    THEN: The fasta sequences are saved in files in the jobs folder for future use in yasara.
    """

    time.sleep(1)
    fasta.organize_sequences2(seq_list)
    return "Files saved"


@shared_task
def assync_save_results(seq_idx, results_dict):
    """
    GIVEN: The user sent sequences for analysis in any of the provided input methods.
    WHEN: The sequence was analysed using the HBM Sequence detector.
    THEN: Save the results in the db.

    :param seq_idx:
    :param results_dict:
    :return: A message
    """

    time.sleep(.5)
    seq_obj = models.Sequence.objects.get(pk=seq_idx)
    for coord, res in results_dict.items():
        res_obj = models.Result_HBM(
            sequence=seq_obj,
            coord_atom=coord,
            ninemer=res["ninemer"],
            net_charge=res['netcharge'],
            disulfide_possible=bool(res['comment'].strip())
        )
        res_obj.save()
    if (seq_obj.status_hbm != models.Sequence.STATUS_FAILED and
       seq_obj.mode == models.Sequence.STRUCTURE_MODE):
        seq_obj.status_hbm = models.Sequence.STATUS_PROCESSED
    seq_obj.save()
    return "Results saved"


@shared_task
def access_wesa(seq_idx):
    """
    GIVEN: There are sequences in the database with pending WESA analysis.
    WHEN: Periodically, scheduled.
    THEN: Check if the results are ready in the WESA repository and update the records. """

    # check if the result is ready
    wesa_result = SeqDHBM.GetResultsFromWESA(seq_idx)
    if wesa_result:
        for result in models.Result_HBM.objects.filter(sequence=seq_idx):
            # check if the coordinating aa is not exposed
            if not wesa_result[int(result.coord_atom[1:])]:
                result.delete()
        seq_obj = models.Sequence.objects.get(pk=seq_idx)

        results = models.Result_HBM.objects.filter(sequence=seq_idx)
        if results:  # Check if there is any ninemer left
            # results.sort(key=lambda x: x.coord_atom[1:])
            table = PrettyTable(["S.no", "Coord. residue", "9mer motif", "Net charge", "Comment", "Kd or strength"])
            for pos, record in enumerate(sorted(results, key=lambda x: x.coord_atom[1:])):
                table.add_row([pos+1,
                               record.coord_atom,
                               record.ninemer,
                               record.net_charge,
                               "Possible S-S Bond" if record.disulfide_possible else " "*17,
                               ""])
            seq_obj.partial_hbm_analysis += "\n" + str(table)
        else:
            # and add the warning if they are all buried
            seq_obj.partial_hbm_analysis += f"\n{'*' * 80}\n"
            seq_obj.partial_hbm_analysis += "\nNOTE : THE SEQUENCE HAS NO SOLVENT ACCESSIBLE COORDINATION RESIDUES !\n"
            seq_obj.partial_hbm_analysis += "*" * 80
            seq_obj.warnings_hbm += "The sequence has no solvent accessible coordination residues"
        seq_obj.status_hbm = models.Sequence.STATUS_PROCESSED
        seq_obj.save()

        # update the partial full report
        seq_obj.jobnum.set_full_hbm_analysis()

        seqs_from_job = models.Sequence.objects.filter(jobnum=seq_obj.jobnum, status_hbm=models.Sequence.STATUS_QUEUED)
        if not seqs_from_job:
            if seq_obj.jobnum.submittedby:
                body = f"You can access the analysis at http://localhost:8000/SeqDHBM/{seq_obj.jobnum.id}"
                e_msg = EmailMessage(
                    subject=f'SeqD-HBM: Your analysis number {seq_obj.jobnum.id} is complete',
                    body=body,
                    from_email='seqdhbm@gmail.com',
                    to=[seq_obj.jobnum.submittedby],
                    headers={'Message-ID': 'foo'},
                )
                e_msg.send(fail_silently=False)


# Scheduled in the file celery.py
@task(name='check_wesa_server')
def check_for_pending_sequences():
    """
    GIVEN: The website is running.
    WHEN: Periodically, scheduled.
    THEN: Reads if there are any pending sequence. Try to read from the
    WESA server and update the records in the database."""

    queryset = models.Sequence.objects.filter(
            status_hbm=models.Sequence.STATUS_QUEUED,
            mode=models.Sequence.WESA_MODE)
    coco = "debugmessage:\n"
    for seq in queryset:
        coco += f"tasks.check_for_pending_sequences {seq.id}\n {seq.seqchain}\n"
        access_wesa.delay(seq.id)
    # DONE manage the status of the seq and jobs (wesa?)
    # DONE Load the sequence and all their results from db
    # DONE Check which of them are buried and delete
    # DONE If no coord aa is left, add the warning
    # TODO If no coord aa is left, the full analysis is not being writen
    # DONE update the partial analysis
    # DONE Check if there are still pending sequences for WESA_tmp
    # DONE If all were processed, update the job to finished, compile the full analysis
    # TEST and send an email to the submitter
    return coco
    # uncomment after debugging
    # return f"There were {len(queryset)} sequences pending"
