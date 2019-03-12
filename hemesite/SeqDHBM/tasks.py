# coding: utf-8

"""Asynchronous functions."""

from __future__ import absolute_import, unicode_literals
import time

from celery import shared_task
from celery.task import task
from django.conf import settings
from django.core.mail import EmailMessage
from prettytable import PrettyTable
from seqdhbm import fasta
from seqdhbm import SeqDHBM

from SeqDHBM import models

__all__ = [
    'assync_organize_seq', 'assync_save_results', 'access_wesa',
    'check_for_pending_sequences'
]

# About Celery
# http://docs.celeryproject.org/en/latest/django/first-steps-with-django.html
# http://docs.celeryproject.org/en/latest/getting-started/first-steps-with-celery.html#first-steps
# hemesitefolder$ celery -A hemesite worker --detach -l debug --concurrency=15
# hemesitefolder$ celery -A hemesite worker --detach -l debug --concurrency=1 -Q wesa,celery
# hemesitefolder$ celery -A hemesite beat --detach -l info --scheduler django_celery_beat.schedulers:DatabaseScheduler
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
    WHEN: Periodically, scheduled in celery beat.
    THEN: Check if the results are ready in the WESA repository and update the records.

    :param seq_idx: the id of the sequence record in the db
    """

    # check if the result is ready
    no_solvent_analysis_msg = [
        f"\n{'*' * 80}\n",
        "NOTE: THE SEQUENCE HAS NO SOLVENT ACCESSIBLE COORDINATION RESIDUES!",
        "\n",
        "*" * 80
    ]
    no_solvent_warn = "\nThe sequence has no solvent " +\
                      "accessible coordination residues"
    wesa_error_msg = [
        "\n",
        f"{'*' * 80}\n ",
        f"Unexpected WESA Error: {e}\n ",
        "Try again after a while"
    ]
    wesa_error_but_results = [
        "\n",
        f"{'*' * 80}\n",
        "Showing the results without solvent accessibility analysis",
        "\n"
    ]
    try:
        wesa_result = SeqDHBM.get_results_from_wesa(seq_idx)
        if wesa_result:
            for result in models.Result_HBM.objects.filter(sequence=seq_idx):
                # check if the coordinating aa is not exposed
                if not wesa_result[int(result.coord_atom[1:])]:
                    result.delete()
            seq_obj = models.Sequence.objects.get(pk=seq_idx)

            results = models.Result_HBM.objects.filter(sequence=seq_idx)
            if results:  # Check if there is any ninemer left
                seq_obj.partial_hbm_analysis += "\n" + result_to_pretty_table(results)
            else:
                # and add the warning if they are all buried
                for line in no_solvent_analysis_msg:
                    seq_obj.partial_hbm_analysis += line
                seq_obj.warnings_hbm += no_solvent_warn
            seq_obj.status_hbm = models.Sequence.STATUS_PROCESSED
            seq_obj.save()
        else:
            return
    except AssertionError as e:
        seq_obj = models.Sequence.objects.get(pk=seq_idx)
        for line in wesa_error_msg:
            seq_obj.partial_hbm_analysis += line
        results = models.Result_HBM.objects.filter(sequence=seq_idx)
        if results:  # Check if there is any ninemer left
            for line in wesa_error_but_results:
                seq_obj.partial_hbm_analysis += line
            seq_obj.partial_hbm_analysis += result_to_pretty_table(results)
        seq_obj.warnings_hbm += f"\n Unexpected WESA Error: {e}\n"
        seq_obj.warnings_hbm += "No solvent accessibility " +\
                                "prediction could be done."
        seq_obj.status_hbm = models.Sequence.STATUS_FAILED
        seq_obj.save()

    # update the partial full report
    seq_obj.jobnum.set_full_hbm_analysis()

    seqs_from_job = models.Sequence.objects.filter(
        jobnum=seq_obj.jobnum,
        status_hbm=models.Sequence.STATUS_QUEUED
    )
    # if no more sequences are pending for the job and the user filled
    #  an email address, send him the analysis
    if not seqs_from_job and seq_obj.jobnum.submittedby:
        site_domain = settings.SITE_DOMAIN
        email_subject = f'SeqD-HBM: Your analysis ' +\
                        f'number {seq_obj.jobnum.id} is complete'
        email_body = f"You can access your results at " + \
                     f"http://{site_domain}/SeqDHBM/" + \
                     f"{seq_obj.jobnum.id}/" + \
                     f"{seq_obj.jobnum.pass_gen()}"
        e_msg = EmailMessage(
            subject=email_subject,
            body=email_body,
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
    THEN: If there are any pending sequence, try to read from the WESA server and
    update the records in the database.

    :return: A message to the Queue Handler
    """

    queryset = models.Sequence.objects.filter(
            status_hbm=models.Sequence.STATUS_QUEUED,
            mode=models.Sequence.WESA_MODE)
    for seq in queryset:
        access_wesa.delay(seq.id)
    return f"There were {len(queryset)} sequences pending"


def result_to_pretty_table(results):
    """
    Formats the result as a pretty table

    :param results: Results object
    :return: A pretty table with the analysis
    """

    table = PrettyTable([
        "S.no",
        "Coord. residue",
        "9mer motif",
        "Net charge",
        "Comment",
        "Kd or strength"])  # TODO Predicted Kd? Strength as another field?
    for pos, record in enumerate(sorted(results, key=lambda x: int(x.coord_atom[1:]))):
        table.add_row([pos + 1,
                       record.coord_atom,
                       record.ninemer,
                       record.net_charge,
                       "Possible S-S Bond" if record.disulfide_possible else " " * 17,
                       ""])
    return str(table)
