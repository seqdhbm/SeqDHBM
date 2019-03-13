# coding: utf-8

from django.conf import settings
# from django.contrib.sites.models import Site
from django.core.mail import EmailMessage
import os
from SeqDHBM import models
import seqdhbm.workflow as wf

__all__ = [
    'handle_uploaded_file', 'clean_raw_seq'
]

APP_URL_PATTERN = ""


def handle_uploaded_file(f, filename):
    """
    WHEN: The user submits a form using a POST method.
    GIVEN: The FILES section of the form is not empty.
    THEN: download the files, so they can be processed.

    :param f: the binary information sent in the form.
    :param filename: the file path to save the file
    """

    with open(filename, 'wb+') as destination:
        for chunk in f.chunks():
            destination.write(chunk)


def clean_raw_seq(input_seq: str) -> str:
    """Removes spaces and linebreaks from a string"""

    _clean = "".join(input_seq.split("\n"))
    _clean = "".join(_clean.split())
    _clean = "".join(_clean.split("\r"))
    return _clean


def process_form(form, submitted_files):
    """

    1. Validate the form data,
    2. Save the information sent by the user,
    3. Analyse the data,
    4. Update the database with the result
    5. Inform the user

    :param form: The http form containing the information.
    :param submitted_files: the files submitted in the request.
    :return: The number of the job and the password to access it.
    """
    ##############################################
    # Prepare and validate the data from the form
    ##############################################
    email = form.cleaned_data['email'] if form.cleaned_data['email'] else ""
    # Create the job in the database
    myjob = models.Job(submittedby=email)
    myjob.save()
    jobfolder = os.path.join(settings.BASE_DIR, "jobs/J%07d/" % myjob.id)
    os.makedirs(jobfolder, exist_ok=True)
    # save the fastafile in disk
    if "fastafiles" in submitted_files:
        fastafile = submitted_files['fastafiles'].name
        fastafilepath = os.path.join(
            jobfolder,
            submitted_files['fastafiles'].name
        )
        handle_uploaded_file(submitted_files['fastafiles'], fastafilepath)
    else:
        fastafile = ""
    # Obtain the manual input and clean it
    rawseq = clean_raw_seq(form.cleaned_data['rawseq'])

    ###########################################################
    # run the workflow seqdhbm analysis for the current inputs
    ###########################################################
    result = wf.workflow(jobfolder=jobfolder,
                         rawseq=rawseq,
                         fastafile=fastafile,
                         job=myjob,
                         mode=form.cleaned_data["mode"]
                         )

    #####################################################
    # Update the database and inform user of his request
    #####################################################
    myjob.set_full_hbm_analysis([x["analysis"] for x in result
                                 if "analysis" in x])
    # if all went well, send the email to the user
    # https://docs.djangoproject.com/en/2.1/topics/email/#django.core.mail.EmailMessage
    email = form.cleaned_data['email']
    password = myjob.pass_gen()
    if email:
        site_domain = settings.SITE_DOMAIN
        body = "Dear user,\n\nYour results for submission " +\
                f"{myjob.id} is available at " + \
               f"http://{site_domain}/SeqDHBM/{myjob.id}/{password}" +\
               "\n\n Thank you for using SeqD-HBM!"
        e_msg = EmailMessage(
            subject=f'SeqD-HBM: your submission {myjob.id}',
            body=body,
            from_email='seqdhbm@gmail.com',
            to=[email],
            headers={'Message-ID': 'foo'},
        )
        e_msg.send(fail_silently=False)
    return myjob.id, password
