from __future__ import absolute_import, unicode_literals
from django.conf import settings
from django.core.mail import EmailMessage
from django.http import HttpResponse, Http404
from django.shortcuts import get_object_or_404, get_list_or_404, redirect
from django.template import loader
from .forms import SeqSubmission
import seqdhbm.workflow as wf
from SeqDHBM import models
import os

# DOING tasks.check_for_pending_sequences
# TODO Ask Ajay to test the wesa
# TODO Model a db for the group's sequences
# TODO write if the sequence was predicted by wesa or not.. show results and red message
# TODO check for errors in the html
# TODO check if the sequence has the minimum size for wesa
# TODO Upgrade the software to use yasara
# TODO Read papers from the Bioinformatics Journals
# TODO setup a different sgbd/web server (GUnicorn / apache)
# TODO Minor adjustments (require the email to be from university, security measures -> requires apache)
# DONE: use the timestamp to protect the results
# DONE NEXT: Database
#            -  WESA queue
#                -  One function to send to wesa
#                -  One function (scheduled) to check if there are pending jobs
#                -  Update fields in tables job and Sequence to control pending
# DONE:
#            -  Save to database via celery ok
#            -  Create jobs ok
#            -  save the submission ok
#            -  save the results ok
#            -  make the results accessible ok
#            -  validate if the email is filled when wesa is being used
# Future: celery has a 'chain' function that might be handy to pipe the
#    homology -> docking -> md process
# Future: start celery as a daemon:
#  http://docs.celeryproject.org/en/latest/userguide/daemonizing.html#daemonizing
# Future: read about security: https://docs.djangoproject.com/en/2.1/topics/security/#user-uploaded-content-security
# Future: check a service to verify if the email belongs to an academic institution


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


def index(request):
    """
    GIVEN: The website is active.
    WHEN: The user accesses the main page or submits a job to be processed.
    THEN: Show the webpage or validate the information and process the job.

    :param request: the http request (POST or GET)
    :return: The home webpage of the site (or a redirect to the results page).
    """
    result = []
    debugging = ""
    # if this is a POST request we need to process the form data
    if request.method == 'POST':
        # create a form instance and populate it with data from the request:
        form = SeqSubmission(request.POST, request.FILES)

        debugging = form.non_field_errors()
        # check whether it's valid:
        if form.is_valid():
            email = form.cleaned_data['email'] if form.cleaned_data['email'] else ""
            # Create the job in the database
            myjob = models.Job(submittedby=email)
            myjob.save()
            jobfolder = os.path.join(settings.BASE_DIR, "jobs/J%07d/" % myjob.id)
            os.makedirs(jobfolder, exist_ok=True)
            # obtain the fastafile
            if "fastafiles" in request.FILES:
                fastafile = request.FILES['fastafiles'].name
                fastafilepath = os.path.join(jobfolder, request.FILES['fastafiles'].name)
                handle_uploaded_file(request.FILES['fastafiles'], fastafilepath)
            else:
                fastafile = ""
            # Obtain the manual input and clean it
            rawseq = "".join(form.cleaned_data['rawseq'].split("\n"))
            rawseq = "".join(rawseq.split())
            rawseq = "".join(rawseq.split("\r"))

            # run the workflow seqdhbm analysis for the current inputs
            result = wf.workflow(jobfolder=jobfolder,
                                 rawseq=rawseq,
                                 fastafile=fastafile,
                                 job=myjob,
                                 mode=form.cleaned_data["mode"]
                                 )

            # Update the full analysis of the job using the
            myjob.set_full_hbm_analysis([x["analysis"] for x in result])
            # if all went well, send the email to the user
            # https://docs.djangoproject.com/en/2.1/topics/email/#django.core.mail.EmailMessage
            email = form.cleaned_data['email']
            password = myjob.pass_gen()
            if email:
                body = "You can access the analysis at http://localhost:8000/SeqDHBM/%d/%s" % (myjob.id, password)
                e_msg = EmailMessage(
                    subject=f'SeqD-HBM: Your analysis number {myjob.id}',
                    body=body,
                    from_email='seqdhbm@gmail.com',
                    to=[email],
                    headers={'Message-ID': 'foo'},
                )
                e_msg.send(fail_silently=False)
            return redirect("%d/%s" % (myjob.id, password))
    # if a GET (or any other method) we'll create a blank form
    else:
        form = SeqSubmission()

    template = loader.get_template('SeqDHBM/index.html')
    context = {
        'results': result,
        'form': form,
        'debug': debugging
    }
    return HttpResponse(template.render(context, request))


def show_result(request, job_id, passw: str = ""):
    """
    GIVEN: The database of jobs (completed or queued)
    WHEN: The user requests information about a job
    THEN: List the number of pending jobs as well as html information for the processed jobs.

    :param request: the http request.
    :param job_id: The identifier of the job.
    :param passw: Security code to the given job
    :return: The formatted webpage with the requested information. The template comes from one http file.
    """
    template = loader.get_template('SeqDHBM/result.html')
    job = get_object_or_404(models.Job, id=job_id)
    seqs = get_list_or_404(models.Sequence, jobnum=job)
    if job.pass_gen() == passw:
        res = {}
        for seq in seqs:
            res[seq] = [seq.seqchain[x:x+70] for x in range(0, len(seq.seqchain), 70)]
        context = {
            "job": job,
            "processing": len([x for x in seqs if x.status_hbm == models.Sequence.STATUS_QUEUED]),
            "result": {x: y for x, y in res.items() if x.status_hbm == models.Sequence.STATUS_PROCESSED},
            "failed": {x: y for x, y in res.items() if x.status_hbm == models.Sequence.STATUS_FAILED},
            "passw": passw
        }
        return HttpResponse(template.render(context, request))
    else:
        raise Http404("Invalid URL (bad secret code and job id combination)")


def show_analysis(request, job_id, passw: str = ""):
    """
    WHEN: The user requests information about a job he submitted.
    GIVEN: The job exists.
    THEN: The full analysis of the job is presented to the user as a text webpage.

    :param request: http request
    :param job_id: The identifier of the job
    :param passw: Security code to the given job
    :return: The formatted webpage with the requested information.
    """
    job = get_object_or_404(models.Job, id=job_id)
    if job.pass_gen() == passw:
        return HttpResponse(job.full_hbm_analysis, content_type="text/plain")
    else:
        raise Http404("Invalid URL (bad secret code and job id combination)")


def hemewf(request, job_id, passw):
    """
    for debugging.
    /runworkflow/job_id/passw

    :param request:
    :return:
    """
    message = ""
    for seq in models.Sequence.objects.all():
        message += "<div>"
        message += "<p>Job num: %d</p>" % seq.jobnum.id
        message += "<p>Sequence: %s</p>" % seq.seqchain
        message += "<p>Sent by: %s</p>" % seq.submittedas
        message += "<p>Prediction: %s</p>" % seq.mode
        message += "<p>Status of motif detection: %s</p>" % seq.status_hbm
        message += "<p>file location: %s</p>" % seq.fasta_file_location
        message += "<p>Warnings: %s</p>" % ("<br>".join(seq.warnings_hbm.split("\n")))
    message += "<h1>Results</h1>"
    for res in models.Result_HBM.objects.all():
        message += "<p>Sequence num: %d</p>" % res.sequence.id
        message += "<p>Coord atom: %s</p>" % res.coord_atom
    # return HttpResponse(message)
    message += "\n\n\n"
    for seq in models.Sequence.objects.filter(jobnum=models.Job.objects.get(pk=27)):
        message += seq.partial_hbm_analysis + "\n"
    job = models.Job.objects.get(pk=27)
    if job.pass_gen() == passw:
        print(job.id)
        print()
        message += f"\n{request.META}"
    return HttpResponse(message, content_type="text/plain")
