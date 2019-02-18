from __future__ import absolute_import, unicode_literals
from django.conf import settings
from django.core.mail import EmailMessage
from django.http import HttpResponse
from django.shortcuts import get_object_or_404, get_list_or_404, redirect
from django.template import loader
from .forms import SeqSubmission
import seqdhbm.workflow as wf
from SeqDHBM import models
import os

# DOING tasks.check_for_pending_sequences
# TODO validate if the email is filled if wesa is used
# TODO NEXT: Database
#            -  WESA queue
#                -  One function to send to wesa
#                -  One function (scheduled) to check if there are pending jobs
#                -  Update fields in tables job and Sequence to control pending
# TODO timesheet - Monday, 14:00 -
# DONE:
#            -  Save to database via celery ok
#            -  Create jobs ok
#            -  save the submission ok
#            -  save the results ok
#            -  make the results accessisble ok
# Future: celery has a 'chain' function that might be handy to pipe the
#    homology -> docking -> md process
# Future: start celery as a daemon:
#  http://docs.celeryproject.org/en/latest/userguide/daemonizing.html#daemonizing
# Future: use the timestamp to protect the results
# Future: read about security: https://docs.djangoproject.com/en/2.1/topics/security/#user-uploaded-content-security
# Future: check a service to verify if the email belongs to an academic institution


# TODO: put this function somewhere appropriate:
def handle_uploaded_file(f, filename):
    """Necessary for receiving files from the user (POST method)."""

    with open(filename, 'wb+') as destination:
        for chunk in f.chunks():
            destination.write(chunk)


def index(request):
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

            result = wf.workflow(jobfolder=jobfolder,
                                 rawseq=rawseq,
                                 fastafile=fastafile,
                                 job=myjob,
                                 mode=form.cleaned_data["mode"]
                                 )

            myjob.set_full_hbm_analysis([x["analysis"] for x in result])

            # if all went well, send the email to the user
            # https://docs.djangoproject.com/en/2.1/topics/email/#django.core.mail.EmailMessage
            email = form.cleaned_data['email']
            if email:
                body = "You can access the analysis at http://localhost:8000/SeqDHBM/%d" % myjob.id
                e_msg = EmailMessage(
                    subject=f'SeqD-HBM: Your analysis number {myjob.id}',
                    body=body,
                    from_email='seqdhbm@gmail.com',
                    to=[email],
                    headers={'Message-ID': 'foo'},
                )
                e_msg.send(fail_silently=False)
            return redirect("%d/" % myjob.id)
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


def show_result(request, job_id):
    template = loader.get_template('SeqDHBM/result.html')
    job = get_object_or_404(models.Job, id=job_id)
    seqs = get_list_or_404(models.Sequence, jobnum=job)
    res = {}
    for seq in seqs:
        res[seq] = [seq.seqchain[x:x+70] for x in range(0, len(seq.seqchain), 70)]
    context = {
        "job": job,
        "processing": len([x for x in seqs if x.status_hbm == models.Sequence.STATUS_QUEUED]),
        "result": {x: y for x, y in res.items() if x.status_hbm == models.Sequence.STATUS_PROCESSED},
        "failed": {x: y for x, y in res.items() if x.status_hbm == models.Sequence.STATUS_FAILED}
    }
    return HttpResponse(template.render(context, request))


def show_analysis(request, job_id):
    job = get_object_or_404(models.Job, id=job_id)
    return HttpResponse(job.full_hbm_analysis, content_type="text/plain")


def hemewf(request):
    message = ""
    for seq in models.Sequence.objects.all():
        message+="<div>"
        message+="<p>Job num: %d</p>"%seq.jobnum.id
        message+="<p>Sequence: %s</p>"%seq.seqchain
        message+="<p>Sent by: %s</p>"%seq.submittedas
        message+="<p>Prediction: %s</p>"%seq.mode
        message+="<p>Status of motif detection: %s</p>"%seq.status_hbm
        message+="<p>file location: %s</p>"%seq.fasta_file_location
        message+="<p>Warnings: %s</p>"%("<br>".join(seq.warnings_hbm.split("\n")))
    message+="<h1>Results</h1>"
    for res in  models.Result_HBM.objects.all():
        message+="<p>Sequence num: %d</p>"%res.sequence.id
        message+="<p>Coord atom: %s</p>"%res.coord_atom
    return HttpResponse(message)
