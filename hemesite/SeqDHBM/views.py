# coding: utf-8

"""Functions that handles requests from users"""

from __future__ import absolute_import, unicode_literals

from django.conf import settings
from django.http import HttpResponse, Http404
from django.shortcuts import get_object_or_404, get_list_or_404, redirect
from django.template import loader

from .forms import SeqSubmission
from .views_utils import process_form
from SeqDHBM.models import Job, Result_HBM, Sequence

# DONE: start celery as a daemon:
# http://docs.celeryproject.org/en/latest/userguide/daemonizing.html#daemonizing
# TODO Validate if the file is fasta before going on
# TODO bottom bar with external links, references (eg: WESA), disclaimer
# TODO Create a github for the app
# TODO Docstrings
# TODO Replace localhost for IP
# TODO unit testing!
# TODO Upgrade the software to use yasara
# TODO Read papers from the Bioinformatics Journals
# TODO setup a different sgbd/web server (GUnicorn / apache)
# TODO Model a db for the group's sequences
# Future: celery has a 'chain' function that might be handy to pipe the
#    homology -> docking -> md process


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
        if form.is_valid():
            job_id, password = process_form(form, request.FILES)
            return redirect("%d/%s" % (job_id, password))
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
    THEN: List the number of pending jobs as well as html information for the
    processed jobs.

    :param request: the http request.
    :param job_id: The identifier of the job.
    :param passw: Security code to the given job
    :return: The formatted webpage with the requested information. The template
    comes from one http file.
    """
    template = loader.get_template('SeqDHBM/result.html')
    job = get_object_or_404(Job, id=job_id)
    seqs = get_list_or_404(Sequence, jobnum=job)
    if job.pass_gen() == passw:
        res = {}
        for seq in seqs:
            res[seq] = [seq.seqchain[x:x+70] 
                        for x in range(0, len(seq.seqchain), 70)]
        context = {
            "job": job,
            "processing": len([x for x in seqs 
                               if x.status_hbm == Sequence.STATUS_QUEUED]),
            "result": {x: y for x, y in res.items() 
                       if x.status_hbm == Sequence.STATUS_PROCESSED},
            "failed": {x: y for x, y in res.items() 
                       if x.status_hbm == Sequence.STATUS_FAILED},
            "passw": passw
        }
        return HttpResponse(template.render(context, request))
    else:
        raise Http404("Invalid URL (bad secret code and job id combination)")


def show_analysis(request, job_id, passw: str = ""):
    """
    WHEN: The user requests information about a job he submitted.
    GIVEN: The job exists.
    THEN: The full analysis of the job is presented to the user as a plain
    text webpage.

    :param request: http request
    :param job_id: The identifier of the job
    :param passw: Security code to the given job
    :return: The formatted webpage with the requested information.
    """
    job = get_object_or_404(Job, id=job_id)
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

    if not settings.DEBUG:
        return Http404()
    message = ""
    job = Job.objects.get(pk=job_id)
    for seq in Sequence.objects.filter(jobnum=job):
        message += "<h1>Sequence</h1>"
        message += "<div>"
        message += "<p>Job num: %d</p>" % seq.jobnum.id
        message += "<p>Sequence: %s</p>" % seq.seqchain
        message += "<p>Sent by: %s</p>" % seq.submittedas
        message += "<p>Prediction: %s</p>" % seq.mode
        message += "<p>Status of motif detection: %s</p>" % seq.status_hbm
        message += "<p>file location: %s</p>" % seq.fasta_file_location
        the_warnings = '<br>'.join(seq.warnings_hbm.split('\n'))
        message += f"<p>Warnings: {the_warnings}</p>"
        message += "<h1>Results</h1>"
        for res in Result_HBM.objects.filter(sequence=seq):
            message += "<p>Sequence num: %d</p>" % res.sequence.id
            message += "<p>Coord atom: %s</p>" % res.coord_atom
    # return HttpResponse(message)
    message += "\n\n\n"
    for seq in Sequence.objects.filter(
            jobnum=Job.objects.get(pk=job_id)):
        # message += seq.partial_hbm_analysis + "\n"
        ...
    job = Job.objects.get(pk=job_id)
    if job.pass_gen() == passw:
        print(job.id)
        print()
        # message += f"\n{request.META}"
    return HttpResponse(message)  # , content_type="text/plain")
