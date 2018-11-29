from django.http import HttpResponse
from django.template import loader
from django.core.mail import EmailMessage
from .forms import SeqSubmission
import seqdhbm.workflow as wf
from SeqDHBM import models
from django.conf import settings
from django.shortcuts import get_object_or_404, get_list_or_404, redirect
import os

# TODO discuss with Ajay about how to protect the results
# Read about Celery
# # http://docs.celeryproject.org/en/latest/django/first-steps-with-django.html
# # http://docs.celeryproject.org/en/latest/getting-started/first-steps-with-celery.html#first-steps
# TODO read about security: https://docs.djangoproject.com/en/2.1/topics/security/#user-uploaded-content-security
# TODO fix a bug in fasta.fasta_to_seq2 when the user gives an non fasta file. ok
# TODO should I have a webpage with the results? can i redirect back to the main page on error?
# TODO NEXT: Database
#            -  Create jobs ok
#            -  save the submission ok
#            -  WESA queue
#            -  save the results ok
#            -  make the results accessisble ok
# TODO check a service to verify if the email belongs to an academic institution
# worklog: Fri Nov16th: 9:30 - 12:00; 13:00 - 18:30
#          Thu Nov22nd: 10:00 - 12:00; 13:00 - 17:00; 17:30 - 20:15
#          Fri Nov23rd: 9:15 - 13:00; 13:45 - 16:15; 16:45 -

# TODO: put this function somewhere appropriate:
def handle_uploaded_file(f, filename):
    with open(filename, 'wb+') as destination:
        for chunk in f.chunks():
            destination.write(chunk)

def index(request):
    result = []
    debugging=""
    # if this is a POST request we need to process the form data
    if request.method == 'POST':
        # create a form instance and populate it with data from the request:
        form = SeqSubmission(request.POST, request.FILES)

        debugging = form.non_field_errors()
        # check whether it's valid:
        if form.is_valid():
            email = form.cleaned_data['email'] if form.cleaned_data['email'] else ""
            # Create the job in the database
            myjob = models.Job(submittedby= email)
            myjob.save()
            jobfolder = os.path.join(settings.BASE_DIR, "jobs/J%07d/"%myjob.id)
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
            result = wf.workflow(jobfolder=jobfolder, rawseq=rawseq, fastafile=fastafile, job=myjob)
            print(result)

            ANALYSIS_HEADER = "*"*100
            ANALYSIS_HEADER += "\nSeqD-HBM : [Seq]uence based [D]etection of [H]eme [B]inding [M]otifs\n"
            ANALYSIS_HEADER += myjob.submission_date.strftime("%A , %B-%d-%Y, %H:%M:%S")
            ANALYSIS_HEADER += "\nJob number %d\n"%myjob.id
            ANALYSIS_HEADER += "Full analysis report\n" +("*"*100)+ "\n\n"

            myjob.full_hbm_analysis = ANALYSIS_HEADER + "\n\n\n".join([x["analysis"] for x in result])
            myjob.save()

            # if all went well, send the email to the user
            # https://docs.djangoproject.com/en/2.1/topics/email/#django.core.mail.EmailMessage
            email = form.cleaned_data['email']
            if email:
                body = "You can access the analysis at http://localhost:8000/SeqDHBM/%d"%myjob.id
                e_msg = EmailMessage(
                    subject='SeqD-HBM: Your analysis number %(job)s'%myjob.id,
                    body=body,
                    from_email='seqdhbm@gmail.com',
                    to=[email],
                    # connection: An email backend instance. Use this parameter
                    # if you want to use the same connection for multiple messages.
                    # If omitted, a new connection is created when send() is called.
                    headers={'Message-ID': 'foo'},
                )
                e_msg.send(fail_silently=False)
            return redirect("%d/"%myjob.id)
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
    print (job)
    res = {}
    for seq in seqs:
        res[seq] = [seq.seqchain[x:x+70] for x in range(0, len(seq.seqchain), 70)]
    context = {
        "job": job,
        "result": {x:y for x,y in res.items() if x.status_hbm not in [models.Sequence.STATUS_FAILED, models.Sequence.STATUS_SKIPPED] },
        "failed": {x:y for x,y in res.items() if x.status_hbm == models.Sequence.STATUS_FAILED}
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


"""    message=""
    rawseq  = request.GET.get('aaseq')
    result_list = wf.workflow(rawseq=rawseq)#fastafile="/home/imhof_team/Public/mauricio/workflow/test.fasta")#
    try:
        for result in result_list:
            message += "<h3>%s</h3>"%result["name"]
            message += "<font face= 'Courier New'>"
            for i in range(0, len(result["seq"]), 70):
                message += "<p>%s</p>"%result["seq"][i:i+70]
            message += "</font>"
            if result["fail"]:
                message += "<font color='FF0000'><b>"
            for warn in result["warnings"]:
                message += "<p>%s</p>"%warn
            if result["fail"]:
                message += "</b></font>"
            if result["result"]:
                message += "<table>"
                message += ""<tr>
                        <th>Num</th>
                        <th>Coord. residue</th>
                        <th>9mer motif</th>
                        <th>Net charge</th>
                        <th>Comment</th></tr>""
                line_num = 0
                for coord_atom, atom_data in sorted(result["result"].items(), key = lambda x: int(x[0][1:])):
                    line_num+=1
                    message += "<tr>"
                    message += "<td>%d</td>"%line_num
                    message += "<td>%s</td>"%coord_atom
                    # TODO should I move the style to a css?
                    message += "<td><font face= 'Courier New'>%s</font></td>"%atom_data["ninemer"]
                    message += "<td>%s</td>"%atom_data["netcharge"]
                    message += "<td>%s</td>"%atom_data["comment"]
                    message += "</tr>"
                message += "</table>"
    except Exception as e:
        message += str(e)

    return HttpResponse((message if message else "You submitted nothing!"))
"""
