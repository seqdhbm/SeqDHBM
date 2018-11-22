from django.http import HttpResponse
from django.template import loader
from django.core.mail import EmailMessage
from .forms import SeqSubmission
import seqdhbm.workflow as wf
import SeqDHBM.models as mod
from django.conf import settings
import os

# TODO read about security: https://docs.djangoproject.com/en/2.1/topics/security/#user-uploaded-content-security
# TODO fix a bug in fasta.fasta_to_seq2 when the user gives an non fasta file. ok
# TODO should I have a webpage with the results? can i redirect back to the main page on error?
# TODO NEXT: Database
#            -  Create jobs ok
#            -  save the submission
#            -  WESA queue
#            -  save the results
#            -  make the results accessisble
# TODO check a service to verify if the email belongs to an academic institution
# worklog: Fri Nov16th: 9:30 - 12:00; 13:00 - 18:30
#          Thu Nov23rd: 10:00 - 12:00; 13:00 - 18:30

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
            # Create the job
            email = form.cleaned_data['email'] if form.cleaned_data['email'] else ""
            myjob = mod.Job(submittedby= email)
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
            rawseq = form.cleaned_data['rawseq']
            rawseq = "".join(rawseq.split("\n"))
            rawseq = "".join(rawseq.split())
            rawseq = "".join(rawseq.split("\r"))
            result = wf.workflow(jobfolder=jobfolder, rawseq=rawseq, fastafile=fastafile)
            # if all went well, send the email to the user
            # https://docs.djangoproject.com/en/2.1/topics/email/#django.core.mail.EmailMessage
            email = form.cleaned_data['email']
            if email:
                body = "Body goes here"
                e_msg = EmailMessage(
                subject='SeqD-HBM: Your job number %(job)s',
                body=body,
                from_email='seqdhbm@gmail.com',
                to=[email],
                # connection: An email backend instance. Use this parameter
                # if you want to use the same connection for multiple messages.
                # If omitted, a new connection is created when send() is called.
                headers={'Message-ID': 'foo'},
                )
                e_msg.send(fail_silently=False)
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

def hemewf(request):
    message=""
    rawseq  = request.GET.get('aaseq')
    result_list = wf.workflow(rawseq=rawseq)#fastafile="/home/imhof_team/Public/mauricio/workflow/test.fasta")#
    try:
        for result in result_list:
            message += "<h3>%s</h3>"%result["name"]
            message += "<font face= 'Courier New'>"
            for i in range(0, len(result["seq"]), 70):
                message += "<p>%s</p>"%result["seq"][i:i+70]
            message += "</font>"
            for warn in result["warnings"]:
                message += "<p>%s</p>"%warn
            if result["result"]:
                message += "<table>"
                message += """<tr>
                        <th>Num</th>
                        <th>Coord. residue</th>
                        <th>9mer motif</th>
                        <th>Net charge</th>
                        <th>Comment</th></tr>"""
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
