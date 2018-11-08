from django.http import HttpResponse
from django.template import loader
from django.core.mail import EmailMessage


import seqdhbm.workflow as wf

def index(request):
    cont = 1#wf.workflow(fastafile="/home/imhof_team/Public/mauricio/workflow/ARO2(376)/ARO2(376).fasta" ,mode="structure")
    template = loader.get_template('SeqDHBM/index.html')
    context = {
        'variablename': [1,2,3],
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

    """
    # https://docs.djangoproject.com/en/2.1/topics/email/#django.core.mail.EmailMessage

    email = EmailMessage(
    subject='Hello',
    body='Body goes here',
    from_email='from@example.com',
    to=['to1@example.com', 'to2@example.com'],
    # connection: An email backend instance. Use this parameter
    # if you want to use the same connection for multiple messages.
    # If omitted, a new connection is created when send() is called.
    headers={'Message-ID': 'foo'},
    )
    email.send(fail_silently=False)
    """
    return HttpResponse((message if message else "You submitted nothing!"))
