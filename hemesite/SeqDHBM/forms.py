# coding: utf-8

from django import forms
from django.conf import settings


class SeqSubmission(forms.Form):
    """Django form to submit data for analysis"""

    emaillabel = "Enter your email address to receive the complete \
    analysis of the submitted sequences. As the structure prediction \
    can take a long time (from 5 to 90 minutes). This is mandatory if\
    you want the structure to be predicted:"
    email = forms.EmailField(label=emaillabel,
                             required=False)
    # pdbfiles = forms.FileField(label="PDB files")
    fastafiles = forms.FileField(label="Fasta files",
                                 required=False)
    rawseq = forms.CharField(label='Enter or paste your sequence data (one sequence with no fasta header):',
                             widget=forms.Textarea(attrs={'cols': '80', 'rows':'5'}),
                             required=False)
    CHOICES = [('structure','Skip'),
               ('wesa','Predict')]
    mode = forms.ChoiceField(label="WESA?:", choices=CHOICES, initial="structure", widget=forms.RadioSelect())

    def clean(self):
        """Checks if the form was filled correctly."""

        cleaned_data = super().clean()
        msg = "No sequence was provided for analysis!"
        big_file_msg = f"Please, keep file size under " +\
                       f"{settings.MAX_FILE_SIZE/1024} KB"
        if not cleaned_data.get("rawseq") and\
           not cleaned_data.get("fastafiles"):
            self.add_error("rawseq", msg)
            self.add_error("fastafiles", msg)
        if cleaned_data.get("fastafiles") and \
                self.cleaned_data["fastafiles"].size > settings.MAX_FILE_SIZE:
            self.add_error("fastafiles", big_file_msg)
        if cleaned_data.get("rawseq") and cleaned_data.get("rawseq")[0]==">":
            self.add_error("rawseq", "Warning: Please remove the fasta header")
            self.add_error("rawseq", "Note: For processing multiple sequences, please use the file upload functionality")
        if (not cleaned_data.get("email")) and cleaned_data.get("mode") == "wesa":
            self.add_error("email", "Email is mandatory when using WESA predictions")
