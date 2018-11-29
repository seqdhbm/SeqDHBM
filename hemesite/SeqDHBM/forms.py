from django import forms

class SeqSubmission(forms.Form):
    emaillabel = "Enter your email address to receive the complete \
    analysis of the submitted sequences. As the structure prediction \
    can take a long time (from 5 to 90 minutes), this is mandatory if\
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
    mode = forms.ChoiceField(label="WESA?:", choices=CHOICES, disabled= True, initial="structure", widget=forms.RadioSelect())

    def clean(self):
        cleaned_data = super().clean()
        msg = "No sequence was provided for analysis!"
        if not cleaned_data.get("rawseq") and\
           not cleaned_data.get("fastafiles"):
            self.add_error("rawseq", msg)
            self.add_error("fastafiles", msg)
        if cleaned_data.get("rawseq") and cleaned_data.get("rawseq")[0]==">":
         self.add_error("rawseq", "Warning: Please remove the fasta header")
         self.add_error("rawseq", "Note: For processing multiple sequences, please use the file upload functionality")
        # i can check if the email is empty for WESA here
        if (not cleaned_data.get("email")) and cleaned_data.get("mode") == "wesa":
            self.add_error("email", "Email is mandatory when using WESA predictions")
