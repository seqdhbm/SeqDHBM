from django import forms

class SeqSubmission(forms.Form):
    emaillabel = "Enter a valid email address to receive the complete \
    results for the submitted sequence(s). Email entry is MANDATORY in \
    cases where WESA solvent accessibility predictions are requested"
    email = forms.EmailField(label=emaillabel,
                             required=False)
    # pdbfiles = forms.FileField(label="PDB files")
    fastafiles = forms.FileField(label="Select a local fasta file. Fasta files can contain multiple \
        sequences placed one below the other provided individual sequences are separated by valid fasta headers",
                                 required=False)
    rawseq = forms.CharField(label='Enter your sequence in upper case 1-letter code without fasta headers. \
        Only the 20 standard amino acids are accepted',
                             widget=forms.Textarea(attrs={'cols': '80', 'rows':'5'}),
                             required=False)
    CHOICES = [('structure','No, don\'t run WESA'),
               ('wesa','Yes, run WESA')]
    message = forms
    mode = forms.ChoiceField(label="Run WESA solvent accessibility prediction to discard buried residues ?", 
                            choices=CHOICES, initial="structure", widget=forms.RadioSelect())
    
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
