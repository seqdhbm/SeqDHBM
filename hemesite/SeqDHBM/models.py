from django.db import models
from django.conf import settings

# Create your models here.
class Job(models.Model):
    submittedby = models.EmailField(max_length=254, blank=True)
    submission_date = models.DateTimeField(auto_now_add=True)
    full_hbm_analysis = models.TextField(blank=True)

class Sequence(models.Model):
    # choices for form of submission
    SUB_MANUAL_INPUT = "M"
    SUB_FASTA_FILE = "F"
    SUB_PDB_FILE = "P"
    SUB_PDB_ID = "I"
    SUBMISSION_FORMS = ((SUB_MANUAL_INPUT, "Manual Input"), (SUB_PDB_ID, "PDB Id"),
                        (SUB_FASTA_FILE, "Fasta File"), (SUB_PDB_FILE, "PDB File"))

    # choices for structure prediction
    WESA_MODE = "W"
    STRUCTURE_MODE = 'S'
    PREDICTION = ((WESA_MODE, "wesa"), (STRUCTURE_MODE, "structure"))

    # status for each of the steps of the process
    STATUS_QUEUED = '1'
    STATUS_PROCESSED = '2'
    STATUS_FAILED = '3'
    STATUS_SKIPPED = '4'
    STATUS = ((STATUS_QUEUED, "Received"),
              (STATUS_PROCESSED, "Processed"),
              (STATUS_FAILED, "Failed"),
              (STATUS_SKIPPED, "Skipped"))
    # Fields
    jobnum = models.ForeignKey(
        Job,
        on_delete=models.CASCADE,
        verbose_name="Part of the job"
    )
    header = models.CharField(max_length=300)
    seqchain = models.TextField()
    submittedas = models.CharField(max_length=1, choices=SUBMISSION_FORMS)
    mode = models.CharField(max_length=1, choices=PREDICTION)
    status_hbm = models.CharField(max_length=1, choices=STATUS)
    partial_hbm_analysis = models.TextField(blank=True)
    # status_homology
    # status_docking
    # status_md_sim
    fasta_file_location = models.FilePathField(
        path=settings.BASE_DIR,
        allow_files=True,
        allow_folders=False
        )
    # pdb_file_location = models.FilePathField(path=BASE_DIR, allow_files=True, allow_folders=False)
    warnings_hbm = models.TextField(blank= True)

class Result_HBM(models.Model): # Result_HBM
    sequence = models.ForeignKey(
        Sequence,
        on_delete=models.CASCADE,
        verbose_name="Belongs to a sequence"
        )
    coord_atom = models.CharField(max_length=6)
    ninemer = models.CharField(max_length=9)
    net_charge = models.CharField(max_length=20)
    disulfide_possible = models.BooleanField()
    # spacer_check =
