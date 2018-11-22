from django.db import models

# Create your models here.
class Job(models.Model):
    submittedby = models.EmailField(max_length=254, blank=True)
    submission_date = models.DateTimeField(auto_now_add=True)
    # full_analysis = models.TextField(blank=True)

class Sequence(models.Model):
    # choices
    MANUAL_INPUT = "M"
    FASTA_FILE = "F"
    PDB_FILE = "P"
    PDB_ID = "I"
    SUBMISSION_FORMS = ((MANUAL_INPUT, "Manual Input"), (PDB_ID, "PDB Id"),
                        (FASTA_FILE, "Fasta File"), (PDB_FILE, "PDB File"))

    WESA_MODE = "W"
    STRUCTURE_MODE = 'S'
    PREDICTION = ((WESA_MODE, "wesa"), (STRUCTURE_MODE, "structure"))

    # processes: to_detect_hbm, to_homology, to_dock
    # status: skipped, queued, processed, failed
    RECEIVED_STATUS = '1'
    PROCESSED_STATUS = '2'
    FAILED_STATUS = '3'
    STATUS = ((RECEIVED_STATUS, "Received"),
              (PROCESSED_STATUS, "Processed"),
              (FAILED_STATUS, "Failed"))
    # Fields
    jobnum = models.ForeignKey(
        Job,
        on_delete=models.CASCADE,
        verbose_name="Part of the job"
    )
    seqchain = models.TextField()
    submittedas = models.CharField(max_length=1, choices=SUBMISSION_FORMS)
    mode = models.CharField(max_length=1, choices=PREDICTION)
    status = models.CharField(max_length=1, choices=STATUS)
    # file_location = models.FilePathField(path=BASE_DIR, allow_files=True, allow_folders=False)
    # warnings blank= True

class Result(models.Model):
    sequence = models.ForeignKey(Sequence, on_delete=models.CASCADE, verbose_name="Belongs to a sequence")
    coord_atom = models.CharField(max_length=6)
    ninemer = models.CharField(max_length=9)
    net_charge = models.SmallIntegerField()
    disulfide_possible = models.BooleanField()
