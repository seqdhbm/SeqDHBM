# YASARA MACRO
# TOPIC:       5. Structure prediction
# TITLE:       Building a homology model
# REQUIRES:    Structure
# AUTHOR:      Elmar Krieger
# LICENSE:     GPL
# DESCRIPTION: This macro builds a homology model using a FASTA sequence of the target, and optionally template structures and alignments

# Parameter section - adjust as needed, but NOTE that some changes only take effect
# if you start an entirely new modeling job, not if you continue an existing one. 
# =================================================================================

# You can either set the target structure by clicking on Options > Macro > Set target,
# by providing it as command line argument (see docs at Essentials > The command line),
# or by uncommenting the line below and specifying it directly.
# MacroTarget='/home/imhof_team/Public/mauricio/IL-1_homology/IL-1_formatted'

# Number of PSI-BLAST iterations
# Set to '1' to search the template with a simple BLAST, not a PSSM-based PDB-BLAST
psiblasts=3

# Maximum PSI-BLAST Evalue allowed for templates
evalue=0.5 

# Maximum number of templates to use (if you provide fewer templates, the others will be picked by YASARA) 
templates=5

# Maximum number of ambiguous alignments to consider per template
alignments=5

# Maximum oligomerization state, build at most tetrameric models
oligostate=4

# Maximum number of unaligned loop residues to add to the termini
termextension=10

# How template profiles are calculated. Set to 'no' for simple sequence profiles.
# If set to 'yes', the more accurate PSSP database is used, which is based on
# twisted structural alignments of related structures. Search docs for PSSP.
structprofile='try'

Name = 'Yami'
Print 'Hi I am (Name)!'
# Normally no change required below this point

# Sanity checks
if MacroTarget==''
  RaiseError "This macro requires a target. Either edit the macro file or click Options > Macro > Set target to choose a target sequence"

Exit
