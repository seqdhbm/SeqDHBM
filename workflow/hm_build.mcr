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
MacroTarget='/home/imhof_team/Public/mauricio/IL-1_homology/IL-1_formatted'

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

# Normally no change required below this point

# Sanity checks
if MacroTarget==''
  RaiseError "This macro requires a target. Either edit the macro file or click Options > Macro > Set target to choose a target sequence"

Clear

# Do we have a target sequence? 
command1='pass'
seqfilename='(MacroTarget).fasta'
seq = FileSize (seqfilename) 
if seq
  # Yes, create a command to specify the sequence file
  command1='SequenceFile (seqfilename)'
  
# Do we have alignments?
alifilename='(MacroTarget)_align.fasta'
ali = FileSize (alifilename) 
if ali
  # Yes, use the alignments instead of the target sequence
  command1='AlignFile (alifilename)'
  
if !seq and !ali
  RaiseError 'Neither a sequence file (seqfilename), nor an alignment file (alifilename) was found, no model can be built'  
  
# Do we have templates?
command2='pass'
for i=001 to 999
  for ext in 'yob','pdb'
    tmpfilename='(MacroTarget)_t(i).(ext)'
    tmp = FileSize (tmpfilename)
    if tmp
      # Load template and name it T001.. 
      tmp = Load(ext) (tmpfilename)
      NameObj (tmp),T(i)
      # Create command to select all templates
      command2='TemplateObj all'
      # Template found, try next
      continue 2
  break

# Build the model
Experiment HomologyModeling
  # Specify either a target sequence or an alignment file
  (command1)
  # Select any templates that have been loaded
  (command2)
  # Set parameters defined at the beginning
  # Number of PSI-BLAST iterations
  PsiBLASTs (psiblasts)
  # The maximum PSI-BLAST Evalue (minimum can be set to avoid easy models)
  EValue Max=(evalue),Min=0
  # Maximum oligomerization state
  OligoState (oligostate)
  # Maximum number of templates to consider, and how many with the same sequence
  Templates (templates),SameSeq=1
  # Maximum number of ambiguous alignments to consider per template
  Alignments (alignments)
  # Maximum number of terminal loop residues
  TermExtension (termextension)
  # Speed of animations: fast/normal/slow
  Animation Normal
  # Number of samples to try per loop (see OptimizeLoop command)
  LoopSamples 50
  # Use structure-based profiles
  StructProfile (structprofile)
  # The common start for result filenames
  ResultFile (MacroTarget)
  # Uncomment below if you do not want to access PDB_REDO
  # PDBREDO No
  # Uncomment below to select certain template residues for deletion, e.g. het-groups
  # DelTemplateRes Hetgroup
  # Uncomment below to rank certain templates first (only affects templates found by YASARA)
  # TemplateList 2ZPY,2HE7,2I1K
  # Uncomment below to specify a list of excluded templates
  # TemplateExList 1CRN,5TIM,1AON
  # Uncomment below to add alignment anchors, see Recipes > Build a homology model > Useful hints
  # Equivalence Target=DHFGRER,Template=FDGFHGR
  # Uncomment below to keep the side-chains of selected conserved residues fixed, e.g. target residues 90 and 105
  # FixModelRes 90 105
  # Uncomment below to keep the side-chains of all conserved residues fixed (only useful if seqID near 100%).
  # FixModelRes conserved
  # Uncomment below to accept only templates that cover the selected target residues, e.g. residues 95 and 110
  # RequireRes 95 110
  # Uncomment below to skip loops that are longer than 9 residues (models will have gaps, no hybrid model built)
  # LoopLenMax 9
  # Uncomment below to build models with reduced accuracy as quickly as possible
  # Speed fast
  # Uncomment below to use your own secondary structure prediction
  # SecStrFile YourModel_secstr.ali
Experiment On
Wait ExpEnd

if runWithMacro and ConsoleMode
  # Exit YASARA if this macro was provided as command line argument in console mode
  Exit
if !ConsoleMode
  # Show homology modeling report
  ShowURL file://(MacroTarget).html
