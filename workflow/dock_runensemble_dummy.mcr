# YASARA MACRO
# TOPIC:       5. Structure prediction
# TITLE:       Docking a ligand to a receptor ensemble with flexible side-chains 
# REQUIRES:    Structure
# AUTHOR:      Elmar Krieger
# LICENSE:     GPL
# DESCRIPTION: This macro predicts the structure of a ligand-receptor complex. Receptor flexibility is considered by creating a receptor ensemble with alternative high-scoring solutions of the side-chain rotamer network. An analysis log file is written at the end. The macro can also continue a docking run that got interrupted, especially if the scene has not been moved or rotated manually during the docking. 

# Parameter section - adjust as needed, but NOTE that some changes only take effect
# if you start an entirely new docking job, not if you continue an existing one. 
# =================================================================================

# You can either set the target structure by clicking on Options > Macro > Set target,
# by providing it as command line argument (see docs at Essentials > The command line),
# or by uncommenting the line below and specifying it directly.
#MacroTarget='/home/imhof_team/Public/mauricio/IL-1_homology/a'

# Docking method, either AutoDockLGA or VINA
method='VINA'

# Size of the receptor ensemble with alternative side-chain rotamers
receptors=5

# Number of docking runs per receptor ensemble member, '0' leaves the choice to YASARA.
# (maximally 999, each run can take up to an hour).
# The total number of docking runs is receptors*receptorruns.
receptorruns=0

# Docking results usually cluster around certain hot spot conformations,
# and the lowest energy complex in each cluster is saved. Two complexes belong to
# different clusters if the ligand RMSD is larger than this minimum [A]:
rmsdmin=5.0

# Set to 0 if you don't want to merge clusters from the receptor ensemble members
# in the end (merging keeps only one cluster within the 'rmsdmin' limit set above).
merged=1

# Set to 1 to keep the ligand completely rigid (alternatively you can provide
# the ligand as a *.yob file and fix certain dihedral angles only).
rigid=0

# A selection of receptor residues to keep flexible, e.g. flexres = 'Lys 91 Leu 100',
# or (if the receptor is not a monomer) flexres = 'Res Lys 91 Mol A or Res Leu 100 Mol B'.
# Alternatively you can provide the receptor as a *.sce or *.yob file with fixed
# atoms, which gives better control (e.g. you can keep only part of a side-chain
# or even a terminal backbone flexible). 
flexres='' 
flexres = ''

# Force field used for charge assignment, not used by VINA.
ForceField AMBER03

# The docking temperature. This currently affects the receptor flexibility only
temp='298K' 

# Flag if a scene with all results should be saved for visualization with dock_play.mcr.
# If you use a large receptor ensemble with large ligands and many docking runs, it may
# exhaust the available memory to show them all together on screen. In this case set
# scesaved to 0. The docking result player will then only show clusters but not all results,
# which can of course still be found in the logfile.
scesaved=1 


# Normally no change required below this point
# ============================================

chargefof = ForceField
if !receptorruns
  # For an optimum speed/accuracy tradeoff, we select at least 20 runs, but more if CPUs are available
  cpus = Processors
  receptorruns=cpus
  while receptorruns<20
    receptorruns=receptorruns+cpus

# Sanity checks
if MacroTarget==''
  RaiseError "This macro requires a target. Either edit the macro file or click Options > Macro > Set target to choose a target structure"
if receptors>999 or receptorruns>999 or receptors*receptorruns>9999
  RaiseError 'Too many docking runs selected, (receptors)/(receptorruns) would take forever'

structlist='receptor','ligand'

Clear
# Do we have a user-supplied ensemble?
ensfilename='(MacroTarget)_receptor_ensemble.sce'
scene = FileSize (ensfilename)
if !scene
  # Do we already have an automatically generated receptor ensemble scene?
  ensfilename='(MacroTarget)_receptor_ensemble(receptors).sce'
scene = FileSize (ensfilename)
if !scene
  # Not yet. Do we have a receptor scene?
  scefilename='(MacroTarget)_receptor.sce'
  scene = FileSize (scefilename)
  if !scene
    # No, load PDB or YOB file of receptor and ligand
    for struct in structlist
      (struct)=0
      for type in 'yob','pdb','sdf'
        filename='(MacroTarget)_(struct).(type)'
        exists = FileSize (filename)
        if exists and (type!='sdf' or struct=='ligand')
          (struct) = Load(type) (filename)
          break
      if not (struct)
        RaiseError 'The (struct) was not found at (filename)'
    # Orient the receptor, create a docking cell that includes the entire receptor 
    # and also has enough empty space around to accomodate the ligand (which has to
    # be removed temporarily, since 'Cell Auto' encloses the entire soup). 
    ligandradius = RadiusObj (ligand)  
    NiceOriObj (receptor)
    # Ligand was only needed to determine the required cell size, will be loaded later
    DelObj (ligand)
    Cell Auto,Extension=(ligandradius*2)
  else
    # Load receptor scene
    LoadSce (scefilename)
    # Verify that the cell is present
    simcell = CountObj SimCell
    if !simcell
      RaiseError 'If you provide a scene, it must contain a simulation cell, but none was found in (scefilename)'
    if Objects>2
      # We can't continue, because the ensemble can only be created for a single object
      RaiseError 'Your scene contains (Objects) objects, while only the receptor and the docking cell are expected.'
  # Warn if the cell is too large
  x,y,z = Cell
  if x>96 or y>96 or z>96
    ShowMessage "A cell axis is longer than 96 A, which may reduce docking accuracy. Consider focusing on the active site."
    Wait ContinueButton
    HideMessage
  # To create a receptor ensemble, the entire receptor must be inside the cell,
  # so we remember the current cell and add it again later
  x,y,z = Cell
  posorilist() = PosOriObj SimCell
  DelObj SimCell
  # Start an energy minimization but stop immediately.
  # This cleans the structure in case it was supplied without hydrogens.
  ForceField (chargefof)
  ShowMessage 'Checking if receptor is ready for docking...'
  Experiment Minimization
  Experiment On
  Experiment Off
  # Create a receptor ensemble with flexibility corresponding to the specified temperature
  ShowMessage 'Creating receptor ensemble with (receptors) members, each with alternative high-scoring side-chain conformations...' 
  Style Ribbon,Stick
  ForceField YASARA2
  Temp (temp)
  # Don't change side-chain rotamers of residues with dative bonds to metals
  FixRes all with arrow to metal
  OptimizeAll Method=SCWALL,Structures=(receptors-1)
  FreeAll
  DelObj SimCell
  NumberObj all,1
  # Add original cell again
  Cell (x),(y),(z)
  PosOriObj SimCell,(posorilist)
  # Optimize the hydrogen bonding network of the ensemble members
  RemoveObj !SimCell 
  for i=1 to receptors
    AddObj (i)
    ShowMessage 'Optimizing hydrogen bonding network of receptor ensemble member (i)/(receptors)...'
    OptHydAll
    NameObj (i),Receptor(i)
    RemoveObj (i)
  AddObj !SimCell
  SaveSce (ensfilename)
  HideMessage
# Load receptor ensemble scene
LoadSce (ensfilename)
# The segment name C*** is used to tag ligand Conformations, and cannot be used for the receptor
hit = ListAtom Segment C???
if hit
  MarkAtom (hit)
  segname = SegAtom (hit)
  RaiseError 'Segment names starting with "C" are unfortunately reserved to identify ligand conformations, please click "Edit > Rename > Segment" to rename the receptor segment "(segname)"' 
# The molecule name 'l' is reserved for the ligand
hit = ListAtom Mol l
if hit
  RaiseError 'The receptor must not contain a molecule named "l", please click Edit > Rename > Molecule and try again'
# Docking is done without periodic boundaries
Boundary Wall
Longrange None  
ForceField (chargefof)
# Count receptors, in case the user provided the ensemble
receptors = CountObj !SimCell
runs=receptors*receptorruns
# Loop over the receptors in the ensemble
RemoveObj !SimCell
for i=1 to receptors
  AddObj (i)
  # Show progress (ShowMessage is used by the docking experiment)
  LabelAll 'Receptor ensemble member (i)/(receptors)',Height=0.3,Color=Yellow,Y=3.0
  # Load the ligand for this receptor ensemble member
  ligand=0
  for type in 'yob','pdb'
    filename='(MacroTarget)_ligand.(type)'
    exists = FileSize (filename)
    if exists 
      ligand = Load(type) (filename)
      break
  if not ligand
    RaiseError 'The ligand was not found at (filename)'
  # Just in case a user-defined ensemble was loaded
  NameObj (i),Receptor(i)
  # Number the ligands sequentially
  NameObj (ligand),Ligand(i)
  # Keep selected side-chains flexible
  if flexres!=''
    FixObj (i)
    FreeAtom Res (flexres) Sidechain and Obj (i)
    FixRes Cys Atom SG with bond to Atom SG or Ala Pro and Obj (i)
  # Do not show secondary structure for protein ligands, this slows things down
  HideSecStrObj Ligand(i)
  ShowObj Ligand(i)
  StickObj Ligand(i)
  # Fix ligand's internal degrees of freedom if requested 
  if rigid
    FixObj Ligand(i)
  # Align ligand with the cell (to end up with the same local atom coordinates independent
  # of cell orientation, which are used to calculate the checksum for the *.adr file) 
  TransferObj Ligand(i),SimCell,Local=Keep
  # Perform the docking
  # DEFINE HERE THE CELL!!
  Cell Auto,extension=17, Shape=Cuboid, Res Tyr 29
  Experiment Docking
    Method (method)
    ReceptorObj (i)
    LigandObj Ligand(i)
    Runs (receptorruns)
    ClusterRMSD (rmsdmin)
    # Result file number must be separated with '_', in case MacroTarget also ends with number
    ResultFile (MacroTarget)_(i)_001
    # Uncomment below to set the number of energy evaluations (AutoDock ga_num_evals):
    # DockPar ga_num_evals 25000000
    # Uncomment below to add 1000 to the random number seed when more than 999 runs are needed
    # DockPar seed 1000
    # Uncomment the two lines below to provide your own atom parameters (VdW radii etc.)
    # GridPar parameter_file /Path/To/Custom/AD4_parameters.dat
    # DockPar parameter_file /Path/To/Custom/AD4_parameters.dat
    # Uncomment below to keep temporary files in the current working directory:
    # TmpFileID 1adb_tmp
  Experiment On
  Wait ExpEnd
  UnlabelAll
  RemoveObj (i) Ligand(i) RecepFlex
  NameObj RecepFlex,RecFlex(i)
if scesaved
  # Save a scene with all receptor and all ligand conformations
  AddObj all
  SaveSce (MacroTarget)
  RemoveObj !SimCell

# Collect the initial results
Console off
for i=1 to receptors
  AddObj (i) Ligand(i)
  for j=001 to receptorruns
    # To pass the docking results back to this macro, the docking experiment
    # stores the binding energy as the B-factor...
    run=(i-1)*receptorruns+j
    resultlist(run),bindnrglist(run) = DockingResult i,j,'Obj (i)','Obj Ligand(i) Segment C(j)'
    if !(run%10)
      ShowMessage 'Analyzing docking results, (100*run/runs)% completed...'
      Wait 1
  RemoveObj (i) Ligand(i)
# Sort the results by binding energy
ShowMessage 'Sorting docking results...'
Wait 1
SortResults 'bindnrglist','resultlist',''
# Save a log file with an initial analysis
RecordLog (MacroTarget)
print 'Ensemble docking result analysis'
print '================================'
print
print 'The ligand was docked (receptorruns) times ["Run"] with (method) against each of the (receptors) receptors ["Rec"] in the ensemble, yielding the following (runs) results, sorted by binding energy:'
print '[More positive energies indicate stronger binding, and negative energies mean no binding]'
print
print 'Rec | Run |Bind.energy[kcal/mol]|Dissoc. constant [pM]| Contacting receptor residues'
print '----+-----+---------------------+---------------------+-----------------------------'
for i=1 to runs
  print (resultlist(i))
print

# Get cell dimensions, needed to move the atoms in the cluster YOb files
cd1,cd2,cd3=Cell
# Delete previous resultlist and bindnrglist
bindnrglist()=0
resultlist()=0
# Now merge the clusters from the receptor ensemble members
clusters=0
for i=1 to receptors
  for j=001 to receptorruns
    filename='(MacroTarget)_(i)_(j).yob'
    exists = FileSize (filename)
    if exists
      ShowMessage 'Merging clusters from receptor ensemble member (i) run (j), (clusters) found so far...'
      Wait 1
      cluobj = LoadYOb (filename)
      UnlabelAll
      # Place the object in the dock cell to aid visualization 
      MoveAtom Obj (cluobj),(-cd/2)
      TransferObj (cluobj),SimCell,Local=Keep  
      result,bindnrg = DockingResult i,j,'Obj (cluobj) Segment !C???','Obj (cluobj) Segment C???'
      added=1
      # Keep only the ligand to save time and memory, otherwise we can quickly exceed one million atoms on screen
      DelRes Segment !C??? Obj (cluobj)
      if clusters and merged
        # Get ligand RMSDs from all other clusters
        rmsdlist() = RMSDAtom Element !H Segment C??? Obj (cluobj),Element !H Obj (join cluobjlist)  
        for k=1 to clusters
          if rmsdlist(k)<rmsdmin
            # Another cluster is close, remember the one with the better binding energy
            added=0
            if bindnrg>bindnrglist(k)
              # Current one has better binding energy, delete old one
              resultlist(k)=result
              bindnrglist(k)=bindnrg
              filenamelist(k)=filename
              SwapObj (cluobjlist(k)),(cluobj)
            DelObj (cluobj)
            break
      if added
        # Found a new cluster
        clusters=clusters+1
        resultlist(clusters)=result
        bindnrglist(clusters)=bindnrg
        cluobjlist(clusters)=cluobj
        filenamelist(clusters)=filename
ShowMessage 'Sorting cluster results...'
Wait 1
DelObj (join cluobjlist)
SortResults 'bindnrglist','resultlist','filenamelist'
# Save the unique clusters in order, mainly needed for the dock_play.mcr
for i=1 to clusters
  obj = LoadYOb (filenamelist(i))
  # Transfer again into the cell
  MoveAtom Obj (obj),(-cd/2)
  TransferObj (obj),SimCell,Local=Keep  
  SaveYOb (obj),(MacroTarget)_(000+i).yob
  ShowMessage 'Saving sorted cluster (i)/(clusters)...'
  Wait 1
  DelObj (obj)
# Delete original clusters on disk
for i=1 to receptors
  for j=001 to receptorruns
    DelFile (MacroTarget)_(i)_(j).yob
if scesaved
  # Bring original objects back, unless that would exhaust memory
  AddAll

# Print cluster log
print 'After clustering the (runs) runs, the following (clusters) distinct complex conformations were found:'
print '[They all differ by at least (rmsdmin) A heavy ligand atom RMSD, and have been saved as YOb files with terminal number "Num"]' 
print
print 'Num | Rec | Run |Bind.energy[kcal/mol]|Dissoc. constant [pM]| Contacting receptor residues'
print '----+-----+-----+---------------------+---------------------+-----------------------------'
for i=1 to clusters
  print '(000+i) | (resultlist(i))'
Style Ribbon
StopLog
Console On
HideMessage
# Exit YASARA if this macro was provided as command line argument in console mode
if runWithMacro and ConsoleMode
  Exit

# EXTRACT DOCKING RESULT
# ======================
# Returns a result string and binding energy, extracted from atomic BFactor and Property.
# 'rec' is the receptor number in the ensemble, 'run' is the docking run for this receptor, 
# 'recsel' selects the receptor, 'ligsel' selects the ligand.
def DockingResult rec,run,recsel,ligsel
  # Get the binding energy from the B-factor...
  bindnrg = BFactorAtom (ligsel)
  # ...and the dissociation constant from the atomic property.
  dissconst = PropAtom (ligsel)
  if bindnrg<=0
    dissconst='               None'
  else
    dissconst= 00000000000000.0000+dissconst
  # Get receptor residues that contact the ligand
  reslist() = ListRes (recsel) with distance<4.0 from (ligsel),Format='RESNAME MOLNAME RESNUM'
  if !count reslist
    reslist='No residues in contact'
  else
    reslist=join reslist
  result='(000+rec) | (000+run) |         (000000.0000+bindnrg) | (dissconst) | (reslist)'
  return result,bindnrg

# SORT LISTS BY VALUES, HIGHEST FIRST
# ===================================
# Lists are passed by reference, see Yanaconda doc section 'Calls to user-defined functions..'
def SortResults vallist,key1list,key2list
  caller (vallist),(key1list),(key2list)

  elements=count (vallist)
  # Create a list of indices to sort
  for i=1 to elements
    indexlist(i)=i
  # Sort the list of indices  
  passes=0
  do
    sorted=1
    for i=1 to elements-1
      if (vallist)(i)<(vallist)(i+1)
        swap=(vallist)(i)
        (vallist)(i)=(vallist)(i+1)
        (vallist)(i+1)=swap
        swap=indexlist(i)
        indexlist(i)=indexlist(i+1)
        indexlist(i+1)=swap
        sorted=0
    if elements>250 and not passes%10
      ShowMessage 'Sorting results, pass (passes)/(elements)..'
      Wait 1
    passes=passes+1
  while not sorted
  # Rebuild the sorted keylists
  for i=1 to 2
    if key(i)list!=''
      for j=1 to elements
        keylist(j)=(key(i)list)(indexlist(j))
      (key(i)list)()=keylist
  HideMessage

