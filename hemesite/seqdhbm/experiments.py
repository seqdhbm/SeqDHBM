
# coding: utf-8

# # Homology function

# In[3]:


imtesting = True # After done testing, remove every condition this appears


# In[16]:


# TODO: Will we ever need to change this configuration parameters?
# TODO: Set it to use the GPU and all CPU threads? yasara.Processors(cputhreads=None, gpu=None)
import logging
logging.basicConfig(level=logging.DEBUG, format='%(levelname)s - %(message)s')

import yasara as y
y.info.mode='txt'
def yasara_homology(filename,            in_psiblasts=3,          in_evalue=0.5,                     in_templates=5,      in_alignments=5,         in_oligostate=4,                    in_termextension=10, in_structprofile='try'):
    """"""
    logging.info("Homology experiment: " + filename)
    y.ExperimentHomologyModeling(sequencefile=filename,                                psiblasts=in_psiblasts,                                evalue=in_evalue,                                 templates=in_templates,                                 alignments=in_alignments,                                oligostate=in_oligostate,                                termextension=in_termextension,                                structprofile=in_structprofile)
    y.Experiment("On")
    y.Wait("ExpEnd")
    logging.info("Exp ended")
    return True


# # Function to Analyze the Quality of the homologies

# In[7]:


import re

h_one = re.compile("This hybrid model with Z-score (.{,10}) was saved as the final one, <.*?>(.*\.yob)")

single_z = re.compile("The model with an overall quality Z-score of (.{,10}) has been saved as ")
single_yob = re.compile("Instead, the model .*? was simply saved as the final one <.*?>(.*?yob).")

no_hyb = re.compile("The hybrid model was discarded, and (.*?yob) was saved as the final model")
no_hyb_table = re.compile("<table><tr><td>Transfer.*?</td></tr></table>")
no_hyb_z = re.compile("<b>(.*?)</b>")


def homology_qa(folder):
    """Obtain z-scores and the final yasara object for a simulation.
    
    """
    global imtesting
    file = folder
    if imtesting:
        folder = "Saved_run/"+folder
    with open(folder+'/'+file+'.html') as f:
        html = f.read() # load html as text
    html = html[html.find("The model ranking"):] # discard the initil analysis, keep the summary
    html = "".join(html.splitlines()) # remove newline symbols
    
    # Try to obtain results for single hybrid model
    re_obj = h_one.search(html)
    if re_obj:
        return re_obj.groups()
    re_obj = single_z.search(html)
    # Try to obtain results for single model
    if re_obj:
        re_obj1 = single_yob.search(html)
        if re_obj1:
            return (re_obj.group(1), re_obj1.group(1))
    re_obj = no_hyb.search(html)
    # Try to obtain results for a multi model tthat couldn't be hybridized
    if re_obj:
        table = no_hyb_table.search(html)
        if table:
            re_obj1 = no_hyb_z.search(table.group(0))
            if re_obj1:
                return (re_obj1.group(1), re_obj.group(1))
    # TODO: if nothing was found, warn the user
    logging.error("Could not parse file '"+folder+"' properly")
    return False

if __name__ == "__main__":
    print (homology_qa("IL-18(189)"))
    """print (homology_qa("IL-18(193)"))
    print (homology_qa("IL-1A(271)"))
    print (homology_qa("IL-1B(268)"))
    print (homology_qa("IL-1B(300)"))
    print (homology_qa("IL-1Ra(143)"))
    print (homology_qa("IL-1Ra(159)"))
    print (homology_qa("IL-1Ra(177)"))
    print (homology_qa("IL-1Ra(180)"))
    print (homology_qa("IL-33(144)"))
    print (homology_qa("IL-33(228)"))
    print (homology_qa("IL-33(270)"))
    print (homology_qa("IL-36Ra(155)"))
    print (homology_qa("IL-37(157)"))
    print (homology_qa("IL-37(178)"))
    print (homology_qa("IL-37(192)"))
    print (homology_qa("IL-37(197)"))
    print (homology_qa("IL-37(218)"))
    print (homology_qa("IL-38(152)"))
    print (homology_qa("IL-6(212)"))
    print (homology_qa("IL-8(99)"))
"""


# In[1]:


if __name__ == "__main__":
    seq_dict = fasta_to_seq("IL-1_formatted.fastaR")
    organize_sequences(seq_dict)
    for folder in seq_dict.keys():
        print(folder)
        yasara_homology(folder+'/'+folder+'.fasta')
        z_score, yas_obj = homology_qa(folder) # TODO: Store and/or print the z-score


# # Using a function to run the Dock Ensembl macro

# In[11]:


imtesting = True # After done testing, remove every condition this appears (keep any 'else' clause)
import yasara as y

y.info.mode='txt' 

# TODO: Adjust the parameters inside the file
# y.CellAuto(extension=17, shape="Cuboid", selection1="Obj 1 Res Pro 41")
# Cell Auto,extension=17, Shape=Cuboid, Res Tyr 29

import shutil
from os import rename
from os.path import isfile

def parametrize_macro(flexres="flexres = ''", celldef = ""):
    """ Creates a parametrized macro from the model."""
    macromodel = "dock_model.mcr"
    macro = "dock_runensemble_dummy.mcr"

    with open(macromodel) as f:
        macrotext = f.read()
    macrotext = flexres.join(macrotext.split("# $$FLEXRES$$"))
    macrotext = celldef.join(macrotext.split("# $$CELL$$"))
    with open(macro, "w") as f:
        f.write(macrotext)


def yasara_docking(folder, receptor_file, method='VINA'):
    global imtesting
    if imtesting:
        folder = "Saved_run/"+folder
    filename = ".".join(receptor_file.split(".")[:-1])
    fileext = receptor_file.split(".")[-1]
    ligmodel = "a_ligand.pdb"
    newlig = folder + "/"+ filename + "_ligand.pdb"
    
    
    parametrize_macro(celldef = "Cell Auto,extension=17, Shape=Cuboid, Res Tyr 29")
    
    # copies the ligand to the receptor's folder with the proper name
    if not isfile(newlig):
        shutil.copyfile(ligmodel, newlig)

    # renames the receptor's yasara object adding _receptor
    oldrec = folder + "/" + receptor_file
    newrec = folder + "/" + filename+ "_receptor." + fileext
    if isfile(oldrec):
        rename(oldrec, newrec)
    else:
        if not isfile(newrec):
            logging.error("Receptor file not found for "+ folder) # TODO: List of warnings and errors
    logging.info("Docking experiment:"+ folder )
    y.ApplyMacro(filename="dock_runensemble_dummy.mcr",                 targets=newrec,                  remove="FromUnderscore",                  newextension="")

    # TODO: Continue from here: read the log file that is stored as <filename>.log
    
    logging.info("Exp "+folder+" Ended")
    
    
# yasara_docking("IL-36Ra(155)", "IL-36Ra(155)_1md6-a01.yob")
# yasara_docking("IL-1Ra(143)", "IL-1Ra(143).yob")


# In[6]:


"""if __name__ == '__main__':
    # assert len(sys.argv) == 2, 'Usage is "python organize_seqs <fastafile>"'
    seq_dict = fasta_to_seq(sys.argv[1])
    organize_sequences(seq_dict)
    for filename in seq_dict.keys():
        yasara_homology(filename)""";


# ### Selection cell 

# In[12]:


if __name__=="__main__":
    import yasara as y
    y.info.mode="txt"
    y.ApplyMacro(filename="/home/imhof_team/Public/mauricio/workflow/dummy.mcr",                 targets="Saved_run/IL-36Ra(155)/IL-36Ra(155)_1md6-a01_receptor.yob",                 remove="FromUnderscore",                 newextension="")


# In[18]:


if __name__=="__main__":
    macromodel = "dock_model.mcr"
    macro = "dock_runensemble_dummy.mcr"

    with open(macromodel) as f:
        macrotext = f.read()
    macrotext = "flexres = ''".join(macrotext.split("# $$FLEXRES$$"))
    macrotext = "Cell Auto,extension=17, Shape=Cuboid, Res Tyr 29".join(macrotext.split("# $$CELL$$"))
    print(a)
    with open(macro, "w") as f:
        f.write(macrotext)


# ### Define the flexible residues

# In[30]:


_aa_dict = {"G": "Gly", "P": "Pro", "A": "Ala", "V": "Val", "L": "Leu",             "I": "Ile", "M": "Met", "C": "Cys", "F": "Phe", "Y": "Tyr",             "W": "Trp", "H": "His", "K": "Lys", "R": "Arg", "Q": "Gln",             "N": "Asn", "E": "Glu", "D": "Asp", "S": "Ser", "T": "Thr"}

def flexchain(sequence, coordatom):
    coordpos = int(coordatom[1:])-1 # coordatom starts in 1, index start in 0
    seqstart = max(coordpos-4, 0)
    seqstart = min(seqstart, len(sequence)-9)
    result= ["Res"]
    for i, j in enumerate(sequence[seqstart:seqstart+9], seqstart+1):
        result+= [_aa_dict[j]]
        result+= [str(i)]
    return " ".join(result)
if__name__ == "__main__":
    flexchain("ARCDEFGHIQKLMN", "C1")

