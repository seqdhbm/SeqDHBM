#!/usr/bin/env python
# coding: utf-8

# In[3]:


OUTPUT_FILE = "output_yeast_structure.csv"
INPUT_FILE = '/home/imhof_team/Public/mauricio/workflow/yeast_type3/yeast_type3_safe2.fasta'

import logging
logging.basicConfig(level=logging.DEBUG, format='%(levelname)s - %(message)s')

#import experiments
import seqdhbm.SeqDHBM as SeqDHBM
import seqdhbm.fasta as fasta
import seqdhbm.pdbfiles as pdbfiles
from SeqDHBM import models
import sys, os
import shutil
import time

def workflow(jobfolder= "J0", fastafile=None, pdbs=[], pdbid="", rawseq="", mode="structure", job=None):
    """ Runs the workflow

    Params:
    - jobfolder: The relative folder where all the files from this run should be stored
    - fasta: One fasta file (no path)
    - pdbfiles: list of pdb filenames (no path)
    - pdbid: string with PDB codes (comma separated)
    - mode: 'wesa' to predict, 'structure' to skip"""
    # TODO: use jobnum to organize submissions

    # I think this variable lost utility
    seq_dict = {}

    seq_list = []
    # this has to come first, otherwise organize_sequences2 will run on unnecessary sequences
    if fastafile:
        try:
            """seq_dict = fasta.fasta_to_seq(fastafile)
            fasta.organize_sequences(seq_dict)"""
            seq_list = fasta.fasta_to_seq2(fastafile, jobfolder)
            fasta.organize_sequences2(seq_dict)

            #logging.debug("fasta:")
            #logging.debug(seq_dict.keys())
        except Exception as e:
            seq_list += [{"seq":"",
                         "name":"Your fasta File",
                         "folder": "",
                         "file": "",
                         "submited_as": "Fasta file",
                         "fail":True,
                         "warnings": ["Error reading the fasta file sent by the user."]}]
            logging.error("Error reading the fasta file sent by the user.")
            logging.error(type(e))
            logging.error(e)

    # Deal with the pdb ids (by downloading them)
    pdbid_dict = {} # stores the pdbfile path to each pdb id
    if pdbid:
        try:
            raise Exception("not prepared to work with the new variables")
            l = pdbfiles.text_to_list(pdbid)
            pdbid_dict = pdbfiles.get_pdb_files(l)
            #logging.debug("pdbid")
            #logging.debug(pdbid_dict)
        except Exception as e:
            logging.error("Error working with the pdb ids sent by the user.")
            logging.error(type(e))
            logging.error(e)


    if rawseq:
        try:
            seq_dict["Your input sequence"] = ("Your input sequence header", rawseq)
            folder = os.path.join(jobfolder, "MI")
            file = "MI1.fasta"
            name = "Your input sequence"
            os.makedirs(folder, exist_ok=True)
            with open(os.path.join(folder,file), "w") as f:
                f.write(">%s\n"%name)
                lines=fasta.break_fasta_sequence(rawseq)
                f.write('\n'.join(lines))
                f.write('\n')
            seq_list += [{"seq":rawseq,
                         "name":name,
                         "folder": folder,
                         "file": file,
                         "submited_as": "Manual input"}]
        except Exception as e:
            seq_list += [{"seq":"",
                         "name":"Your input sequence",
                         "folder": "",
                         "file": "",
                         "submited_as": "Manual input",
                         "fail":True,
                         "warnings": ["Error reading the manual input provided by the user."]}]
            logging.error("Error while saving the manual input to a file.")
            logging.error(type(e))
            logging.error(e)

    # Arrange the pdb files (user-submitted or downloaded)
    ## Put the structure file in folders - for docking and MD
    ## get the primary sequence for motif search
    for file in pdbs:
        try:
            raise Exception("not prepared to work with the new variables")
            try:
                folder = ".".join(file.split(".")[:-1]) # remove extension
            except:
                folder = file
            try:
                os.makedirs(folder, exist_ok=True)
                # TODO: organize the files into folders -- I should use an auto increment to create the folders? :D
                shutil.move(os.path.join(folder, file), folder+"/"+file)
                pdbid_dict[folder] = folder+"/"+file
            except:
                logging.error("Could not move file"+ file + "\nSkiping it!") # TODO: List of warnings and errors
        except Exception as e:
            logging.error("Error while working with the pdb files sent by the user.")
            logging.error(type(e))
            logging.error(e)


    # TODO: convert the pdbs into fasta
    for file in pdbid_dict.values():
        seq = pdbfiles.structure_to_fasta(file)
        if seq:
            seq_dict[file] = (">header", seq)
    # Run the motif check
    #logging.debug("*"*30)
    #logging.debug(seq_dict)

    #spacerinfo = {}

    hbm_result = {}
    cnt_progress = 0
    seq_dict = {}
    for item in seq_list:
        try:
            # TODO: SAVE HERE?
            seq_obj = None # model object
            if job:
                if (item["submited_as"] == "Manual input"):
                    _sub = models.Sequence.SUB_MANUAL_INPUT
                elif (item["submited_as"] == "Fasta file"):
                    _sub = models.Sequence.SUB_FASTA_FILE
                elif (item["submited_as"] == "PDB file"):
                    _sub = models.Sequence.SUB_PDB_FILE
                elif (item["submited_as"] == "PDB id"):
                    _sub = models.Sequence.SUB_PDB_ID
                _stat = models.Sequence.STATUS_FAILED if "fail" in item else models.Sequence.STATUS_QUEUED
                seq_obj = models.Sequence(jobnum=job, seqchain=item['seq'], submittedas=_sub,
                                          mode=models.Sequence.WESA_MODE if mode=="wesa" else models.Sequence.STRUCTURE_MODE,
                                          header = item["name"], status_hbm=_stat,
                                          fasta_file_location=os.path.join(item["folder"], item["file"]))
                seq_obj.save()
            # test if this is a record of an error
            if not "fail" in item:
                analysis = SeqDHBM.SpotCoordinationSite({">"+item["name"]: item["seq"]}, mode)
                analysedSeq = analysis[">" + item["name"]]
                # TODO: check fail in item
                if job and ("fail" in analysedSeq) and (analysedSeq["fail"]):
                    seq_obj.status_hbm = models.Sequence.STATUS_FAILED
                else:
                    seq_obj.status_hbm = models.Sequence.STATUS_PROCESSED

                # Output? Save in our user output?
                item["result"] = analysedSeq["result"]
                for coord, res in analysedSeq["result"].items():
                    res_obj = models.Result_HBM(
                        sequence = seq_obj,
                        coord_atom = coord,
                        ninemer = res["ninemer"],
                        net_charge = res['netcharge'],
                        disulfide_possible = bool(res['comment'].strip())
                    )
                    res_obj.save()

                item["warnings"] = analysedSeq["warnings"]
                seq_obj.warnings_hbm = "\n".join(analysedSeq["warnings"])
                seq_obj.save()
                cnt_progress+=1

        except Exception as e:
            logging.error("****")
            logging.error("Sequence name:%s \n %s"%(item["name"], item["file"]))
            logging.error(type(e))
            logging.error(e.args)
            logging.error(e)
            logging.error("****")
        finally:
            logging.debug("Got the results for %s"%item["name"])
            logging.debug("%d out of %d complete"%(cnt_progress, len(seq_dict)))

    return seq_list
