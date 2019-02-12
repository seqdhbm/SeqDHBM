#!/usr/bin/env python
# coding: utf-8

# In[3]:


from . import fasta, pdbfiles, SeqDHBM
import logging
import os
from SeqDHBM import models, tasks
import shutil
import sys
import time


def workflow(jobfolder= "J0", fastafile=None, pdbs=None, pdbid="", rawseq="", mode="structure", job=None):
    """ Runs the workflow

    Params:
    - jobfolder: The relative folder where all the files from this run should be stored
    - fasta: One fasta file (no path)
    - pdbfiles: list of pdb filenames (no path)
    - pdbid: string with PDB codes (comma separated)
    - mode: 'wesa' to predict, 'structure' to skip"""

    # I think this variable lost utility
    seq_dict = {}

    if pdbs is None:
        pdbs = []
    seq_list = []
    # this has to come first, otherwise organize_sequences2 will run on unnecessary sequences
    if fastafile:
        try:
            seq_list = fasta.fasta_to_seq2(fastafile, jobfolder)
            print(tasks.assync_organize_seq.delay(seq_list))
            # fasta.organize_sequences2(seq_list)

        except Exception as e:
            seq_list += [{"seq": "",
                          "name": "Your fasta File",
                          "folder": "",
                          "file": "",
                          "submitted_as": "Fasta file",
                          "fail": True,
                          "warnings": ["Error reading the fasta file sent by the user."]}]
            logging.error("Error reading the fasta file sent by the user.")
            logging.error(type(e))
            logging.error(e)

    # Deal with the pdb ids (by downloading them)
    pdbid_dict = {}  # stores the pdbfile path to each pdb id
    if pdbid:
        try:
            """l = pdbfiles.text_to_list(pdbid)
                        pdbid_dict = pdbfiles.get_pdb_files(l)
                        #logging.debug("pdbid")
                        #logging.debug(pdbid_dict)"""
            raise Exception("not prepared to work with the new variables")
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
            with open(os.path.join(folder, file), "w") as f:
                f.write(">%s\n" % name)
                lines=fasta.break_fasta_sequence(rawseq)
                f.write('\n'.join(lines))
                f.write('\n')
            seq_list += [{"seq": rawseq,
                          "name": name,
                          "folder": folder,
                          "file": file,
                          "submitted_as": "Manual input"}]
        except Exception as e:
            seq_list += [{"seq": "",
                          "name": "Your input sequence",
                          "folder": "",
                          "file": "",
                          "submitted_as": "Manual input",
                          "fail": True,
                          "warnings": ["Error reading the manual input provided by the user."]}]
            logging.error("Error while saving the manual input to a file.")
            logging.error(type(e))
            logging.error(e)

    # Arrange the pdb files (user-submitted or downloaded)
    # - Put the structure file in folders - for docking and MD
    # - get the primary sequence for motif search
    for file in pdbs:
        try:
            raise Exception("not prepared to work with the new variables")
            try:
                folder = ".".join(file.split(".")[:-1])  # remove extension
            except:
                folder = file
            try:
                os.makedirs(folder, exist_ok=True)
                # TODO: organize the files into folders -- I should use an auto increment to create the folders? :D
                shutil.move(os.path.join(folder, file), folder+"/"+file)
                pdbid_dict[folder] = folder+"/"+file
            except:
                logging.error("Could not move file" + file + "\nSkipping it!")  # TODO: List of warnings and errors
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
    # logging.debug("*"*30)
    # logging.debug(seq_dict)

    # spacerinfo = {}

    hbm_result = {}
    cnt_progress = 0
    seq_dict = {}
    for item in seq_list:
        # try:
            # TODO: SAVE HERE?
            seq_obj = None  # model object
            if job:
                if item["submitted_as"] == "Manual input":
                    _sub = models.Sequence.SUB_MANUAL_INPUT
                elif item["submitted_as"] == "Fasta file":
                    _sub = models.Sequence.SUB_FASTA_FILE
                elif item["submitted_as"] == "PDB file":
                    _sub = models.Sequence.SUB_PDB_FILE
                elif item["submitted_as"] == "PDB id":
                    _sub = models.Sequence.SUB_PDB_ID
                _stat = models.Sequence.STATUS_FAILED if "fail" in item else models.Sequence.STATUS_QUEUED
                seq_obj = models.Sequence(jobnum=job,
                                          seqchain=item['seq'],
                                          submittedas=_sub,
                                          mode=models.Sequence.WESA_MODE if mode == "wesa" else models.Sequence.STRUCTURE_MODE,
                                          header=item["name"],
                                          status_hbm=_stat,
                                          fasta_file_location=os.path.join(
                                            item["folder"],
                                            item["file"])
                                          )
                seq_obj.save()

            # Test if there is an input error
            if "fail" not in item:
                analysis = SeqDHBM.SpotCoordinationSite(
                    {">"+item["name"]: item["seq"]},
                    seq_obj.id,
                    mode
                    )
                analysed_seq = analysis[">" + item["name"]]
                # TODO: SET TO 'RECEIVED' NOW AND TO 'PROCESSED' AFTER THE RESULTS ARE SAVED
                if ("fail" in analysed_seq) and (analysed_seq["fail"]):
                    seq_obj.status_hbm = models.Sequence.STATUS_FAILED
                elif analysed_seq["mode"] != mode:
                    seq_obj.mode = (models.Sequence.STRUCTURE_MODE
                                    if analysed_seq["mode"] == "structure"
                                    else models.Sequence.WESA_MODE)
                seq_obj.partial_hbm_analysis = analysed_seq["analysis"]

                item["result"] = analysed_seq["result"]
                item["analysis"] = analysed_seq["analysis"]
                item["warnings"] = analysed_seq["warnings"]
                seq_obj.warnings_hbm = "\n".join(analysed_seq["warnings"])
                seq_obj.save()
                tasks.assync_save_results.delay(seq_obj.id, analysed_seq["result"])
                cnt_progress += 1

    """    except Exception as e:
            logging.error("****")
            logging.error("Sequence name:%s \n %s"%(item["name"], item["file"]))
            logging.error(type(e))
            logging.error(e.args)
            logging.error(e)
            logging.error("****")
        finally:
            logging.debug("Got the results for %s"%item["name"])
            logging.debug("%d out of %d complete"%(cnt_progress, len(seq_dict)))"""
    return seq_list
