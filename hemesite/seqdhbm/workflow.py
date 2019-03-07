# coding: utf-8

"""Organize the submission in files/database and call the analysis"""

import logging
import os

from SeqDHBM import models, tasks

from . import fasta, pdbfiles, SeqDHBM


def manage_fasta_files(fastafile: str, jobfolder: str):
    """
    Reads the sequences from the fasta file and save them in the appropriate folder.

    :param fastafile: File name
    :param jobfolder: Where the individual sequences will be saved in disk
    :return: a list with each sequence as a dictionary.
    """

    try:
        seq_list = fasta.fasta_to_seq2(fastafile, jobfolder)
        if seq_list:
            print(tasks.assync_organize_seq.delay(seq_list))
        else:
            raise Exception("The file doesn't contain fasta sequences.")
        return seq_list
    except Exception as e:
        logging.error("Error reading the fasta file sent by the user.")
        logging.error(type(e))
        logging.error(e)
        return [{"seq": "",
                 "name": "Your fasta File",
                 "folder": "",
                 "file": "",
                 "submitted_as": "Fasta File",
                 "fail": True,
                 "warnings": ["Error reading the fasta file sent by the user.",
                              e]
                 }]


def manage_pdb_ids():
    """l = pdbfiles.text_to_list(pdbid)
    pdbid_dict = pdbfiles.get_pdb_files(l)"""
    # TODO save the fasta file and store the sequence in the dictionary
    raise Exception("Not implemented")


def manage_raw_sequence(rawseq, jobfolder):
    folder = os.path.join(jobfolder, "MI")  # Manual Input
    file = "MI1.fasta"
    name = "Your input sequence"
    os.makedirs(folder, exist_ok=True)
    with open(os.path.join(folder, file), "w") as f:
        f.write(">%s\n" % name)
        lines = fasta.break_fasta_sequence(rawseq)
        f.write('\n'.join(lines))
        f.write('\n')
    return [{"seq": rawseq,
             "name": name,
             "folder": folder,
             "file": file,
             "submitted_as": "Manual Input"}]


def manage_pdb_files():
    """
    Arrange the pdb files (user-submitted or downloaded)
     - Put the structure file in folders - for docking and MD
     - get the primary sequence for motif search

    :return:
    """

    # TODO implement when needed
    """
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
            logging.error("Could not move file" + file + "\nSkipping it!")"""
    raise Exception("Not implemented")


def workflow(jobfolder: str = "J0",
             fastafile: str = None,
             pdbs: list = None,
             pdbid: str = "",
             rawseq: str = "",
             mode: str = "structure",
             job: models.Job = None) -> list:
    """
     Runs the workflow:
     1) Reads the inputs from the user from all possible sources (fasta file, fasta sequence (and also pdb files,
     pdb ids in the future).
     2) Save the fasta sequences read into files and to the database.
     3) For each sequence: run the analysis, save the results in the database.

    :param jobfolder: The relative folder where all the files from this run should be stored.
    :param fastafile: One fasta file (no path).
    :param pdbs: list of pdb filenames (no path).
    :param pdbid: string with PDB codes (comma separated).
    :param rawseq: an amino acid sequence: one letter, only 20 standard, no spaces, no linebreaks.
    :param mode: 'wesa' (to predict the structure) or 'structure' (when structure is assumed to be known).
    :param job: the job object, if exists.
    :return: A list of results in the form of a list of dictionaries. Each dictionary contains one of
    the submitted sequences. The fields in the dictionary are:
    # seq: str - The amino acid sequence,
    # name: str - The header,
    # folder: str - the path to the file,
    # file: str - the name of the file,
    # submitted_as: str - How it was submitted (Manual Input, Fasta File, PDB file, PDB id),
    # fail: bool - if the analysis failed,
    # warnings: list[str] - Warnings and errors
    # result: dict - The results as a dictionary, where the keys are the coordinating amino acids and the values are
    a dictionary with '9mer', 'charge', 'disulfide brige', ...
    # analysis: list[str] - A long text description of the steps of the analysis.
    """

    if pdbs is None:
        pdbs = []

    seq_list = []  # The sequences and their analysis results

    # this has to come first, otherwise organize_sequences2 will run on unnecessary sequences
    if fastafile:
        seq_list += manage_fasta_files(fastafile, jobfolder)

    # Deal with the pdb ids (by downloading them)
    pdbid_dict = {}  # stores the pdbfile path to each pdb id
    if pdbid:
        try:
            manage_pdb_ids()
        except Exception as e:
            logging.error("Error working with the pdb ids sent by the user.")
            logging.error(type(e))
            logging.error(e)

    if rawseq:
        try:
            seq_list += manage_raw_sequence(rawseq, jobfolder)
        except Exception as e:
            seq_list += [{"seq": rawseq,
                          "name": "Your input sequence",
                          "folder": "",
                          "file": "",
                          "submitted_as": "Manual Input",
                          "fail": True,
                          "warnings": ["Error reading the manual input provided by the user.",
                                       e]
                          }]
            logging.error("Error while saving the manual input to a file.")
            logging.error(type(e))
            logging.error(e)

    if pdbs:
        try:
            # TODO: Work with pdbs?
            seq_list += manage_pdb_files()
        except Exception as e:
            logging.error("Error while working with the pdb files sent by the user.")
            logging.error(type(e))
            logging.error(e)

    # TODO: convert the pdbs into fasta
    for file in pdbid_dict.values():
        seq = pdbfiles.structure_to_fasta(file)
        if seq:
            seq_list += [{}]  # Implement

    # Run the motif check
    for item in seq_list:
        seq_obj = None  # model object
        if job:
            _sub = None
            for sub_form_id, sub_form_desc in models.Sequence.SUBMISSION_FORMS:
                if item["submitted_as"] == sub_form_desc:
                    _sub = sub_form_id
            _stat = models.Sequence.STATUS_FAILED if "fail" in item else models.Sequence.STATUS_QUEUED
            seq_obj = models.Sequence(
                jobnum=job,
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

        # Test if there is no input error
        if "fail" not in item:
            # Then run the analysis
            analysis = SeqDHBM.spot_coordination_site(
                {">"+item["name"]: item["seq"]},
                seq_obj.id,
                mode
                )
            # and update with the results
            analysed_seq = analysis[">" + item["name"]]
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
    return seq_list
