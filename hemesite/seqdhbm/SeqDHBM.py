# coding: utf-8

import json
import logging
import os
import re
import time
from typing import Dict, List
import unittest

from django.conf import settings
from django.test import TestCase
from prettytable import PrettyTable

from seqdhbm import fasta


class Result(dict):
    def __init__(self, ninemer, netcharge, comment, kd):
        dict.__init__(
            self,
            ninemer=ninemer,
            netcharge= netcharge,
            comment=comment,
            kd=kd
        )


class Analysis:
    def __init__(
            self,
            mode: str,
            result: Dict[str, Result] = None,
            warnings: List[str] = None,
            analysis: str = '',
            fail: bool = False
    ):
        self.result = result if result else dict()
        self.warnings = warnings if warnings else list()
        self.mode = mode
        self.analysis = analysis
        self.fail = fail


def sequence_validity_check(current_sequence_list: str) -> int:
    """Takes a list of characters and checks if each element in the list is
    part of the 20 standard amino acids

    Return a flag value based on the number of errors in the input
    If all characters are valid then the flag should be 0
    """
    standard_aa_list = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L',
                        'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']  # List of the 20 standard amino acids
    flag = 0
    for aa in current_sequence_list:
        if aa not in standard_aa_list:
            flag += 1
    return flag


def mylog(lvl, buffer: list, line: str):
    """Stores messages to produce the full analysis and logs them, when necessary.
    Should save into a file and print for the standard output.

    :param lvl: The log level (debug, info, warning, ...)
    :param buffer: The list being used to save the analysis
    :param line: The message to be stored/logged
    """

    # FOR MORE OUTPUT IN DJANGO CONSOLE, UNCOMMENT THE LINE BELLOW
    # logging.log(lvl, line)
    buffer += [line]


def build_ninemer_slice(seq: str, ind: int) -> str:
    """Return the 9mer sequence corresponding to the given sequence and to the index of the coordinating amino acid.
    If the index points to the terminals of the sequence. add an appropriate number of amino acids
    add a filler character (_).

    :param seq: The complete sequence chain
    :param ind: The index of the coordinating a.a.
    :return: A 9mer.
    """

    # get the 4 aa before the coordinating. If at the beginning of the seq, add '_'
    ninemer = seq[max(ind-4, 0):ind].rjust(4, "_") + seq[ind:ind+5].ljust(5, "_")

    assert len(ninemer) == 9, "The ninemer could not be built properly"
    return ninemer


def check_basic_adjacent(ninemer_sequence: str) -> str:
    """Check if there are basic amino acids in the 9mer other than the coordinating AA

    :param ninemer_sequence: 9mer sequence
    :return: 'y' if a positively charged amino acid was found, otherwise 'n'
    """
    basic_aa_list = ['R', 'K', 'H']  # List of basic and positively charged amino acids
    for i in ninemer_sequence[:4] + ninemer_sequence[5:]:
        if i in basic_aa_list:
            return "y"
    return "n"


def get_net_charge(ninemer_sequence:str):
    """Compute the net charge on the motif

    :param ninemer_sequence: The 9mer
    :return: The charge value of the 9mer
    """
    basic_aa_list = ['R', 'K', 'H']  # List of basic and positively charged amino acids
    neg_aa_list = ['D', 'E']  # List of negatively charged amino acids
    net_charge = 0
    for i in range(len(ninemer_sequence)):
        if ninemer_sequence[i] in basic_aa_list:
            net_charge += 1
        if ninemer_sequence[i] in neg_aa_list:
            net_charge -= 1
    return net_charge


def ship_seq_to_wesa2(seq: str, seq_pk: int):
    """Simply ships the sequence to wesa using an standardized jobname.

    :param seq: The sequence chain.
    :param seq_pk: The identifier of the sequence
    """

    wesa_dir = os.path.join(settings.BASE_DIR, "wesa")
    inp_file = os.path.join(wesa_dir, f"seq_for_wesa{seq_pk}.txt")
    with open(inp_file, 'w') as f:  # Open a text file named seq_for_wesa.txt and write the sequence into it
        f.write(seq)
    jobname = f"SeqDHBM2WESA{seq_pk}"
    email = "seqdhbm@gmail.com"
    html = os.path.join(wesa_dir, jobname+'.html')
    wesa_tup = ['python2 ',
                os.path.join(wesa_dir, 'WESA-submit.py'),
                jobname,
                email,
                inp_file,
                '>',
                html
                ]  # This is stored as a tuple
    wesa_str = ' '.join(wesa_tup)  # Convert the list to string
    os.system(wesa_str)  # Send the job to the WESA server


def get_results_from_wesa(seq_pk):
    """Checks the WESA server if the results are ready.

    :param seq_pk: The identifier of the sequence
    :return: When ready, a dictionary with the position
    of each aa (keys) and the surface accessibility,
     Otherwise returns false
    """

    wesa_dir = os.path.join(settings.BASE_DIR, "wesa")
    jobname = f"SeqDHBM2WESA{seq_pk}"
    html = os.path.join(wesa_dir, jobname+'.html')
    download = os.path.join(wesa_dir, jobname+'.out')
    with open(html) as f:
        html_content = f.read().splitlines()
    html_content = "".join(html_content)
    # First look for some error in the html content (submission error)
    match = re.search("<h1>WESA Error</h1><p>(.*?)</p>", html_content)
    if match:  # Raise an assertion error if any error was found
        raise AssertionError(match.group(1))
    wesa_wget_tup = 'wget -O', download, '-F -i', html
    wesa_wget_str = ' '.join(wesa_wget_tup)  # Convert the tuple to string
    os.system(wesa_wget_str)  # Download the output as a .out file
    time.sleep(1)
    with open(download) as f:
        lines = f.read()
    if "Index of /wesa/out/" in lines:
        return False
    elif ("Job done" in lines) or ("Online version of the prediction is" in lines):
        pass
    else:
        logging.error(
            "The WESA Team might have changed the content of the results"
        )
        logging.error(f"Check the file SeqDHBM2WESA{seq_pk}.out" +
                      " and adapt the program.")
        return False
    # Split in lines and remove leading and trailing spaces
    lines = [x.strip() for x in lines.splitlines()]
    # find the lines containing surface accessibility information
    regex = re.compile("^([0-9]+)+ [A-Z] ")
    # boolean dictionary about surface accessibility info
    exposed_aa = dict()
    for line in lines:
        match = regex.search(line)
        if match:
            exposed_aa[int(match.group(1))] = (line[-7:-6] == '1')
    return exposed_aa


def disulfide_bond_check(ninemer, cys_cnt):
    """Checks if the ninemer contains a non-coordinating Cys a.a. and
    there is another Cys in the protein.

    :param ninemer: The 9mer being analysed
    :param cys_cnt: The number of Cys in the whole chain
    :return: The text 'Possible S-S Bond', if the disulfide bridge is possible.
    """
    return "Possible S-S Bond" if (cys_cnt > 1 and ("C" in ninemer[:4]+ninemer[5:])) else " "*17


def dummy_function_for_ajay(sequence: str) -> float:
    """"""
    if sequence:
        return 0.
    else:
        return 1.


def spot_coordination_site(
        fasta_dict: Dict['str', 'str'],
        seq_id: int,
        mode: str
) -> Dict[str, Analysis]:
    """Analyse the sequence for transient binding sites.
    1. Find Cys, His or Tyr
    2. Check if the adjacency is basic
    3. If the solvent accessibility has to be predicted, submit the sequence
    to the WESA service
    4. Compile the results

    :param fasta_dict: The sequences to be analysed in a dictionary formatted
    as {header: sequence}
    :param seq_id: The database id of the sequence
    :param mode: 'wesa' for solvent accessibility prediction. 'structure' to
    skip this check
    :return: The predicted coodinating sites.
    """

    # TODO refactor this long function

    filtered_dict = {}  # To store only valid sequences
    max_seq_length_for_wesa = 2000
    usrout = []

    output = {}
    out = {}
    """ the results of the check. Key is the fasta header and the values are another
         dictionary with:
         results:{coordinating_atom: {ninemer, netcharge, spacercheck}
         warnings: []
         analysis ""
         fail: true/false
         } """
    """ the results of the check. Key is the fasta header and the values are a
         named tuple with:
         results:{coordinating_atom: Results
         warnings: []
         analysis ""
         fail: true/false
         } """

    ###############################
    # Check for sequence validity #
    ###############################
    for header, sequence in fasta_dict.items():
        current_sequence = sequence

        if len(current_sequence) < 9:
            mylog(logging.INFO, usrout, "--------------------------")
            mylog(logging.INFO, usrout, "SEQUENCE VALIDITY CHECK:")
            mylog(logging.INFO, usrout, "--------------------------")
            mylog(logging.INFO, usrout, "The sequence is too short (minimum: 9)")
            mylog(logging.INFO, usrout, "-" * 80)
            output[header] = {"result": {},
                              "analysis": "\n".join(usrout),
                              "fail": True,
                              "mode": mode,
                              "warnings": ["The sequence is too short (minimum: 9 amino acids)"]}
            out[header] = Analysis(
                mode=mode,
                analysis="\n".join(usrout),
                fail=True,
                warnings=["The sequence is too short (minimum: 9 amino acids)"]
            )
        elif sequence_validity_check(current_sequence) == 0:  # Pass the list to the SequenceValidityCheck function
            filtered_dict[header] = sequence  # Add only the valid sequences to the new dictionary
        else:
            mylog(logging.INFO, usrout, "--------------------------")
            mylog(logging.INFO, usrout, "SEQUENCE VALIDITY CHECK:")
            mylog(logging.INFO, usrout, "--------------------------")
            mylog(logging.INFO, usrout, "The following sequence(s) contains invalid characters other " +
                                        "than the 20 standard amino acids")
            mylog(logging.INFO, usrout, "-" * 80)
            output[header] = {"result": {},
                              "analysis": "\n".join(usrout),
                              "fail": True,
                              "mode": mode,
                              "warnings": ["The sequence contains invalid characters " +
                                           "other than the 20 standard amino acids"]}
            out[header] = Analysis(
                mode=mode,
                analysis="\n".join(usrout),
                fail=True,
                warnings=["The sequence contains invalid characters " +
                          "other than the 20 standard amino acids"]
            )

    # Add information about the sequences with error to the output analysis.
    if len(fasta_dict) != len(filtered_dict):
        for header, sequence in fasta_dict.items():
            if header not in filtered_dict:
                mylog(logging.INFO, usrout, "%s\n%s" % (header[1:], "\n".join(fasta.break_fasta_sequence(sequence))))

    if len(filtered_dict) == 0:  # if there are no valid sequences...
        return out  # output

    #################################
    # Work with the valid sequences #
    #################################
    seq_with_coord = 0
    for header, current_sequence in filtered_dict.items():  # Iterate through the dictionary
        initial_ninemer_dict = {}  # An empty dictionary that will contain the residue number
        # and the associated 9 mer as key value pair eg. {H150: 'XXXXHXXXX'}
        basic_adj_list = []  # Initialize an empty list to collect the flags when the check for basic adjacent
        # amino acids is done. Values will be either yes or no.
        dict_for_net_charge_calc = {}  # 9mer motifs that PASS the basic adjacency test are populated here
        pass_charge_out_dict = {}
        # for comparison with the WESA output
        output[header] = {
            "warnings": [],
            "result": {},
            "analysis": "",
            "mode": mode,
            "fail": False}  # initialize the informative dictionary
        out[header] = Analysis(
            mode=mode,
        )

        ###############################
        # Check for coodinating sites #
        ###############################
        coord_site_count, cys_count, his_count, tyr_count = precheck_coordinating_sites(current_sequence,
                                                                                        header,
                                                                                        usrout)

        if coord_site_count:
            seq_with_coord += 1
            mylog(logging.INFO, usrout, "NOTE : Heme coordination check PASS!")
            mylog(logging.INFO, usrout, f"Length of sequence: {len(current_sequence)}")
            mylog(logging.INFO, usrout, f"Total number of potential coordination " +
                                        f"sites found: {coord_site_count}")
            mylog(logging.INFO, usrout, "Number of potential CYS based sites: "+str(cys_count))
            mylog(logging.INFO, usrout, "Number of potential HIS based sites: "+str(his_count))
            mylog(logging.INFO, usrout, "Number of potential TYR based sites: "+str(tyr_count)+"\n")

            coord_site_index_list = ([pos for pos, char in enumerate(current_sequence)
                                      if char in('H', 'C', 'Y')])  # List with the indices of the coordination sites
            for i, pos in enumerate(coord_site_index_list):  # Iterate through the list containing the
                # indices of the coordination sites
                pos_on_sequence = pos + 1  # Position on the sequence = index + 1 since Python index starts at 0
                coord_aa = current_sequence[pos]  # Variable stores the value on the index
                # which is the actual amino acid residue
                formatted_coord_site = coord_aa+str(pos_on_sequence)  # Variable stores the amino acid
                # and index in the format eg. H151
                ninemer = build_ninemer_slice(current_sequence, pos)  # Call the BuildNinemerSlice function
                initial_ninemer_dict[formatted_coord_site] = ninemer  # Populate the initial_NinemerDict {}
                # in the format eg. {H12: XXXXHXXXX}

            if len(coord_site_index_list) == len(initial_ninemer_dict):  # Final record of the initial_NinemerDict{}
                mylog(logging.INFO, usrout, "NOTE : Successfully built 9mer motifs for all the " +
                      str(coord_site_count) + " potential coordination sites !")
                mylog(logging.INFO, usrout, "-" * 50)

            #######################
            # Check for CP motifs #
            #######################
            cp_motif = {}  # This dictionary will skip basic adjacency
            for aa, ninemer in initial_ninemer_dict.items():
                # check if the coordinating aa is Cysteine
                if aa[0] == 'C':
                    pos = int(aa[1:]) - 1  # the position on string
                    dict_for_net_charge_calc[aa] = ninemer
                    # If it is a CP motif, also store in a dictionary
                    if pos < len(current_sequence)-1 and current_sequence[pos+1] == 'P':
                        cp_motif[aa] = ninemer

            ###############################
            # Check for basic amino acids #
            ###############################
            mylog(logging.INFO, usrout, "ADJACENT BASIC AMINO ACID CHECK:")
            mylog(logging.INFO, usrout, "-" * 35)
            mylog(logging.INFO, usrout, "NOTE : Screening all potential motifs for adjacent basic amino acids")

            for key, ninemer in initial_ninemer_dict.items():
                basic_adj_list.append(check_basic_adjacent(ninemer))  # this list contains values "y" or "n"
                if check_basic_adjacent(ninemer) == "y":
                    dict_for_net_charge_calc[key] = ninemer

            count_basic_adj_ninemer = (basic_adj_list.count("y"))

            if count_basic_adj_ninemer > 0:
                mylog(logging.INFO, usrout, "NOTE : %d out of the %d 9mers PASS the adjacent basic amino acids check" %
                      (count_basic_adj_ninemer, coord_site_count))
                if len(dict_for_net_charge_calc)-count_basic_adj_ninemer:
                    mylog(logging.INFO,
                          usrout,
                          "NOTE : An additional %d " % (len(dict_for_net_charge_calc) - count_basic_adj_ninemer) +
                          "motif coordinated by Cystein will be checked for charge")
                mylog(logging.INFO, usrout, "NOTE : These %d " % len(dict_for_net_charge_calc) +
                                            "9mers will be used in the next step to check for positive net charge")
            elif dict_for_net_charge_calc:
                mylog(logging.INFO, usrout, "NOTE : None of the %d " % coord_site_count +
                                            "9mers passed the adjacent basic amino acids check")
                mylog(logging.INFO, usrout, "NOTE : %d 9mers coordinated by Cysteine " % len(dict_for_net_charge_calc) +
                                            "will be used in the next step to check for net charge")
            else:
                mylog(logging.INFO, usrout, "-" * 50)
                mylog(logging.INFO, usrout, "NOTE : None of the 9mers have valid motifs")
                mylog(logging.INFO, usrout, "NOTE : The likelyhood of any part of this sequence binding" +
                                            " or coordinating heme is very low !!!")
                mylog(logging.INFO, usrout, "NOTE : Checking the next valid sequence")
                output[header]["warnings"] += ["None of the 9mers have valid motifs.",
                                               "The likelihood of any part of this sequence binding" +
                                               " or coordinating heme is very low"]
                out[header].warnings.extend([
                    "None of the 9mers have valid motifs.",
                    "The likelihood of any part of this sequence binding"
                    + " or coordinating heme is very low"
                ])
                output[header]["fail"] = False
                out[header].fail = False

                mylog(logging.INFO, usrout, "-" * 50)
                mylog(logging.INFO, usrout, "-" * 50)
                continue

            # NET CHARGE CHECK
            pass_charge_dict = net_charge_check(cp_motif, dict_for_net_charge_calc, pass_charge_out_dict, usrout)

            if len(pass_charge_dict) > 1:
                mylog(logging.INFO, usrout, "NOTE : Additional coordination site check: PASS")
                mylog(logging.INFO, usrout, "NOTE : Proceeding to check spacer length between coordination sites")
            elif len(pass_charge_dict) == 1:
                mylog(logging.INFO, usrout, "NOTE : No additional coordination sites found !")
                mylog(logging.INFO, usrout, "NOTE : Spacer check will be skipped !")
            else:
                mylog(logging.INFO, usrout, "NOTE : No coordination sites passed the net charge check !")
                mylog(logging.INFO,
                      usrout,
                      "NOTE : It is very likely that this sequence does not bind/coordinate heme !")
                output[header]["warnings"] += [
                    "No coordination sites passed the net charge check.",
                    "It is very likely that this sequence does not bind/coordinate heme"
                ]
                out[header].warnings.extend([
                    "No coordination sites passed the net charge check.",
                    "It is very likely that this sequence "
                    + "does not bind/coordinate heme"
                ])
                continue

            # (Re-add when necessary)
            # spacer_dict = spacer_check(
            #    pass_charge_dict, pass_charge_index_list, spacer_dict
            # )

            # WESA supports only sequences with length up to 2000
            if (mode == "wesa") and (len(current_sequence) > max_seq_length_for_wesa):
                inform_wesa_incompatibility(len(current_sequence), header, max_seq_length_for_wesa, out, usrout)

            if (mode == "structure") or (len(current_sequence) > max_seq_length_for_wesa):
                # structure
                save_results_structure_mode(cys_count, header, out, pass_charge_dict, pass_charge_out_dict, usrout)
            else:
                # WESA
                mylog(logging.INFO, usrout, "Sending sequence to the WESA server for solvent accessibility prediction")
                ship_seq_to_wesa2(current_sequence, seq_id)
                output[header]["result"] = {}
                # out[header]. result ja ta inicializado
                for coord in sorted(pass_charge_dict.keys(), key=lambda x: int(x[1:])):
                    out[header].result[coord] = Result(
                        ninemer=pass_charge_dict[coord],
                        netcharge=pass_charge_out_dict[coord],
                        comment=disulfide_bond_check(
                            pass_charge_dict[coord],
                            cys_count
                        ),
                        kd=dummy_function_for_ajay(pass_charge_dict[coord])
                    )
                    output[header]["result"][coord] = {"ninemer": pass_charge_dict[coord],
                                                       "netcharge": pass_charge_out_dict[coord],
                                                       "comment": disulfide_bond_check(pass_charge_dict[coord],
                                                                                       cys_count)}
            output[header]["analysis"] = "\n".join(usrout)
            out[header].analysis = "\n".join(usrout)

        else:  # When there are no coordinating sites in the sequence
            no_coord_site_warning(header, out, usrout)
            return out
    return out


def precheck_coordinating_sites(current_sequence: str, header: str, usrout: list):
    """counts C, H and Y. Add that information to the full analysis"""
    mylog(logging.INFO, usrout, "-" * 80)
    mylog(logging.INFO, usrout, "WORKING ON THE VALID SEQUENCE: " + header[1:])
    mylog(logging.INFO, usrout, "-" * 80)
    mylog(logging.INFO, usrout, "\n".join(fasta.break_fasta_sequence(current_sequence)))
    mylog(logging.INFO, usrout, "-" * 30)
    mylog(logging.INFO, usrout, "HEME COORDINATION SITE CHECK:")
    mylog(logging.INFO, usrout, "-" * 30)
    cys_count = current_sequence.count('C')
    his_count = current_sequence.count('H')
    tyr_count = current_sequence.count('Y')
    coordination_site_count = cys_count + his_count + tyr_count
    return coordination_site_count, cys_count, his_count, tyr_count


def net_charge_check(cp_motif, dict_for_net_charge_calc, pass_charge_out_dict, usrout):
    mylog(logging.INFO, usrout, "-" * 35)
    mylog(logging.INFO, usrout, "NET CHARGE CHECK:")
    mylog(logging.INFO, usrout, "-" * 35)
    mylog(logging.INFO, usrout, "NOTE : Checking net charge on the individual 9mer motifs")
    mylog(logging.INFO, usrout, "NOTE : An individual motif PASSES this check:")
    mylog(logging.INFO, usrout, "       if it's net charge is positive")
    mylog(logging.INFO, usrout, "       if it's neutral, but it is a CYS BASED motif")
    mylog(logging.INFO, usrout, "       It's a CP coordinated motif")
    pass_charge_dict = {}  # contains 9 mers that have passed the net positive charge check
    for d in dict_for_net_charge_calc:
        charge = get_net_charge(dict_for_net_charge_calc[d])
        if charge > 0:
            pass_charge_dict[d] = dict_for_net_charge_calc[d]
            pass_charge_out_dict[d] = "+" + str(charge)
        elif charge == 0 and d[0] == "C":
            pass_charge_dict[d] = dict_for_net_charge_calc[d]
            pass_charge_out_dict[d] = str(charge) + "(CYS motif)"
        elif d in cp_motif:
            pass_charge_dict[d] = dict_for_net_charge_calc[d]
            pass_charge_out_dict[d] = str(charge) + "(CP motif)"
    mylog(
        logging.INFO,
        usrout,
        (f"NOTE : {len(pass_charge_dict)} out of " +
         f"{len(dict_for_net_charge_calc)} PASS the net charge check!")
    )
    # Additional coordination site check
    mylog(logging.INFO, usrout, "-" * 50)
    mylog(logging.INFO, usrout, "ADDITIONAL COORDINATION SITE CHECK")
    mylog(logging.INFO, usrout, "-" * 50)
    mylog(
        logging.INFO,
        usrout,
        "NOTE : Checking the presence of additional coordination sites"
    )
    return pass_charge_dict


def inform_wesa_incompatibility(current_sequence_length, header, max_seq_length_for_wesa, output, usrout):
    mylog(logging.WARNING, usrout, "-" * 50)
    mylog(logging.WARNING, usrout, "Warning:")
    msg1 = f"Structure prediction by the WESA server" + \
           f" supports sequences with length up to" + \
           f" {max_seq_length_for_wesa} amino acids"
    msg2 = f"This sequence has {current_sequence_length} a.a. " + \
           f"and is not supported."
    for msg in [msg1, msg2]:
        # output[header]["warnings"] += [msg]
        output[header].warnings.append(msg)
        mylog(logging.WARNING, usrout, msg)
    output[header]["mode"] = "structure"


def save_results_structure_mode(
        cys_count,
        header,
        output,
        pass_charge_dict,
        pass_charge_out_dict,
        usrout
):
    # All checks done, preparing tabular summary
    """Prepare output lists for PrettyTable"""
    mylog(logging.INFO, usrout, "*" * 80)
    mylog(logging.INFO, usrout, "TABULAR SUMMARY")
    mylog(logging.INFO, usrout, "NOTE : Please use the available structure "
                                "information to only consider those motifs that are \"surface exposed\"")
    table = PrettyTable([
        "S.no",
        "Coord. residue",
        "9mer motif",
        "Net charge",
        "Comment",
        "Kd or strength"
    ])
    sr_no = 0
    # output[header]["result"] = {}
    for coord in sorted(pass_charge_dict.keys(), key=lambda x: int(x[1:])):
        bond_check = disulfide_bond_check(
            pass_charge_dict[coord],
            cys_count
        )
        r = Result(
            ninemer=pass_charge_dict[coord],
            netcharge=pass_charge_out_dict[coord],
            comment=bond_check,
            kd=dummy_function_for_ajay(pass_charge_dict[coord])
        )
        output[header].result[coord] = r
        sr_no += 1
        table.add_row([
            sr_no,
            coord,
            pass_charge_dict[coord],
            pass_charge_out_dict[coord],
            bond_check,
            ""
        ])
    mylog(logging.INFO, usrout, str(table))


def no_coord_site_warning(header, output, usrout):
    mylog(logging.INFO, usrout, "-" * 80)
    msg1 = "NOTE : The sequence does not have potential heme coordination sites(C, H, Y) !!"
    msg2 = "It is very likely that this sequence does not bind/coordinate heme"
    for msg in [msg1, msg2]:
        # output[header]["warnings"] += [msg]
        output[header].warnings.append(msg)
        mylog(logging.INFO, usrout, msg)


def spacer_check(pass_charge_dict, pass_charge_index_list):
    """Checks for spacer distance.

    :param pass_charge_dict:
    :param pass_charge_index_list:
    :return: Dict containing "y" or "n" values to check if two coordination sites are
    separated by a spacer of at least 2
    """

    spacer_dict = {}
    for head in pass_charge_dict:  # Iterate through this dictionary
        # append to the pass_charge_index_list only the number from the
        # dictionary header eg. in C9 extract 9 and converts it to 'int'
        pass_charge_index_list.append(int(head[1:]))
    pass_charge_index_list.sort()
    if len(pass_charge_index_list) > 1:
        for i in range(len(pass_charge_index_list)):
            if i == 0:
                if abs((pass_charge_index_list[i] - pass_charge_index_list[i + 1])) > 2:
                    spacer_dict[pass_charge_index_list[i]] = "YES"
                else:
                    spacer_dict[pass_charge_index_list[i]] = "NO"
            elif i == len(pass_charge_index_list) - 1:
                if (abs(pass_charge_index_list[i] - pass_charge_index_list[i - 1])) > 2:
                    spacer_dict[pass_charge_index_list[i]] = "YES"
                else:
                    spacer_dict[pass_charge_index_list[i]] = "NO"
            else:
                if ((abs(pass_charge_index_list[i] - pass_charge_index_list[i + 1])) > 2 and
                        (abs(pass_charge_index_list[i] - pass_charge_index_list[i - 1])) > 2):
                    spacer_dict[pass_charge_index_list[i]] = "YES"
                else:
                    spacer_dict[pass_charge_index_list[i]] = "NO"
    else:
        spacer_dict[int(pass_charge_index_list[0])] = "NO"  # position of the only valid aa in the chain is NO
    return spacer_dict
