from django.test import TestCase
# from unittest import TestCase
from seqdhbm import SeqDHBM

# Create your tests here.


class TestSeqdhbm(TestCase):

    def test_sequence_validity_check(self):

        valid_sequences = [
            "MVLSEGEWQLVLHVWAKVEADVAGHGQDILIRLFKSHPETLEKFDRFKHLKTEAEMKASEDLKKAGVTVLTALGAILKKKGHHEAELKPL",
            "AALEPTDSGAPSAIVMFPVGEKPNPKGAAMKPVVFNHLIHEKKIADCETCHHTGDPVSCSTCHTVEGKAEGDYITLDRAM",
            "ETVLAEATVAPVSPMLAPYKVVIDALADKYEPSDFTHRRHLTSLMESIKDDKLAQAFHDKPEILCATCHHRSPLSLTPPK"
        ]
        invalid_sequences = [
            "MVLSEGEWQLVLHVWAKVEADVAGHGQDILIRLFKSHPETLEKFDRFKHLKT EAEMKASEDLKKAGVTVLTALGAILKKKGHHEAELKPL",
            "MVLSEGEWQLVLHVWAKVEADVAGHGQDILIRLFKSHPETLEKFDRFKHLKTEAEMKASEDLKKAGVTVLTALGAILKKKGHHEAELKP\nL",
            "AALEPTDSGAPSAIVBMF",
            "afhdkpeilcatchhr",
            "SHPETLEKFDR,F"
        ]
        for i in valid_sequences:
            self.assertFalse(SeqDHBM.sequence_validity_check(i))
        for i in invalid_sequences:
            self.assertTrue(SeqDHBM.sequence_validity_check(i))

    def test_build_ninemer_slice(self):
        valid_analysis = [
            ("MVLSEGEWQLVLHVWAKVEADVAGHGQDILIRLFKSHPETLEKFDRFKHLKTEAEMKASEDLKKAGVTVLTALGA", 21, "VEADVAGHG"),
            ("AALEPTDSGAPSAIVCMF", 13, "APSAIVCMF"),
            ("ETVLAEATVAPVSPMLAPYKVVI", 4, "ETVLAEATV"),
            ("EGKAEGDYITLDRAM", 2, "__EGKAEGD"),
            ("ILKKKGHHEAELKPL", 12, "EAELKPL__"),
            ("CAR", 0, "____CAR__")
        ]
        for in_seq, ind_, expected in valid_analysis:
            self.assertEqual(SeqDHBM.build_ninemer_slice(in_seq, ind_), expected)

    def test_check_basic_adjacent(self):
        valid_analysis = [
            ("QLVLHVWAK", "y"),
            ("RIIKYEFIL", "y"),
            ("QDDLCPLDG", "n"),
            ("LYLSCVLKD", "y"),
            ("LYLSHCVLD", "n"),
        ]
        for input_, expected in valid_analysis:
            self.assertEqual(SeqDHBM.check_basic_adjacent(input_), expected)

    def test_check_net_charge(self):
        valid_analysis = [
            ("KSFYHVSYG", 2),
            ("RIIKYEFIL", 1),
            ("QDDLCPLDG", -3),
            ("LYLSCVLKD", 0),
        ]
        for input_, expected in valid_analysis:
            self.assertEqual(SeqDHBM.get_net_charge(input_), expected)

    def test_precheck_coordinating_sites(self):
        valid_analysis = [
            {"seq": "XNEGDAAKGEKEFNKCKACHMIQAPDGTDIKGGKTGPNLYGVVGRKIASEEGFKYGEGILEVAEKNPDLTWTEANLIEYV",
             "site_count": 6, "C": 2, "H": 1, "Y": 3},
            {"seq": "HHEAELKPLAQSHATKHKIPIKYLEFISEAIIHVLHSRHPGDFGADAQGAMNKALELFRKDIAAKYKELGYQG",
             "site_count": 10, "C": 0, "H": 7, "Y": 3},
            {"seq": "MVLSEGEWQLVLHVWAKVEADVAGHGQDILIRLFKSHPETLEKFDRVKHLKTEAEMKASEDLKKHGVTVLTALGAILKKK",
             "site_count": 5, "C": 0, "H": 5, "Y": 0},
            {"seq": "XNEGDAAKGEKEFNKCKACMIQAPDGTDIKGGKTGPNLYGVVGRKIASEEGFKYGEGILEVAEKNPDLTWTEANLIEYV",
             "site_count": 5, "C": 2, "H": 0, "Y": 3},
        ]
        for expected in valid_analysis:
            actual_usrout = []
            site_cnt, c, h, y = SeqDHBM.precheck_coordinating_sites(expected["seq"],
                                                                    ">header",
                                                                    actual_usrout)
            self.assertEqual(site_cnt, expected["site_count"])
            self.assertEqual(c, expected["C"])
            self.assertEqual(h, expected["H"])
            self.assertEqual(y, expected["Y"])
            self.assertIn("WORKING ON THE VALID SEQUENCE: header", actual_usrout)

