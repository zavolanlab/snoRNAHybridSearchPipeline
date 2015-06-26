import unittest
from Bio.Seq import Seq
from snoRNA import (snoRNA,
                    CD_snoRNA,
                    IncompatibleStrandAndCoordsException,
                    WrongCDBoxOrderExeption,
                    WrongCDBoxPlacementException,
                    ToManyBoxesException,
                    NoBoxException)


class TestsnoRNA(unittest.TestCase):

    """Test snoRNA class"""

    def setUp(self):
        self.snor = snoRNA(snor_id="test_snoRNA",
                      organism="hg",
                      chrom="chr1",
                      start=100,
                      end=170,
                      strand="+",
                      sequence="A"*70,
                      snor_type="CD",
                           modified_sites="18S:U55,28S:A44")

    def test_creation_of_snoRNA(self):
        snor = snoRNA(snor_id="test_snoRNA",
                      organism="hg",
                      chrom="chr1",
                      start=100,
                      end=170,
                      strand="+",
                      sequence="A"*70,
                      snor_type="CD")
        #
        # test arguments
        #
        self.assertEqual(snor.snor_id, "test_snoRNA")
        self.assertEqual(snor.organism, "hg")
        self.assertEqual(snor.position.chrom, "chr1")
        self.assertEqual(snor.position.start, 100)
        self.assertEqual(snor.position.end, 170)
        self.assertEqual(snor.position.strand, "+")
        self.assertEqual(str(snor.sequence), "A"*70)
        self.assertEqual(snor.snor_type, "CD")
        #
        # test keyword arguments (all should be None)
        # this also test if the snoRNA has these attributes
        #
        self.assertEqual(snor.alias, None)
        self.assertEqual(snor.gene_name, None)
        self.assertEqual(snor.accession, None)
        self.assertEqual(snor.modified_sites, None)
        self.assertEqual(snor.host_gene, None)
        self.assertEqual(snor.host_id, None)
        self.assertEqual(snor.organization, None)
        self.assertEqual(snor.note, None)



    # I decided to delete this error message
    # def test_coordinates_and_sequence_length(self):
    #     with self.assertRaises(IncompatibleStrandAndCoordsException):
    #         snoRNA(snor_id="test_snoRNA",
    #                organism="hg",
    #                chrom="chr1",
    #                start=100,
    #                end=170,
    #                strand="+",
    #                sequence="A"*60,
    #                snor_type="CD")
    def test_modification_positions(self):
        self.assertEqual(self.snor.modified_sites["18S"][0][0], 55)
        self.assertEqual(self.snor.modified_sites["28S"][0][0], 44)

    def test_get_bed_string(self):
        bed = "chr1\t100\t170\tsnoRNA:test_snoRNA\tNone\t+\n"
        self.assertEqual(self.snor.get_bed_string(), bed)

    def test_get_fasta_string(self):
        fasta = ">test_snoRNA\n" + "A"*70 + "\n"
        self.assertEqual(self.snor.get_fasta_string(), fasta)


class TestCD_snoRNA(unittest.TestCase):

    """Test CD_snoRNA class"""

    def setUp(self):
        self.snor = CD_snoRNA(snor_id="test_snoRNA",
                   organism="hg",
                   chrom="chr1",
                   start=100,
                   end=164,
                   strand="+",
                   sequence="C"*5 + "T"*20 + "C"*10 + "A"*20 + "CATG" + "T"*5,
                   snor_type="CD",
                   d_boxes="26..29,56..59",
                   c_boxes="49..55,5..13")

    def test_for_non_dbox(self):
        with self.assertRaises(NoBoxException):
            cdsnor = CD_snoRNA(snor_id="test_snoRNA",
                               organism="hg",
                               chrom="chr1",
                               start=100,
                               end=170,
                               strand="+",
                               sequence="A"*70,
                               snor_type="CD",
                               d_boxes=None)

    def test_for_one_dbox_error(self):
        with self.assertRaises(WrongCDBoxPlacementException):
            cdsnor = CD_snoRNA(snor_id="test_snoRNA",
                               organism="hg",
                               chrom="chr1",
                               start=100,
                               end=170,
                               strand="+",
                               sequence="A"*70,
                               snor_type="CD",
                               d_boxes="20..23")

    def test_for_one_dbox(self):
        cdsnor = CD_snoRNA(snor_id="test_snoRNA",
                           organism="hg",
                           chrom="chr1",
                           start=100,
                           end=170,
                           strand="+",
                           sequence="A"*70,
                           snor_type="CD",
                           d_boxes="60..63")
        self.assertEqual(cdsnor.d_box, (60, 63))
        self.assertEqual(cdsnor.dprime_box, None)

    def test_for_two_dbox(self):
        cdsnor = CD_snoRNA(snor_id="test_snoRNA",
                           organism="hg",
                           chrom="chr1",
                           start=100,
                           end=170,
                           strand="+",
                           sequence="A"*70,
                           snor_type="CD",
                           d_boxes="20..23, 60..63")
        self.assertEqual(cdsnor.d_box, (60, 63))
        self.assertEqual(cdsnor.dprime_box, (20, 23))

    def test_for_two_dbox_error(self):
        with self.assertRaises(WrongCDBoxOrderExeption):
            CD_snoRNA(snor_id="test_snoRNA",
                               organism="hg",
                               chrom="chr1",
                               start=100,
                               end=170,
                               strand="+",
                               sequence="A"*70,
                               snor_type="CD",
                               d_boxes="60..63, 20..23")
        raised = False
        try:
            cdsnor = CD_snoRNA(snor_id="test_snoRNA",
                               organism="hg",
                               chrom="chr1",
                               start=100,
                               end=170,
                               strand="+",
                               sequence="A"*70,
                               snor_type="CD",
                               d_boxes="60..63, 20..23",
                               switch_boxes=True)
        except WrongCDBoxOrderExeption:
          raised = True
        self.assertFalse(raised, "WrongCDBoxOrderExeption raised with swith_boxes=True")
        self.assertEqual(cdsnor.d_box, (60, 63))
        self.assertEqual(cdsnor.dprime_box, (20, 23))

    def test_for_many_dbox_error(self):
        with self.assertRaises(ToManyBoxesException):
            cdsnor = CD_snoRNA(snor_id="test_snoRNA",
                               organism="hg",
                               chrom="chr1",
                               start=100,
                               end=170,
                               strand="+",
                               sequence="A"*70,
                               snor_type="CD",
                               d_boxes="60..63, 20..23, 30..33")

    def test_for_one_cbox(self):
        cdsnor = CD_snoRNA(snor_id="test_snoRNA",
                           organism="hg",
                           chrom="chr1",
                           start=100,
                           end=170,
                           strand="+",
                           sequence="A"*70,
                           snor_type="CD",
                           d_boxes="60..63",
                           c_boxes="30..36")
        self.assertEqual(cdsnor.c_box, (30, 36))
        self.assertEqual(cdsnor.cprime_box, None)

    def test_for_two_cbox(self):
        cdsnor = CD_snoRNA(snor_id="test_snoRNA",
                           organism="hg",
                           chrom="chr1",
                           start=100,
                           end=170,
                           strand="+",
                           sequence="A"*70,
                           snor_type="CD",
                           d_boxes="60..63",
                           c_boxes="50..56,30..36")
        self.assertEqual(cdsnor.c_box, (30, 36))
        self.assertEqual(cdsnor.cprime_box, (50, 56))

    def test_for_two_cbox_error(self):
        with self.assertRaises(WrongCDBoxOrderExeption):
            CD_snoRNA(snor_id="test_snoRNA",
                               organism="hg",
                               chrom="chr1",
                               start=100,
                               end=170,
                               strand="+",
                               sequence="A"*70,
                               snor_type="CD",
                               d_boxes="60..63",
                               c_boxes="30..36, 50..56")
        raised = False
        try:
            cdsnor = CD_snoRNA(snor_id="test_snoRNA",
                               organism="hg",
                               chrom="chr1",
                               start=100,
                               end=170,
                               strand="+",
                               sequence="A"*70,
                               snor_type="CD",
                               d_boxes="60..63",
                               c_boxes="30..36, 50..56",
                               switch_boxes=True)
        except WrongCDBoxOrderExeption:
          raised = True
        self.assertFalse(raised, "WrongCDBoxOrderExeption raised with swith_boxes=True")
        self.assertEqual(cdsnor.c_box, (30, 36))
        self.assertEqual(cdsnor.cprime_box, (50, 56))

    def test_for_many_cbox_error(self):
        with self.assertRaises(ToManyBoxesException):
            cdsnor = CD_snoRNA(snor_id="test_snoRNA",
                               organism="hg",
                               chrom="chr1",
                               start=100,
                               end=170,
                               strand="+",
                               sequence="A"*70,
                               snor_type="CD",
                               d_boxes="60..63",
                               c_boxes="50..56, 20..26, 30..36")

    def test_for_plexy_string(self):
        cdsnor1 = CD_snoRNA(snor_id="test_snoRNA",
                   organism="hg",
                   chrom="chr1",
                   start=100,
                   end=170,
                   strand="+",
                   sequence="A"*70,
                   snor_type="CD",
                   d_boxes="61..64",
                   c_boxes=None)
        plexy1 = ">test_snoRNA-_NNNNNNN_10_AAAA_61\n%s\n" % ("A"*70)
        self.assertEqual(cdsnor1.get_plexy_string(), plexy1)

        cdsnor2 = CD_snoRNA(snor_id="test_snoRNA",
                   organism="hg",
                   chrom="chr1",
                   start=100,
                   end=170,
                   strand="+",
                   sequence="A"*70,
                   snor_type="CD",
                   d_boxes="30..33,60..63",
                   c_boxes=None)
        plexy2 = ">test_snoRNA-_NNNNNNN_50_AAAA_30_NNNNNNN_10_AAAA_60\n%s\n" \
                 % ("A"*70)
        self.assertEqual(cdsnor2.get_plexy_string(), plexy2)

        cdsnor3 = CD_snoRNA(snor_id="test_snoRNA",
                   organism="hg",
                   chrom="chr1",
                   start=100,
                   end=170,
                   strand="+",
                   sequence="A"*70,
                   snor_type="CD",
                   d_boxes="60..63",
                   c_boxes="21..27")
        plexy3 = ">test_snoRNA-_AAAAAAA_21_AAAA_60\n%s\n" \
                 % ("A"*70)
        self.assertEqual(cdsnor3.get_plexy_string(), plexy3)

        cdsnor4 = CD_snoRNA(snor_id="test_snoRNA",
                   organism="hg",
                   chrom="chr1",
                   start=100,
                   end=170,
                   strand="+",
                   sequence="A"*70,
                   snor_type="CD",
                   d_boxes="30..33,60..63",
                   c_boxes="21..27")
        plexy4 = ">test_snoRNA-_NNNNNNN_50_AAAA_30_AAAAAAA_21_AAAA_60\n%s\n" \
                 % ("A"*70)
        self.assertEqual(cdsnor4.get_plexy_string(), plexy4)

        cdsnor5 = CD_snoRNA(snor_id="test_snoRNA",
                   organism="hg",
                   chrom="chr1",
                   start=100,
                   end=170,
                   strand="+",
                   sequence="A"*70,
                   snor_type="CD",
                   d_boxes="60..63",
                   c_boxes="49..55,21..27")
        plexy5 = ">test_snoRNA-_AAAAAAA_21_AAAA_60\n%s\n" \
                 % ("A"*70)
        self.assertEqual(cdsnor5.get_plexy_string(), plexy5)

        cdsnor6 = CD_snoRNA(snor_id="test_snoRNA",
                   organism="hg",
                   chrom="chr1",
                   start=100,
                   end=170,
                   strand="+",
                   sequence="A"*70,
                   snor_type="CD",
                   d_boxes="30..33,60..63",
                   c_boxes="49..55,21..27")
        plexy6 = ">test_snoRNA-_AAAAAAA_49_AAAA_30_AAAAAAA_21_AAAA_60\n%s\n" \
                 % ("A"*70)
        self.assertEqual(cdsnor6.get_plexy_string(), plexy6)

    def test_interaction_region_d_box(self):
        inreg = self.snor.get_d_interaction_region(length=20)
        self.assertEqual(type(inreg), type(Seq("AAAA")))
        self.assertEqual(str(inreg), "A"*20)

    def test_interaction_region_dprime_box(self):
        inreg = self.snor.get_dprime_interaction_region(length=20)
        self.assertEqual(type(inreg), type(Seq("AAAA")))
        self.assertEqual(str(inreg), "T"*20)



if __name__ == '__main__':
    unittest.main()
