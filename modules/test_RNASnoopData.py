import unittest
from RNASnoopData import RNASnoopStructure


class TestRNASnoopStructure(unittest.TestCase):

    """Test the interation structure of snoRNAs"""

    def setUp(self):
        # (t_i_gap, U_gap, i_b_gap, i_t_gap, gap_right, stem_length,
        # stem_asymmetry);
        self.structures = (('<<<|.<<<<<<<&..................(((((>>>>>>>((((((((((..........))).)))))........)).>>>...)))))', [1, 2, 3, 0, 0, 39, -9]),
                           ('<<.|.<<<&>>>.(((..((((.(((((((((((....)))..))))..)))).))))...)))>>.........', [0, 3, -1, 1, 0, 51, -5]),
                           ('<.<.|.<<<<&(>>>>..(((....(((((((.(((((.......)))))))))))))))..>.>.)......', [2, 3, 1, 2, 0, 42, 5]),
                           ('<<<|.<<<<<<<&..................(((((>>>>>>>((((((((((..........))).)))))........)).>>>...)))))', [1, 2, 3, 0, 0, 39, -9]),
                           ('<<.|.<<<&>>>.(((..((((.(((((((((((....)))..))))..)))).))))...)))>>.........', [0, 3, -1, 1, 0, 51, -5]),
                           ('<|.<<<<<<<<<<&>>>>>>>>>>.((.(..(((..(((((((((((........))))))))))).))).)))>...........', [0, 2, -1, 1, 0, 49, 3]),
                           ('<<<.<<<<.<.|.<<&.>>.(((..((((.(((((((((((....)))..))))..)))).))))...)))>.>>>>.>>>.', [0, 3, -1, 1, 0, 51, -5]),
                           ('<<<.<<<|.<.<<<&...(..>>>.>.(((..(((..(((((((((((........))))))))))).)))..)))>>>.>>>..).', [0, 2, 2, 1, 1, 49, 1]),
                           ('<<<|.<<<<&((..>>>>..((..((((((.(((.(((......))).))).))))))..))>>>.......)).', [0, 2, 7, 2, 0, 42, 0]),
                           ('<<.<<|.<<<<&>>>>..(((.(((...........((((.(((((.((....)).)))))...))))))).)))..>>.>>......', [2, 2, -1, 2, 0, 57, 9]),
                           ('<<.<|.<<<<<<&((((((...>>>>>>(((((.(((.(((......))).))).))))).>.>>)))))).......', [1, 2, 0, 0, 0, 32, 0]),
                           ('<<<<<<|.<&(.........>..((((....((.((((.(((((.((....)).)))))...)))).)).)).))>>>>>>)....', [0, 2, 0, 2, 0, 52, 0]),
                           ('<<.<<.|.<<<<<&>>>>>.(((((.((((.(((..............))).)))).))).)).>>.>>............', [1, 3, -1, 1, 0, 43, -1]),
                           ('<<<<|.<<&(.((.....>>.((.(((.....(((((..(((((.......)))))..)))))..))).))>>>>)).)....', [0, 2, 0, 1, 0, 50, 3]),
                           ('<<<<.<.|.<.<<<<<&..((.............>>>>>.>.(.(((((((..(((.....)))))))))).)>.>>>>..)).', [0, 3, 2, 1, 1, 31, 2]),
                           ('<<.<<<<|.<<<<&(((.....>>>>..(((......(((((..(((((.......)))))..)))))......)))>>>>.>>.)))', [0, 2, 1, 2, 0, 49, 0]),
                           ('<<<<<<<<|.<<<<.<&...(((>.>>>>..(((((.((.........))...)))))>>>>>>>>)))', [0, 2, 0, 2, 1, 27, -2]),
                           ('<<<.<|.<<<&.....((..........>>>.(((((((......(((((((((((((.((....)).))))))))))))))))))))>.>>>.))...', [0, 2, 1, 1, 0, 56, 6]),
                           ('<<<<<|.<.<<<<<<<<<<<.<<&(.>>.>>>>>>>>>>>.>..(((..((((.((....)).)))).)))>>>>>......)....', [0, 2, 6, 2, 2, 27, 1]),
                           ('<<<<<<.<|.<<<<&.((>>>>.(((...........((((((((((.........))))))..))))))).>.>>>>>>.))......', [1, 2, 1, 1, 0, 48, 9]),
                           ('<<<.|.<&(.>..((..(((...((((((((....))))))))...)))....))>>>........)....', [0, 3, 8, 2, 0, 42, -2]),
                           ('<<<|.<<<<&.((.....>>>>.((.((.((.((((((((((.........))))))..)))).))...))))>>>))......', [0, 2, 0, 1, 0, 50, -3]))

    def test_gap_right(self):
        for structure, output in self.structures:
            model = RNASnoopStructure(structure)
            self.assertEqual(model.gap_right, output[4])

    def test_U_gap(self):
        for structure, output in self.structures:
            model = RNASnoopStructure(structure)
            self.assertEqual(model.U_gap, output[1])

    def test_t_i_gap(self):
        models = [RNASnoopStructure(s[0]) for s in self.structures]
        for model, (structure, output) in zip(models, self.structures):
            self.assertEqual(model.t_i_gap, output[0])

    def test_i_b_gap(self):
        models = [RNASnoopStructure(s[0]) for s in self.structures]
        for model, (structure, output) in zip(models, self.structures):
            self.assertEqual(model.i_b_gap, output[2])

    def test_i_t_gap(self):
        for structure, output in self.structures:
            model = RNASnoopStructure(structure)
            self.assertEqual(model.i_t_gap, output[3])

    def test_stem_length(self):
        for structure, output in self.structures:
            model = RNASnoopStructure(structure)
            self.assertEqual(model.stem_length, output[5])

    # def test_stem_asymmetry(self):
    #     for structure, output in self.structures:
    #         model = RNASnoopStructure(structure)
    #         self.assertEqual(model.stem_asymmetry, output[6])

if __name__ == '__main__':
    unittest.main()

