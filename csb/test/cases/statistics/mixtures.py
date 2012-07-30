from numpy import array, linspace

from csb import test
from csb.bio.io.wwpdb import LegacyStructureParser
from csb.statistics import mixtures


@test.functional
class TestMixtures(test.Case):

    w_ref_segments = array([
        2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
        2, 2, 2, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2,
        2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
        2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 1, 1,
        1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
        1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
        1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
        2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
        2, 2, 2, 2])

    w_ref_conformers = array([2, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 1, 1])

    def _ake_ensemble_coords(self):

        pdbfile = self.config.getTestFile('ake-xray-ensemble-ca.pdb')
        ensemble = LegacyStructureParser(pdbfile).parse_models()
        X = array([model.get_coordinates(['CA'], True) for model in ensemble])

        self.assertEqual(X.shape, (16, 211, 3))

        self._ake_ensemble_coords = lambda: X

        return X
    
    def testSegmentMixture(self):

        self._testMixture(mixtures.SegmentMixture, self.w_ref_segments)

    def testConformerMixture(self):

        self._testMixture(mixtures.ConformerMixture, self.w_ref_conformers, 14./16.)

    def _testMixture(self, cls, w_ref, min_overlap=0.9, repeats=5):

        X = self._ake_ensemble_coords()
        K = len(set(w_ref))

        # non-randomized heuristic with BIC
        m = cls.new(X)
        overlap = m.overlap(w_ref)

        self.assertTrue(overlap >= min_overlap, 'mixture not reproduced with heuristic')

        # annealing (randomized initialization)
        m = cls(X, K, False)
        for _ in range(repeats):
            m.randomize_scales()
            m.anneal(linspace(2.0, 0.1, 10))

            overlap = m.overlap(w_ref)
            if overlap >= min_overlap:
                break
        else:
            self.assertTrue(False, 'mixture not reproduced with annealing')


if __name__ == '__main__':

    test.Console()

# vi:expandtab:smarttab
