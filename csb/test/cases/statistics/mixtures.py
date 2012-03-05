
from numpy import array, argmax

from csb import test

from csb.bio.io.wwpdb import LegacyStructureParser
from csb.statistics import mixtures

def get_coords(chain, what=['CA']):
    coords = []
    for residue in chain.residues:
        if not residue.has_structure:
            continue
        try:
            for atom_kind in what:
                coords.append(residue.atoms[atom_kind].vector)
        except csb.pyutils.ItemNotFoundError:
            continue
    return coords

def get_ensemble_coords(ensemble, what=['CA']):
    X = [get_coords(model.first_chain) for model in ensemble.models]
    return array(X)


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
        X = get_ensemble_coords(ensemble)

        self.assertEqual(X.shape, (16, 211, 3))

        self._ake_ensemble_coords = lambda: X

        return X
    
    def testSegmentMixture(self):

        self._testSegmentMixture(mixtures.SegmentMixture)

    def testSegmentMixture2(self):
        
        self._testSegmentMixture(mixtures.SegmentMixture2)

    def _testSegmentMixture(self, cls):

        X = self._ake_ensemble_coords()

        m = cls(X.shape[1], X.shape[0], 3)

        # algorithms with randomized initialization, try multiple times
        for _ in xrange(5):

            # fit parameters
            m.em(X)

            # we want at least 90% overlap
            same = m.overlap(self.w_ref_segments)
            if same >= 0.9:
                break
        else:
            self.assertTrue(False, 'segmentation not reproduced within 5 iterations')

        self._testHeuristic(cls, self.w_ref_segments)

    def _testHeuristic(self, cls, w_ref):

        X = self._ake_ensemble_coords()

        m = cls.from_coords(X)
        same = m.overlap(w_ref)

        self.assertTrue(same > 0.8, 'mixture not reproduced')

    def testConformerMixture(self):

        X = self._ake_ensemble_coords()

        m = mixtures.ConformerMixture(X.shape[1], X.shape[0], 3)

        # algorithms with randomized initialization, try multiple times
        for _ in xrange(10):

            # fit parameters
            m.em(X)

            # we want at least 14 out of 16
            same = m.overlap(self.w_ref_conformers)
            if same >= 14./16.:
                break
        else:
            self.assertTrue(False, 'conformers not reproduced within 10 iterations')

        self._testHeuristic(mixtures.ConformerMixture, self.w_ref_conformers)

if __name__ == '__main__':

    test.Console()

# vi:expandtab:smarttab
