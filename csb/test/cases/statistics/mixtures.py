
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

        self._testMixture(mixtures.SegmentMixture, self.w_ref_segments)

    def testSegmentMixture2(self):
        
        self._testMixture(mixtures.SegmentMixture2, self.w_ref_segments)

    def testConformerMixture(self):

        self._testMixture(mixtures.ConformerMixture, self.w_ref_conformers, 14./16.)

    def _testMixture(self, cls, w_ref, min_overlap=0.9, repeats=5):

        X = self._ake_ensemble_coords()
        K = len(set(w_ref))

        # algorithms with randomized initialization, try multiple times
        for _ in xrange(repeats):
            m = cls.from_coords(X, K, randomize=1)

            overlap = m.overlap(w_ref)
            if overlap >= min_overlap:
                break
        else:
            self.assertTrue(False, 'mixture not reproduced within %d iterations' % (repeats))

        # non-randomized heuristic with BIC
        m = cls.from_coords(X)
        overlap = m.overlap(w_ref)

        self.assertTrue(overlap >= min_overlap, 'mixture not reproduced with heuristic')

if __name__ == '__main__':

    test.Console()

# vi:expandtab:smarttab
