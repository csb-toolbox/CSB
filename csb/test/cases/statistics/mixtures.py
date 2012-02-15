
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

            # positions of maxima in membership matrix
            w = argmax(m.Z, 1)

            # position numbers might be permutated, so count equal pairs
            ww = zip(w, self.w_ref_segments)
            same = sum(sorted(ww.count(i) for i in set(ww))[-m.K:])

            # we want at least 90% overlap
            if float(same) / X.shape[1] >= 0.9:
                break
        else:
            self.assertTrue(False, 'segmentation not reproduced within 5 iterations')

    def testConformerMixture(self):

        X = self._ake_ensemble_coords()

        m = mixtures.ConformerMixture(X.shape[1], X.shape[0], 3)

        # algorithms with randomized initialization, try multiple times
        for _ in xrange(10):

            # fit parameters
            m.em(X)

            # positions of maxima in membership matrix
            w = argmax(m.Z, 1)

            # position numbers might be permutated, so count equal pairs
            ww = zip(w, self.w_ref_conformers)
            same = sum(sorted(ww.count(i) for i in set(ww))[-m.K:])

            # we want at least 14 out of 16
            if X.shape[0] - same < 3:
                break
        else:
            self.assertTrue(False, 'conformers not reproduced within 10 iterations')

if __name__ == '__main__':

    test.Console()

# vi:expandtab:smarttab
