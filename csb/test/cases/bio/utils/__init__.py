
import sys
import numpy
import multiprocessing

import csb.test as test
import csb.bio.utils as cbu
import csb.io


X1 = numpy.array([
    [ 0.,  0.,  0.],
    [ 1.,  0.,  0.],
    [ 0.,  1.,  0.]])

X2 = numpy.array([
    [ 0.,  0.,  0.],
    [ 1.,  2.,  0.],
    [-2., -1.,  0.]])

X3 = numpy.array([
    [ 0.,  0.,  0.],
    [ 2., -1.,  0.],
    [-1.,  2.,  0.]])

RZ = numpy.array([
    [ 0.,  1.,  0.],
    [-1.,  0.,  0.],
    [ 0.,  0.,  1.]])

X4 = numpy.array([
    [ 0.,  0.,  0.],
    [ 1.,  0.,  0.],
    [ 0.,  1.,  0.],
    [ 1.,  1.,  0.]])

X5 = numpy.array([
    [   0.,    0.,    0.],
    [ 100.,    0.,    0.],
    [   0.,  100.,    0.],
    [  50.,   50.,    0.]])

X6 = numpy.array([
    [   0.,    0.,    0.],
    [ 100.,    0.,    0.],
    [   0.,  100.,    0.],
    [  60.,   60.,    0.]])

            
@test.regression
class Regressions(test.Case):

    def _timeoutTest(self):
        cbu.tm_superimpose([[1, 1, 1]], [[1, 1, 1]])
        
    def _multiprocessingTest(self):        
        return True
    
    def _runProcess(self, target, timeout=1.0):
                
        p = multiprocessing.Process(target=target)
        p.start()
        p.join(timeout=timeout)
        
        return p
    
    @test.skip("n/a on this platform", sys.platform.startswith('win'))        
    def testTMSuperimpose(self):
        """
        @see: [CSB 0000058]
        """
        try:
            self._runProcess(target=self._multiprocessingTest)
        except:
            self.skipTest("may produce a false positive")
                
        p = self._runProcess(target=self._timeoutTest, timeout=5.0)
        
        if p.is_alive():
            p.terminate()
            self.fail('timeout expired')
            

@test.functional
class TestUtils(test.Case):

    def assertArrayEqual(self, first, second, eps=1e-7):
        diff = numpy.asarray(first) - numpy.asarray(second)
        self.assertTrue((abs(diff) < eps).all())

    def testFit(self):
        R, t = cbu.fit(X1, X2)
        Y = numpy.dot(X2, R.T) + t

        self.assertArrayEqual(R, RZ)
        self.assertArrayEqual(t, [0., 0., 0.])
        self.assertArrayEqual(Y, X3)

    def testWFit(self):
        w = numpy.array([1., 1., 0.])
        R, t = cbu.wfit(X1, X2, w)                              #@UnusedVariable

        d = 5.0**0.5
        self.assertArrayEqual(t, [-d / 2.0 + 0.5, 0., 0.])

    def testFitWellordered(self):
        R, t = cbu.fit_wellordered(X5, X6, 10, 1.0)             #@UnusedVariable

        self.assertArrayEqual(t, [0., 0., 0.])

    def testRmsd(self):
        rmsd = cbu.rmsd(X1, X2)

        self.assertAlmostEqual(rmsd, (4./3.)**0.5)

    def testWrmsd(self):
        w = numpy.array([1., 1., 0.])
        rmsd = cbu.wrmsd(X1, X2, w)

        d = 5.0**0.5
        self.assertAlmostEqual(rmsd, d / 2.0 - 0.5)

    def testTorsionRmsd(self):
        rmsd = cbu.torsion_rmsd(X1[:,:2], X1[:,:2])

        self.assertAlmostEqual(rmsd, 0.0)

    def testTmScore(self):
        score = cbu.tm_score(X1, X3)

        self.assertAlmostEqual(score, 0.4074, 4)

    def testTmSuperimpose(self):
        R, t, score = cbu.tm_superimpose(X1, X2)            #@UnusedVariable

        self.assertAlmostEqual(score, 0.4074, 4)

    def testCenterOfMass(self):
        com = cbu.center_of_mass(X4)

        self.assertArrayEqual(com, [0.5, 0.5, 0.0])

    def testRadiusOfGyration(self):
        gyradius = cbu.radius_of_gyration(X4)

        s2 = 2.0**0.5
        self.assertArrayEqual(gyradius, s2 / 2.0)

    def testSecondMoments(self):
        sm = cbu.second_moments(X1)

        # TODO: correct?
        sm_test = numpy.array([
            [ 2./3., -1./3., 0.    ],
            [-1./3.,  2./3., 0.    ],
            [ 0.,     0.,    0.    ]])
        self.assertArrayEqual(sm, sm_test)

    def testInertiaTensor(self):
        it = cbu.inertia_tensor(X1)

        # TODO: correct?
        it_test = numpy.array([
            [ 2./3.,  1./3., 0.    ],
            [ 1./3.,  2./3., 0.    ],
            [ 0.,     0.,    4./3. ]])
        self.assertArrayEqual(it, it_test)

    def testDistanceMatrix(self):
        d = cbu.distance_matrix(X1)

        s2 = 2.0**0.5
        d_test = [
            [ 0., 1., 1. ],
            [ 1., 0., s2 ],
            [ 1., s2, 0. ]]
        self.assertArrayEqual(d, d_test)

    def testDistance(self):
        d = cbu.distance(X1, X2)

        self.assertEqual(d.shape, (len(X1),))
        self.assertArrayEqual(d[:2], [0., 2.])

    def testRmsdCur(self):
        rmsd = cbu.rmsd_cur(X1, X2)

        self.assertAlmostEqual(rmsd, 2.0)

if __name__ == '__main__':

    test.Console()

# vi:expandtab:smarttab

