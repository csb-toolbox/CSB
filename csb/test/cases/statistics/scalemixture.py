import numpy
import csb.test as test
import csb.statistics.pdf

from csb.bio.io.wwpdb import LegacyStructureParser
from csb.bio.utils import probabilistic_fit, average_structure
from csb.statistics.scalemixture import ScaleMixture, GammaPrior, InvGammaPrior
from csb.statistics.scalemixture import GammaPosteriorMAP, InvGammaPosteriorMAP


@test.functional
class TestScaleMixture(test.Case):

    def testRandom(self):
               
        mixture = ScaleMixture(scales=[1., 1., 1., 1.],
                               prior=GammaPrior(), d=3)

        mixture.random()
        samples = mixture.random(10000)

        mu = numpy.mean(samples)
        var = numpy.var(samples)
        
        self.assertWithinDelta(0.0, mu, delta=1e-1)
        self.assertWithinDelta(1., var, delta=1e-1)

    def testLogProb(self):

        x = numpy.linspace(-5, 5, 1000)

        normal = csb.statistics.pdf.Normal()

        mixture = ScaleMixture(scales=[1., 1., 1., 1.],
                               prior=None, d=1)

        px = mixture(x)
        gx = normal(x)

        for i in range(len(px)):
            self.assertWithinDelta(px[i], 4 * gx[i], delta=1e-1)

    def testGamma(self):
        """
        The posterior of a gaussian scale mixture with gamma prior
        is a Student's t distribution, with parameters alpha and beta.

        Give enough samples, we shoud be able to estimate these parameters
        """
        pdbfile = self.config.getTestFile('ake-xray-ensemble-ca.pdb')
        ensemble = LegacyStructureParser(pdbfile).parse_models()
        X = numpy.array(ensemble[0].list_coordinates(['CA'], True))
        Y = numpy.array(ensemble[13].list_coordinates(['CA'], True))

               
        mixture = ScaleMixture(scales=X.shape[0],
                               prior=GammaPrior(), d=3)

        from csb.bio.utils import fit

        R, t = fit(X, Y)
        #numpy.random.seed(100)
        # gibbs sampling cycle
        for i in range(200):
            # apply rotation
            data = numpy.sum((X - numpy.dot(Y, numpy.transpose(R)) - t) ** 2, axis= -1) ** (1. / 2)
            # sample scales
            mixture.estimate(data)
            # sample rotations
            R, t = probabilistic_fit(X, Y, mixture.scales)

        self.assertEqual(mixture.scales.shape, (211,))
        
        R_opt = numpy.eye(3)
        t_opt = numpy.zeros((3,))
        
        for i in range(3):
            self.assertWithinDelta(t[i], t_opt[i], delta=2.)
            for j in range(3):
                self.assertWithinDelta(R_opt[i, j], R[i, j], delta=1e-1)


    def testGammaMAP(self):
        """
        The posterior of a gaussian scale mixture with gamma prior
        is a Student's t distribution, with parameters alpha and beta.

        Give enough samples, we shoud be able to estimate these parameters
        """
        pdbfile = self.config.getTestFile('ake-xray-ensemble-ca.pdb')
        ensemble = LegacyStructureParser(pdbfile).parse_models()
        X = numpy.array(ensemble[0].list_coordinates(['CA'], True))
        Y = numpy.array(ensemble[13].list_coordinates(['CA'], True))

        prior = GammaPrior()
        prior.estimator = GammaPosteriorMAP()
        mixture = ScaleMixture(scales=X.shape[0],
                               prior=prior, d=3)

        from csb.bio.utils import fit, wfit
        R, t = fit(X, Y)
        #numpy.random.seed(100)
        # gibbs sampling cycle
        for i in range(200):
            # apply rotation
            data = numpy.sum((X - numpy.dot(Y, numpy.transpose(R)) - t) ** 2, axis= -1) ** (1. / 2)
            # sample scales
            mixture.estimate(data)
            # sample rotations
            R, t = wfit(X, Y, mixture.scales)

        self.assertEqual(mixture.scales.shape, (211,))
        
        R_opt = numpy.eye(3)
        t_opt = numpy.zeros((3,))
        
        for i in range(3):
            self.assertWithinDelta(t[i], t_opt[i], delta=2.)
            for j in range(3):
                self.assertWithinDelta(R_opt[i, j], R[i, j], delta=1e-1)


    def testInvGamma(self):
        """
        The posterior of a gaussian scale mixture with gamma prior
        is a Student's t distribution, with parameters alpha and beta.

        Give enough samples, we shoud be able to estimate these parameters
        """
        pdbfile = self.config.getTestFile('ake-xray-ensemble-ca.pdb')
        ensemble = LegacyStructureParser(pdbfile).parse_models()
        X = numpy.array(ensemble[0].list_coordinates(['CA'], True))
        Y = numpy.array(ensemble[13].list_coordinates(['CA'], True))

               
        mixture = ScaleMixture(scales=X.shape[0],
                               prior=InvGammaPrior(), d=3)

        from csb.bio.utils import fit

        R, t = fit(X, Y)
        #numpy.random.seed(100)
        # gibbs sampling cycle
        for i in range(200):
            # apply rotation
            data = numpy.sum((X - numpy.dot(Y, numpy.transpose(R)) - t) ** 2, axis= -1) ** (1. / 2)
            # sample scales
            mixture.estimate(data)
            # sample rotations
            R, t = probabilistic_fit(X, Y, mixture.scales)
        
        self.assertEqual(mixture.scales.shape, (211,))
        
        R_opt = numpy.eye(3)
        t_opt = numpy.zeros((3,))
        
        for i in range(3):
            self.assertWithinDelta(t[i], t_opt[i], delta=2.)
            for j in range(3):
                self.assertWithinDelta(R_opt[i, j], R[i, j], delta=1e-1)


    def testInvGammaMAP(self):
        """
        The posterior of a gaussian scale mixture with gamma prior
        is a Student's t distribution, with parameters alpha and beta.

        Give enough samples, we shoud be able to estimate these parameters
        """
        pdbfile = self.config.getTestFile('ake-xray-ensemble-ca.pdb')
        ensemble = LegacyStructureParser(pdbfile).parse_models()
        X = numpy.array(ensemble[0].list_coordinates(['CA'], True))
        Y = numpy.array(ensemble[13].list_coordinates(['CA'], True))


               
        prior = InvGammaPrior()
        prior.estimator = InvGammaPosteriorMAP()
        mixture = ScaleMixture(scales=X.shape[0],
                               prior=prior, d=3)

        from csb.bio.utils import fit, wfit

        R, t = fit(X, Y)
        #numpy.random.seed(100)
        # gibbs sampling cycle
        for i in range(200):
            # apply rotation
            data = numpy.sum((X - numpy.dot(Y, numpy.transpose(R)) - t) ** 2, axis= -1) ** (1. / 2)
            # sample scales
            mixture.estimate(data)
            # sample rotations
            R, t = wfit(X, Y, mixture.scales)
        
        self.assertEqual(mixture.scales.shape, (211,))
        
        R_opt = numpy.eye(3)
        t_opt = numpy.zeros((3,))
        
        for i in range(3):
            self.assertWithinDelta(t[i], t_opt[i], delta=2.)
            for j in range(3):
                self.assertWithinDelta(R_opt[i, j], R[i, j], delta=1e-1)


    def testEnsemble(self):
        """
        The posterior of a gaussian scale mixture with gamma prior
        is a Student's t distribution, with parameters alpha and beta.

        Give enough samples, we shoud be able to estimate these parameters
        """
        pdbfile = self.config.getTestFile('ake-xray-ensemble-ca.pdb')
        ensemble = LegacyStructureParser(pdbfile).parse_models()

        X = numpy.array([model.list_coordinates(['CA'], True) for model in ensemble])
        x_mu = average_structure(X)
        n =X.shape[1]
        m =X.shape[0]
        R = numpy.zeros((m,3,3))
        t = numpy.ones((m,3))

          
        prior = GammaPrior()

        mixture = ScaleMixture(scales=n, prior = prior, d=3)
                               

        from csb.bio.utils import fit, wfit

        for i in range(m):
            R[i,:,:], t[i,:] = fit(x_mu, X[i])
        
        # gibbs sampling cycle
        for j in range(200):
            # apply rotation
            data = numpy.array([numpy.sum((x_mu - numpy.dot(X[i], numpy.transpose(R[i])) - t[i]) **2, -1)**0.5
                                for i in range(m)]).T
            # sample scales
            mixture.estimate(data)
            # sample rotations
            for i in range(m):
                R[i,:,:], t[i,:] = wfit(x_mu, X[i], mixture.scales)


        self.assertEqual(mixture.scales.shape, (211,))
        
        R_opt = numpy.eye(3)
        t_opt = numpy.zeros((3,))
        for k in range(m):
            for i in range(3):
                self.assertWithinDelta(t[k,i], t_opt[i], delta=2.)
                for j in range(3):
                    self.assertWithinDelta(abs(R[k,i, j]), R_opt[i, j], delta=0.15)



            
if __name__ == '__main__':
    
    test.Console()
        
