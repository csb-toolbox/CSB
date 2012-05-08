import numpy
import csb.test as test
import csb.statistics.pdf

from csb.bio.io.wwpdb import LegacyStructureParser
from csb.bio.utils import probabilistic_fit
from csb.statistics.scalemixture import ScaleMixture, GammaPrior, InvGammaPrior


@test.functional
class TestScaleMixture(test.Case):

    def test_random(self):
               
        mixture = ScaleMixture(scales = [1.,1.,1.,1.],
                               prior = GammaPrior(), d = 3)

        s = mixture.random()
        samples = mixture.random(10000)

        mu = numpy.mean(samples)
        var = numpy.var(samples)
        
        self.assertWithinDelta(0.0, mu, delta=1e-1)
        self.assertWithinDelta(1., var, delta=1e-1)


    def test_logprob(self):

        x = numpy.linspace(-5,5,1000)

        normal = csb.statistics.pdf.Normal()

        mixture = ScaleMixture(scales = [1.,1.,1.,1.],
                               prior = None, d = 1)

        px = mixture(x)
        gx = normal(x)

        for i in range(len(px)):
            self.assertWithinDelta(px[i], 4 * gx[i], delta=1e-1)

        
    
    def test_gamma(self):
        """
        The posterior of a gaussian scale mixture with gamma prior
        is a Student's t distribution, with parameters alpha and beta.

        Give enough samples, we shoud be able to estimate these parameters
        """
        pdbfile = self.config.getTestFile('ake-xray-ensemble-ca.pdb')
        ensemble = LegacyStructureParser(pdbfile).parse_models()
        X = numpy.array(ensemble[0].list_coordinates(['CA'], True))
        Y = numpy.array(ensemble[13].list_coordinates(['CA'], True))

               
        mixture = ScaleMixture(scales = X.shape[0],
                               prior = GammaPrior(), d = 3)

        from csb.bio.utils import rmsd, fit

        R,t = fit(X,Y)
        
        # gibbs sampling cycle
        for i in range(200):
            # apply rotation
            data = numpy.sum((X - numpy.dot(Y, numpy.transpose(R)) - t)**2, axis=-1)**(1./2)
            # sample scales
            mixture.estimate(data)
            # sample rotations
            R,t = probabilistic_fit(X,Y, mixture.scales)

        self.assertEqual(mixture.scales.shape, (211,))
        
        R_opt = numpy.eye(3)
        t_opt = numpy.zeros((3,))
        
        for i in range(3):
            self.assertWithinDelta(t[i], t_opt[i], delta = 2.)
            for j in range(3):
                self.assertWithinDelta(R_opt[i,j], R[i,j], delta = 1e-1)


    def test_invgamma(self):
        """
        The posterior of a gaussian scale mixture with gamma prior
        is a Student's t distribution, with parameters alpha and beta.

        Give enough samples, we shoud be able to estimate these parameters
        """
        pdbfile = self.config.getTestFile('ake-xray-ensemble-ca.pdb')
        ensemble = LegacyStructureParser(pdbfile).parse_models()
        X = numpy.array(ensemble[0].list_coordinates(['CA'], True))
        Y = numpy.array(ensemble[13].list_coordinates(['CA'], True))

               
        mixture = ScaleMixture(scales = X.shape[0],
                               prior = InvGammaPrior(), d = 3)

        from csb.bio.utils import rmsd, fit

        R,t = fit(X,Y)
        # gibbs sampling cycle
        for i in range(200):
            # apply rotation
            data = numpy.sum((X - numpy.dot(Y, numpy.transpose(R)) - t)**2, axis=-1)**(1./2)
            # sample scales
            mixture.estimate(data)
            # sample rotations
            R,t = probabilistic_fit(X,Y, mixture.scales)
        
        self.assertEqual(mixture.scales.shape, (211,))
        
        R_opt = numpy.eye(3)
        t_opt = numpy.zeros((3,))
        
        for i in range(3):
            self.assertWithinDelta(t[i], t_opt[i], delta = 2.)
            for j in range(3):
                self.assertWithinDelta(R_opt[i,j], R[i,j], delta = 1e-1)


            
if __name__ == '__main__':
    
    test.Console()
        
