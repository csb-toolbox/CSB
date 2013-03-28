import numpy
import csb.test as test

from csb.statistics.ars import ARS, Gauss


@test.functional
class TestARS(test.Case):

    def testNormal(self):
        mu = 5.
        sigma = 1.
        ars = ARS(Gauss(mu, sigma))
        ars.initialize([mu - 1., mu + 1.1], z0=-10., zmax=30)
        samples = numpy.array([ars.sample() for i in range(10000)])

        self.assertAlmostEqual(mu, numpy.mean(samples), delta=0.5)
        self.assertAlmostEqual(sigma, numpy.std(samples), delta=0.5)
    
           
if __name__ == '__main__':
    
    test.Console()
