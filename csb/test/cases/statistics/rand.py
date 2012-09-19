import numpy
import warnings

import csb.test as test

from csb.numeric import exp, log_sum_exp, log
from csb.statistics.rand import truncated_gamma, truncated_normal, sample_from_histogram
from csb.statistics.pdf import Normal
from csb.statistics import density


@test.functional
class TestRand(test.Case):

    def testTruncatedGamma(self):
        alpha = 2.
        beta = 1.
        x_min = 0.1
        x_max = 5.

        x = truncated_gamma(10000, alpha, beta, x_min, x_max)

        self.assertTrue((x <= x_max).all())
        self.assertTrue((x >= x_min).all())

        hy, hx = density(x, 100)
        hx = 0.5 * (hx[1:] + hx[:-1])
        hy = hy.astype('d')

        with warnings.catch_warnings(record=True) as warning:
            warnings.simplefilter("always")            
            
            hy /= (hx[1] - hx[0]) * hy.sum()
            
            self.assertLessEqual(len(warning), 1)
            
            if len(warning) == 1:
                warning = warning[0]
                self.assertEqual(warning.category, RuntimeWarning)
                self.assertTrue(str(warning.message).startswith('divide by zero encountered'))            

        x = numpy.linspace(x_min, x_max, 1000)
        p = (alpha - 1) * log(x) - beta * x
        p -= log_sum_exp(p)
        p = exp(p) / (x[1] - x[0])

    def testTruncatedNormal(self):

        mu = 2.
        sigma = 1.
        x_min = -1.
        x_max = 5.

        x = truncated_normal(10000, mu, sigma, x_min, x_max)

        self.assertAlmostEqual(numpy.mean(x), mu, delta=1e-1)
        self.assertAlmostEqual(numpy.var(x), sigma, delta=1e-1)

        self.assertTrue((x <= x_max).all())
        self.assertTrue((x >= x_min).all())

        hy, hx = density(x, 100)
        hx = 0.5 * (hx[1:] + hx[:-1])
        hy = hy.astype('d')
        
        with warnings.catch_warnings(record=True) as warning:
            warnings.simplefilter("always")        
            
            hy /= (hx[1] - hx[0]) * hy.sum()
            
            self.assertLessEqual(len(warning), 1)
            
            if len(warning) == 1:
                warning = warning[0]
                self.assertEqual(warning.category, RuntimeWarning)
                self.assertTrue(str(warning.message).startswith('divide by zero encountered'))
            
        x = numpy.linspace(mu - 5 * sigma, mu + 5 * sigma, 1000)

        p = -0.5 * (x - mu) ** 2 / sigma ** 2
        p -= log_sum_exp(p)
        p = exp(p) / (x[1] - x[0])
    
    

    def testSampleFromHistogram(self):
        mu = 5.
        sigma = 1.

        normal = Normal(mu, sigma)

        x = normal.random(10000)
        hx, p = density(x, 100)

        samples = hx[sample_from_histogram(p, n_samples=10000)]

        self.assertAlmostEqual(mu, numpy.mean(samples), delta=0.5)
        self.assertAlmostEqual(sigma, numpy.std(samples), delta=0.5)


if __name__ == '__main__':

    test.Console()
