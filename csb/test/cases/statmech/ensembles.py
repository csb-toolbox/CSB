import numpy

import csb.numeric
import csb.test as test

from csb.statmech.ensembles import BoltzmannEnsemble, TsallisEnsemble, CompositeEnsemble


@test.functional
class TestEnergy(test.Case):

    def testBoltzmann(self):
        e = numpy.linspace(-50, 1000, 1000)

        be = BoltzmannEnsemble(beta=1,)
        te = be.energy(e)
        
        for i in range(len(e)):
            self.assertEqual(e[i], te[i])

        be = BoltzmannEnsemble(beta=0.001,)
        te = be.energy(e)
        
        for i in range(len(e)):
            self.assertEqual(e[i] * 0.001, te[i])

    def testTsallis(self):
        e = numpy.linspace(-50, 1000, 1000)

        tsallis = TsallisEnsemble(q=1.,)
        te = tsallis.energy(e)
        
        for i in range(len(e)):
            self.assertEqual(e[i], te[i])

        tsallis = TsallisEnsemble(q=1.1, e_min= -50.)
        te = tsallis.energy(e)
        q = 1.1
        ee = q / (q - 1.) * csb.numeric.log(1 + (q - 1) * (e + 50.)) - 50
        
        for i in range(len(e)):
            self.assertAlmostEqual(ee[i], te[i], delta=1e-5)


    def testComposite(self):
        e1 = numpy.linspace(-50, 1000, 1000)
        e2 = numpy.linspace(-30, 3000, 1000)

        q = 1.1
        beta = 0.1
        ee = q / (q - 1.) * csb.numeric.log(1 + (q - 1) * (e1 + 50.)) - 50
        ee += e2 * beta

        ce = CompositeEnsemble([TsallisEnsemble(q=q, e_min= -50.),
                                BoltzmannEnsemble(beta=beta,)])

        cee = ce.energy([e1, e2])
        
        for i in range(len(e1)):
            self.assertAlmostEqual(ee[i], cee[i], delta=1e-5)

                
if __name__ == '__main__':
    
    test.Console()

