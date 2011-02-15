import csb.test as test
import numpy as np

@test.unit
class TestMath(test.Case):

    def testDihedralAngles(self):
        from csb.math import dihedral_angle

        a = np.array([1, 0., 0.])
        b = np.array([0, 0., 0.])
        c = np.array([0, 1., 0.])
        d = np.array([0, 1., -1.])

        self.assertEqual(dihedral_angle(a, b, c, d), 90.0)
        self.assertEqual(dihedral_angle(a, b, c, d + a), 45.0)
        self.assertEqual(dihedral_angle(a, b, c, a), 0.0)
        self.assertEqual(dihedral_angle(a, b, c, -d), -90.0)

    def testLog(self):
        from csb.math import log, LOG_MAX, LOG_MIN
        from numpy import log as ref_log

        x = np.linspace(LOG_MIN, LOG_MAX, 1000000)
        
        self.assertTrue((log(x) == ref_log(x)).all() )
        self.assertEqual(log(10 * LOG_MAX), log(LOG_MAX))
        self.assertEqual(log(0.1 * LOG_MIN), log(LOG_MIN))


    def testExp(self):
        from csb.math import exp, EXP_MAX, EXP_MIN
        from numpy import exp as ref_exp
        
        x = np.linspace(EXP_MIN, EXP_MAX, 1000000)

        self.assertTrue((exp(x) == ref_exp(x)).all() )
        self.assertEqual(exp(EXP_MAX + 10.), exp(EXP_MAX))
        self.assertEqual(exp(10. * EXP_MIN), exp(EXP_MIN))

    def testPolar(self):
        from csb.math import polar, from_polar

        rand = np.random.random((10000,3))
        result = np.array(map(polar, map(from_polar, rand)))

        eps = 1e-8
        self.assertTrue((np.abs(rand - result) < eps).all())
        self.assertTrue((np.abs(polar(np.array([0.,1.0])) -\
                                np.array([1., np.pi * 0.5])) < eps).all())

    def testFromPolar(self):
        from csb.math import polar, from_polar

        rand = np.random.random((10000,3))
        result = np.array(map(from_polar, map(polar, rand)))

        eps = 1e-8
        self.assertTrue((np.abs(rand - result) < eps).all())
        self.assertTrue((np.abs(from_polar(np.array([1., np.pi * 0.5])) -\
                                np.array([0., 1.0])) < eps).all())

    def testRadian2Degree(self):
        from csb.math import degree2radian, radian2degree

        rand = np.random.random(10000) * 2 * np.pi
        eps = 1e-8
        
        self.assertTrue((np.abs(rand - 
                                degree2radian(radian2degree(rand)))\
                         < eps).all())


        
    def testDegree2Radian(self):
        from csb.math import degree2radian, radian2degree

        rand = np.random.random(10000) * 2 * np.pi
        eps = 1e-8
        
        self.assertTrue((np.abs(rand -
                                degree2radian(radian2degree(rand)))\
                         < eps).all())

        
        
        
if __name__ == '__main__':
    test.Console()
