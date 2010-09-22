import unittest
import numpy as np

class TestMath(unittest.TestCase):


    def test_dihedral_angles(self):
        from csb.math import dihedral_angle

        
        a = np.array([1,0.,0.])
        b = np.array([0,0.,0.])
        c = np.array([0,1.,0.])
        d = np.array([0,1.,-1.])

        self.assertEqual(dihedral_angle(a,b,c,d), 90.0)
        self.assertEqual(dihedral_angle(a,b,c,d + a), 45.0)
        self.assertEqual(dihedral_angle(a,b,c,a), 0.0)
        self.assertEqual(dihedral_angle(a,b,c,-d), -90.0)

    def test_log(self):
        from csb.math import log, LOG_MAX, LOG_MIN
        from numpy import log as ref_log

        x = np.linspace(LOG_MIN, LOG_MAX, 1000000)
        
        self.assertTrue((log(x) == ref_log(x)).all() )
        self.assertEqual(log(10 * LOG_MAX),log(LOG_MAX))
        self.assertEqual(log(0.1 * LOG_MIN),log(LOG_MIN))


    def test_exp(self):
        from csb.math import exp, EXP_MAX, EXP_MIN
        from numpy import exp as ref_exp
        
        x = np.linspace(EXP_MIN, EXP_MAX, 1000000)

        self.assertTrue((exp(x) == ref_exp(x)).all() )
        self.assertEqual(exp(EXP_MAX + 10.),exp(EXP_MAX))
        self.assertEqual(exp(10. * EXP_MIN),exp(EXP_MIN))

    def test_polar(self):
        from csb.math import polar, from_polar

        rand = np.random.random(10000)
        eps = 1e-8
        self.assertTrue((np.abs(rand - from_polar(polar(rand))) < eps).all())


    def test_from_polar(self):
        from csb.math import polar, from_polar

        rand = np.random.random(10000)
        eps = 1e-8
        self.assertTrue((np.abs(rand - from_polar(polar(rand))) < eps).all())

    def test_radian2degree(self):
        from csb.math import degree2radian, radian2degree

        rand = np.random.random(10000) * 2 * np.pi
        eps = 1e-8
        
        self.assertTrue((np.abs(rand - degree2radian(radian2degree(rand))) < eps).all())

        
    def test_degree2radian(self):
        from csb.math import degree2radian, radian2degree

        rand = np.random.random(10000) * 2 * np.pi
        eps = 1e-8
        
        self.assertTrue((np.abs(rand - degree2radian(radian2degree(rand))) < eps).all())

        
        
        
if __name__ == '__main__':
    unittest.main()
