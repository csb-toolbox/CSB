import csb.test as test
import numpy as np


@test.unit
class TestMath(test.Case):

    def testDihedralAngles(self):
        from csb.numeric import dihedral_angle

        a = np.array([1, 0., 0.])
        b = np.array([0, 0., 0.])
        c = np.array([0, 1., 0.])
        d = np.array([0, 1., -1.])

        self.assertEqual(dihedral_angle(a, b, c, d), 90.0)
        self.assertEqual(dihedral_angle(a, b, c, d + a), 45.0)
        self.assertEqual(dihedral_angle(a, b, c, a), 0.0)
        self.assertEqual(dihedral_angle(a, b, c, -d), -90.0)

    def testLog(self):
        from csb.numeric import log, LOG_MAX, LOG_MIN
        from numpy import log as ref_log

        x = np.linspace(LOG_MIN, LOG_MAX, 1000000)
        
        self.assertTrue((log(x) == ref_log(x)).all())
        self.assertEqual(log(10 * LOG_MAX), log(LOG_MAX))
        self.assertEqual(log(0.1 * LOG_MIN), log(LOG_MIN))


    def testExp(self):
        from csb.numeric import exp, EXP_MAX, EXP_MIN
        from numpy import exp as ref_exp
        
        x = np.linspace(EXP_MIN,
                        EXP_MAX, 100000)
        
        self.assertTrue((exp(x) == ref_exp(x)).all())
        self.assertEqual(exp(EXP_MAX + 10.), exp(EXP_MAX))
        self.assertEqual(exp(10. * EXP_MIN), exp(EXP_MIN))

    def testPolar(self):
        from csb.numeric import polar, from_polar

        rand = np.random.random((10000, 3))
        result = np.array(list(map(polar, map(from_polar, rand))))

        eps = 1e-8
        self.assertTrue((np.abs(rand - result) < eps).all())
        self.assertTrue((np.abs(polar(np.array([0., 1.0])) - \
                                np.array([1., np.pi * 0.5])) < eps).all())

    def testFromPolar(self):
        from csb.numeric import polar, from_polar

        rand = np.random.random((10000, 3))
        result = np.array(list(map(from_polar, map(polar, rand))))

        eps = 1e-8
        self.assertTrue((np.abs(rand - result) < eps).all())
        self.assertTrue((np.abs(from_polar(np.array([1., np.pi * 0.5])) - \
                                np.array([0., 1.0])) < eps).all())

    def testRadian2Degree(self):
        from csb.numeric import degree2radian, radian2degree

        rand = np.random.random(10000) * 2 * np.pi
        eps = 1e-8
        
        self.assertTrue((np.abs(rand - 
                                degree2radian(radian2degree(rand)))\
                         < eps).all())

    def testDegree2Radian(self):
        from csb.numeric import degree2radian, radian2degree

        rand = np.random.random(10000) * 2 * np.pi
        eps = 1e-8
        
        self.assertTrue((np.abs(rand - 
                                degree2radian(radian2degree(rand)))\
                         < eps).all())


@test.unit
class TestNumeric(test.Case):

    def testTrapezoidal(self):
        from csb.numeric import trapezoidal, exp
        
        x = np.linspace(-10., 10, 1000)
        y = exp(-0.5 * x * x) / np.sqrt(2 * np.pi)

        self.assertAlmostEqual(trapezoidal(x, y), 1.0, 10)

    def testLogTrapezoidal(self):
        from csb.numeric import log_trapezoidal, log
        
        x = np.linspace(-100., 100, 1000)
        y = -0.5 * x * x - log(np.sqrt(2 * np.pi))

        self.assertTrue(abs(log_trapezoidal(y, x)) <= 1e-8)

    def testTrapezoidal2D(self):
        from csb.numeric import trapezoidal_2d, exp
        from numpy import pi
        
        xx = np.linspace(-10., 10, 500)
        yy = np.linspace(-10., 10, 500)

        X, Y = np.meshgrid(xx, yy)
        x = np.array(list(zip(np.ravel(X), np.ravel(Y))))        

        # mean = np.zeros((2,))
        cov = np.eye(2)
        mu = np.ones(2)
        # D = 2
        q = np.sqrt(np.clip(np.sum((x - mu) * np.dot(x - mu, np.linalg.inv(cov).T), -1), 0., 1e308))
        f = exp(-0.5 * q ** 2) / ((2 * pi) * np.sqrt(np.abs(np.linalg.det(cov))))
        f = f.reshape((len(xx), len(yy)))
        I = trapezoidal_2d(f) * (xx[1] - xx[0]) * (yy[1] - yy[0])

        self.assertTrue(abs(I - 1.) <= 1e-8)

    def testSimpson2D(self):
        from csb.numeric import simpson_2d, exp
        from numpy import pi

        xx = np.linspace(-10., 10, 500)
        yy = np.linspace(-10., 10, 500)

        X, Y = np.meshgrid(xx, yy)
        x = np.array(list(zip(np.ravel(X), np.ravel(Y))))        

        # mean = np.zeros((2,))
        cov = np.eye(2)
        mu = np.ones(2)
        # D = 2
        q = np.sqrt(np.clip(np.sum((x - mu) * np.dot(x - mu, np.linalg.inv(cov).T), -1), 0., 1e308))
        f = exp(-0.5 * q ** 2) / ((2 * pi) * np.sqrt(np.abs(np.linalg.det(cov))))
        f = f.reshape((len(xx), len(yy)))
        I = simpson_2d(f) * (xx[1] - xx[0]) * (yy[1] - yy[0])

        self.assertTrue(abs(I - 1.) <= 1e-8)

    def testLogTrapezoidal2D(self):
        from csb.numeric import log_trapezoidal_2d, log
        from numpy import pi

        xx = np.linspace(-10., 10, 500)
        yy = np.linspace(-10., 10, 500)

        X, Y = np.meshgrid(xx, yy)
        x = np.array(list(zip(np.ravel(X), np.ravel(Y))))        

        # mean = np.zeros((2,))
        cov = np.eye(2)
        mu = np.ones(2)
        # D = 2
        q = np.sqrt(np.clip(np.sum((x - mu) * np.dot(x - mu, np.linalg.inv(cov).T), -1), 0., 1e308))
        f = -0.5 * q ** 2 - log((2 * pi) * np.sqrt(np.abs(np.linalg.det(cov))))
        f = f.reshape((len(xx), len(yy)))

        logI = log_trapezoidal_2d(f, xx, yy)

        self.assertTrue(abs(logI) <= 1e-8)
        
    
if __name__ == '__main__':
    test.Console()

