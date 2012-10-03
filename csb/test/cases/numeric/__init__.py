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

@test.functional
class InvertibleMatrixTest(test.Case):

    
    d = 100
    # Find some random invertible matrix
    m_general = np.random.uniform(low=-2, high=2, size=(d, d))
    while np.linalg.det(m_general) == 0:
        m_general = np.random.uniform(low=-2, high=2, size=(d, d))
        
    
    m_general_inv = np.linalg.inv(m_general)
    m_diagonal = np.diag(np.random.uniform(low=-2, high=2, size=d))
    m_diagonal_inv = np.linalg.inv(m_diagonal)
    m_2unity = np.eye(d) * 2.
    m_2unity_inv = np.linalg.inv(m_2unity)

    def invertMatrix(self, matrix):
        from csb.numeric import InvertibleMatrix

        for i in range(1000):
            testmatrix = InvertibleMatrix(matrix)
            dummy = testmatrix.inverse

    @test.skip("machine-specific")
    def testDiagonalTweak(self):
        
        self.assertFasterThan(0.75, self.invertMatrix, self.m_diagonal)

    @test.skip("machine-specific")
    def testUnityMultipleTweak(self):

        self.assertFasterThan(0.5, self.invertMatrix, self.m_2unity)
        
    def testCaching(self):
        from csb.numeric import InvertibleMatrix

        # Initialize with matrix
        testmatrix = InvertibleMatrix(self.m_general)
        # Check if it worked
        self.assertListEqual(testmatrix._matrix.flatten().tolist(), 
                             self.m_general.flatten().tolist())
        # Because the testmatrix.inverse hasn't been accessed yet, testmatrix._inverse should be None
        self.assertEqual(testmatrix._inverse_matrix, None)
        # Now we access testmatrix.inverse, which onyl now actually calculates the inverse
        self.assertListEqual(testmatrix.inverse.flatten().tolist(), 
                             self.m_general_inv.flatten().tolist())

        # Let's change testmatrix via testmatrix.__imul__
        testmatrix *= 2.
        # Check if that worked
        self.assertListEqual(testmatrix._matrix.flatten().tolist(),
                             (2.0 * self.m_general.flatten()).tolist())
        # This operation should not have changed the testmatrix._inverse_matrix field, as
        # we didn't access testmatrix.inverse again
        self.assertListEqual(testmatrix._inverse_matrix.flatten().tolist(), 
                             self.m_general_inv.flatten().tolist())
        # New we access testmatrix.inverse, which calculates the inverse and updates the field
        self.assertListEqual(testmatrix.inverse.flatten().tolist(), 
                             (self.m_general_inv / 2.0).flatten().tolist())

        # The same again for testmatrix.__idiv__
        testmatrix /= 2.
        # Check if that worked
        self.assertListEqual(testmatrix._matrix.flatten().tolist(),
                             (self.m_general.flatten()).tolist())
        # This operation should not have changed the testmatrix._inverse_matrix field, as
        # we didn't access testmatrix.inverse again
        self.assertListEqual(testmatrix._inverse_matrix.flatten().tolist(), 
                             (self.m_general_inv / 2.0).flatten().tolist())
        # New we access testmatrix.inverse, which calculates the inverse and updates the field
        self.assertListEqual(testmatrix.inverse.flatten().tolist(), 
                             self.m_general_inv.flatten().tolist())
        
        # Initialize with inverse matrix
        testmatrix = InvertibleMatrix(inverse_matrix=self.m_general_inv)
        # Let's see if that worked, e.g. if the testmatrix._inverse_matrix field has been
        # set correctly
        self.assertListEqual(testmatrix._inverse_matrix.flatten().tolist(), 
                             self.m_general_inv.flatten().tolist())
        # Check if the property returns what it's supposed to be
        self.assertListEqual(testmatrix.inverse.flatten().tolist(), 
                             self.m_general_inv.flatten().tolist())
        # We didn't call testmatrix.__getitem__() yet, so testmatrix._matrix should be None
        self.assertEqual(testmatrix._matrix, None)
        # To include the numerical error
        invinv = np.linalg.inv(self.m_general_inv)
        # Now we access testmatrix by its __getitem__ method, which calculates the
        # testmatrix._matrix field from the testmatrix._inverse_matrix by inversion
        for i in range(len(testmatrix)):
            self.assertListEqual(testmatrix[i].tolist(), invinv[i].tolist())

        testmatrix = InvertibleMatrix(inverse_matrix=self.m_general_inv)
        # Let's change testmatrix via testmatrix.__imul__
        testmatrix *= 2.
        # That shouldn't have changed the testmatrix._matrix field (which currently
        # should be None), but the testmatrix._inverse_matrix field by a factor of 1/2.0 = 0.5
        self.assertEqual(testmatrix._matrix, None)
        self.assertListEqual(testmatrix._inverse_matrix.flatten().tolist(),
                             (self.m_general_inv / 2.0).flatten().tolist())
        # Now we access testmatrix by __getitem__, which calculates the matrix
        # from the inverse and updates the field testmatrix._matrix
        invinv *= 2.0
        for i in range(len(testmatrix)):
            self.assertListEqual(testmatrix[i].tolist(), invinv[i].tolist())

        # The same again for testmatrix.__idiv__
        testmatrix = InvertibleMatrix(inverse_matrix=self.m_general_inv)
        testmatrix /= 2.
        # Check if testmatrix._matrix is None and if the testmatrix._inverse field
        # has been multiplied by a factor of 2.0
        self.assertEqual(testmatrix._matrix, None)
        self.assertListEqual(testmatrix.inverse.flatten().tolist(), 
                             (self.m_general_inv * 2.0).flatten().tolist())
        # All that is supposed to leave testmatrix._matrix with None:
        self.assertEqual(testmatrix._matrix, None)
        # Now we access testmatrix by __getitem__ again, which calculates the matrix from
        # its inverse and updates the field
        invinv /= 4.0
        for i in range(len(testmatrix)):
            self.assertListEqual(testmatrix[i].tolist(), invinv[i].tolist())
        
        
if __name__ == '__main__':
    test.Console()

