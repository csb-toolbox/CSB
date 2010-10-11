import csb.test as test
import numpy as np

@test.unit
class TestNumeric(test.Case):

    def testTrapezoidal(self):
        from csb.numeric import trapezoidal
        from csb.math import exp
        x = np.linspace(-10.,10,1000)
        y = exp(-0.5 * x * x)/np.sqrt(2 * np.pi)

        self.assertEqual(trapezoidal(x,y),1.0)

    def testLogTrapezoidal(self):
        from csb.numeric import log_trapezoidal
        from csb.math import log
        x = np.linspace(-100.,100,1000)
        y = -0.5 * x * x - log(np.sqrt(2 * np.pi))

        self.assertTrue(abs(log_trapezoidal(y,x))<= 1e-8)

    def testTrapezoidal2D(self):
        from csb.numeric import trapezoidal_2d
        from csb.math import exp
        from numpy import pi
        
        xx = np.linspace(-10.,10,500)
        yy = np.linspace(-10.,10,500)

        X,Y = np.meshgrid(xx,yy)
        x = np.array(zip(np.ravel(X),np.ravel(Y)))        

        mean = np.zeros((2,))
        cov = np.eye(2)
        mu = np.ones(2)
        D = 2
        q = np.sqrt(np.clip(np.sum((x-mu)*np.dot(x-mu,np.linalg.inv(cov).T),-1),0.,1e308))
        f =  exp(- 0.5 * q**2) /((2*pi)* np.sqrt(np.abs(np.linalg.det(cov))))
        f = f.reshape((len(xx),len(yy)))
        I = trapezoidal_2d(f)* (xx[1] - xx[0]) * (yy[1] - yy[0])

        self.assertTrue(abs( I - 1.)<= 1e-8)



    def testSimpson2D(self):
        from csb.numeric import simpson_2d
        from csb.math import exp
        from numpy import pi

        xx = np.linspace(-10.,10,500)
        yy = np.linspace(-10.,10,500)

        X,Y = np.meshgrid(xx,yy)
        x = np.array(zip(np.ravel(X),np.ravel(Y)))        

        mean = np.zeros((2,))
        cov = np.eye(2)
        mu = np.ones(2)
        D = 2
        q = np.sqrt(np.clip(np.sum((x-mu)*np.dot(x-mu,np.linalg.inv(cov).T),-1),0.,1e308))
        f =  exp(- 0.5 * q**2) /((2*pi)* np.sqrt(np.abs(np.linalg.det(cov))))
        f = f.reshape((len(xx),len(yy)))
        I = simpson_2d(f)* (xx[1] - xx[0]) * (yy[1] - yy[0])

        self.assertTrue(abs( I - 1.)<= 1e-8)


    def testLogTrapezoidal2D(self):
        from csb.numeric import log_trapezoidal_2d
        from csb.math import log
        from numpy import pi

        xx = np.linspace(-10.,10,500)
        yy = np.linspace(-10.,10,500)

        X,Y = np.meshgrid(xx,yy)
        x = np.array(zip(np.ravel(X),np.ravel(Y)))        

        mean = np.zeros((2,))
        cov = np.eye(2)
        mu = np.ones(2)
        D = 2
        q = np.sqrt(np.clip(np.sum((x-mu)*np.dot(x-mu,np.linalg.inv(cov).T),-1),0.,1e308))
        f =  - 0.5 * q**2  - log((2*pi)* np.sqrt(np.abs(np.linalg.det(cov))))
        f = f.reshape((len(xx),len(yy)))

        logI = log_trapezoidal_2d(f,xx,yy)

        self.assertTrue(abs( logI )<= 1e-8)
        
    
if __name__ == '__main__':
    test.Console()

