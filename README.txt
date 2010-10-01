CSB is a Python 2.6+ class library for computational structural biology.
Developed by the Computational Structural Biology group, Max Planck Institute
for Developmental Biology, Tuebingen, Germany.



Installation

1. Install the latest NumPy distribution - http://numpy.scipy.org
 
2. Install CSB:

   a. Standard way
   
      Unpack the archive. The standard way to install it is simply to run
      the setup.py file in this distribution:

          python setup.py install
	
   b. Custom installation
 
      After installing NumPy there are no more dependencies, so custom 
      installations are straightforward.
 	   
      After unpacking the archive you need to let your Python installation
      know about the "csb" package. This is done by pointing your PYTHONPATH
      environment variable to the directory where "csb" is located.  		

3. Test

Running the CSB test suite may be useful if you make any modifications to
the source code, or if you simply want to check if your installation works.

All CSB tests are executed with the CSB Test Console. A typical way to run 
the console is:	  

	python csb/test/app.py "csb.test.cases.*"	

For help try:

	python csb/test/app.py -h



Documentation

Epydoc API documentation can be found in the "docs/api" directory in the
distribution package. Simply open docs/api/index.html in your browser.



Contact

michael.habeck@tuebingen.mpg.de
