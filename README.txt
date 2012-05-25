Computational Structural Biology Toolbox
========================================

CSB is a python library and application framework, which can be used
to solve problems in the field of Computational Structural Biology. If
you are a bioinformatician, software engineer or a researcher working
in this field, chances are you may find something useful here. Our
package consists of a few major components:

    1. Core class library -- object-oriented, granular, with an emphasis
	   on design and clean interfaces. 
	   
	2. Application framework -- console applications ("protocols"),
	   which consume objects from the core library in order to build
	   something executable (and hopefully useful).
	   
	3. Test framework -- makes sure the library *actually* works.


Compatibility
-------------

In short: requires python 2.6 or higher.

CSB is being developed on Linux, under python 2.7. However, compatibility
is a design goal and the package works on any platform, on any modern python
interpreter since version 2.6 (yes, that includes python 3 support out of
the box). If you found any issues on a platform/interpreter different from
our development environment, please let us know. We will kindly apologize,
fix the problem and write a regression test case to make sure that the
problem will be solved once and forever.


Installation 
------------

Is your interpreter *Compatible*? If yes, proceed with the installation
of a few trivial packages, on which CSB depends:

    1. numpy -- mandatory (numpy.scipy.org)
    2. scipy -- mandatory (scipy.org)

    3. matplotlib and wxPython -- optional, needed only if you want to
       use csb.io.plotting
       
    4. unittest2 -- needed only if you are using python 2.6
    5. argparse -- may be needed if you are using python 2.6

If you are installing CSB on Windows, numpy, scipy and matplotlib can be
installed from binary MSI/exe packages, the rest can be installed from
PyPi. On Debian-like Linux distributions use your package manager::

	$ apt-get install *package*

where *package* is one of: python-numpy, python-scipy, python-matplotlib.

Finally, install CSB itself::

	$ python setup.py install
	
CSB has been just installed at the following location::

	$ python -c "import csb, os; print(os.path.abspath(os.path.join(csb.__path__[0], '..')))"

Let's call this path *$LIB*.
			

Testing
-------

Running the CSB test suite may be useful if you made any modifications to
the source code, or if you simply want to check if your installation works.

All CSB tests are executed with the csb.test.Console. A typical way to run 
the console is::

	$ python $LIB/csb/test/app.py "csb.test.cases.*"
	
or just::

	$ python $LIB/csb/test/app.py 		

For help try::

	$ python $LIB/csb/test/app.py -h	

For more details on our test framework, including guidelines for writing
unit test, please refer to the API documentation, package "csb.test".


Running CSB Applications
------------------------

CSB is bundled with a number of executable console csb.apps. Each app
provides a standard command line interface. To run any app, try::

	$ python $LIB/csb/apps/*appname*.py --help
	
where *appname* is the name of the application.	For more details on
our app framework, including guidelines for writing new applications,
please refer to the API documentation, package "csb.apps".


Documentation
-------------

A quick introduction/tutorial is available in the root package "csb".
Try::

	$ python -c "import csb; help(csb)"

Detailed API documentation can be found in the "docs/api" directory in the
distribution package. Simply open docs/api/index.html in your browser.


Contact
-------

This library is developed by Michael Habeck's Computational Structural
Biology Group, at the Max Planck Institute for Developmental Biology,
Tuebingen, Germany:

	http://www.eb.tuebingen.mpg.de/research/research-groups/michael-habeck.html
	
The source code of the library is hosted at the following SVN repository:

	https://svn.tuebingen.mpg.de/agbs/projects/CSB
		
Send any requests and feedback to:
	
	michael.habeck@tuebingen.mpg.de


License
-------

CSB is an open source library, distributed under OSI-approved MIT license.
::

    Copyright (c) 2012 Michael Habeck
    
    Permission is hereby granted, free of charge, to any person obtaining
    a copy of this software and associated documentation files (the
    "Software"), to deal in the Software without restriction, including
    without limitation the rights to use, copy, modify, merge, publish,
    distribute, sublicense, and/or sell copies of the Software, and to
    permit persons to whom the Software is furnished to do so, subject to
    the following conditions:
    
    The above copyright notice and this permission notice shall be
    included in all copies or substantial portions of the Software.
    
    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
    EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
    MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
    IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY
    CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
    TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
    SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

