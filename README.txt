Computational Structural Biology Toolbox
========================================

CSB is a python library and application framework, which can be used
to solve problems in the field of structural bioinformatics. If
you are a bioinformatician, software engineer or a researcher working
in this field, chances are you may find something useful here. Our
package consists of a few major components:

    1. Core class library -- object-oriented, granular, with an emphasis
       on design and clean interfaces. 
       
    2. Application framework -- console applications ("protocols"),
       which consume objects from the core library in order to build
       something executable (and hopefully useful).
       
    3. Test framework -- ensures that the library *actually* works.


Compatibility
-------------

In short: CSB requires python 2.6 or higher.

CSB is being developed on Linux, under python 2.7. However, compatibility
is a design goal and the package works on any platform, on any modern python
interpreter since version 2.6 (that includes python 3 support out of
the box). If you found any issues on a platform/interpreter different from
our development environment, please let us know.


Installation 
------------

Full installation instructions can be found in the INSTALL file packaged with
this release, as well as on the project's web site:

    http://csb.codeplex.com/documentation

Here we provide only a brief summary of the installation procedure.
First, make sure all required dependencies are installed:
    
    1. numpy, scipy -- required
    2. matplotlib and wxPython -- optional, needed only if you want to use csb.io.plots
    3. unittest2, argparse -- needed only if you are running python 2.6
    
Next, install CSB::     

    $ python setup.py install
    
CSB has just been installed at the following location::

    $ python -c "import csb, os; print(os.path.abspath(os.path.join(csb.__path__[0], '..')))"

Let us call this path *$LIB*.
            

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
unit tests, please refer to the API documentation, package "csb.test".


Running CSB Applications
------------------------

CSB is bundled with a number of executable console csb.apps. Each app
provides a standard command line interface. To run any app, try::

    $ python $LIB/csb/apps/appname.py --help
    
where *appname* is the name of the application. For more details on
our app framework, including guidelines for writing new applications,
please refer to the API documentation, package "csb.apps".


Documentation
-------------

The project's web site at `CodePlex <http://csb.codeplex.com>`_ contains
online documentation and samples. Be sure to check
    
    http://csb.codeplex.com/documentation

Detailed API documentation can be found in the "docs/api" directory in the
distribution package (docs/api/index.html). Many packages contain
introductory module level documentation and samples/tutorials. These are also
available in the HTML docs, but a quick way to access them is by using
the built-in python help system. For example, for a general introduction
see the module documentation of the root package::

    $ python -c "import csb; help(csb)"

If you are interested in a specific package, such as cs.bio.sequence,
try::    
    
    $ python -c "import csb.bio.sequence; help(csb.bio.sequence)"


Contact
-------

CSB is developed by Michael Habeck's Computational Structural Biology
`research group <http://www.eb.tuebingen.mpg.de/research/research-groups/michael-habeck.html>`_.
    
For complete source code, contribution, support or bug reports please visit
our web site at CodePlex:
  
    `csb.codeplex.com <http://csb.codeplex.com>`_
    

License
-------

CSB is open source and distributed under OSI-approved MIT license.
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

