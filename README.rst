Computational Structural Biology Toolbox
========================================

CSB is a python library and application framework, which can be used
to solve problems in the field of structural bioinformatics. If
you are a bioinformatician, software engineer or a researcher working
in this field, chances are you may find something useful here. Our
package consists of a few major components:

1. Core class library - object-oriented, granular, with an emphasis
   on design and clean interfaces.

2. Application framework - console applications ("protocols"),
   which consume objects from the core library in order to build
   something executable (and hopefully useful).

3. Test framework - ensures that the library *actually* works.


Installation 
------------
CSB is being developed on Linux. However, compatibility
is a design goal and the package works on any platform, on any modern Python
interpreter. If you find any issues on a platform/interpreter different from
our development environment, please let us know.

CSB and all of its dependencies can be installed with pip::

    $ pip install csb

See http://csb-toolbox.github.io/installation for more details.


Running CSB Applications
------------------------

CSB is bundled with a number of executable console csb.apps. Each app
provides a standard command line interface. To run any app, try::

    $ csb-app --help
    
where *csb-app* is the name of the application, such as ``csb-hhfrag``.
For more details on our app framework, including guidelines for writing new
applications, please refer to the API documentation, package "csb.apps".


Documentation
-------------

The project's web site at `GitHub <http://github.com/csb-toolbox>`_ contains
online documentation and samples. Visit us at:
    
http://csb-toolbox.github.io

Detailed API documentation can be found in the "docs/api" directory in the
distribution package (docs/api/index.html). This documentaiton is also hosted
on our web site:

https://csb-toolbox.github.io/api-docs/

Many packages contain introductory module level documentation and samples/tutorials.
These are also available in the HTML docs, but a quick way to access them is by using
the built-in python help system. For example, for a general introduction
see the module documentation of the root package::

    $ python -c "import csb; help(csb)"

If you are interested in a specific package, such as cs.bio.sequence,
try::    
    
    $ python -c "import csb.bio.sequence; help(csb.bio.sequence)"


Contact
-------

CSB is developed by Michael Habeck's Computational Structural Biology
`research group <http://www.stochastik.math.uni-goettingen.de/index.php?id=172>`_.
    
For complete source code, contributions, support or bug reports please visit
us on GitHub:
  
http://github.com/csb-toolbox/
    

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
    
------------

.. image:: https://travis-ci.org/csb-toolbox/CSB.svg?branch=master
   :target: https://travis-ci.org/csb-toolbox

