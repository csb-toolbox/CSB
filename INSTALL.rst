

System Requirements
-------------------

CSB is distributed as a pure python package, which renders it cross-platform.
Installation is therefore a very trivial task. The source code is also kept
compatible across python versions, so we guarantee that every official release
package will work on both python 2 and 3 without any change. We test all nightly
builds against python versions 2.7 and 3.2, but technically CSB supports python 2.6
and all subsequent versions except probably 3.0.


User Installation
-----------------

Install csb from PyPi::

    $ pip install csb

pip will take care of all necessary dependencies. However, on Windows this might
fail. If this is the case, install numpy and scipy manually and rerun the pip
installation shown above. The easiest way to get numpy and scipy on Windows is 
to install Anaconda::

    $ conda install numpy scipy
    $ pip install csb


Developer Installation
----------------------

ABOUT DEPENDENCIES

CSB depends on 2 very well-known and readily available python packages:

* numpy - required (numpy.scipy.org)
* scipy - required (scipy.org)
* matplotlib and wxPython - optional, needed only if you want to use csb.io.plots

On python 2.6 you will need these two as well:

* unittest2
* argparse

both of which are standard, simple and available for download from PyPi (pypi.python.org).

Epydoc is required to build the API documentation.


INSTALLATION

First get the source code from github::

    $ git clone https://github.com/csb-toolbox/CSB.git
    $ cd CSB

Then install the project in editable mode::

    $ pip install --editable .[dev]

pip will take care of all runtime and development dependencies and make the source
code location importable. No need to modify your $PYTHONPATH.

If dependencies cannot be resolved on your system - this is often the case on
Windows - the best way to go from here is to install Anaconda. Then use the conda
package manager to install all missing packages::

    $ conda install --file requirements.txt
    $ pip install --editable .[dev]



Testing
-------

Running the CSB test suite may be useful in order to check if your installation works.
All CSB tests are executed with the csb.test.Console. A typical way to run the console is::

    $ csb-test "csb.test.cases.*"


API Documentation
-----------------

CSB comes with API docs in HTML format. Simply navigate to the docs/api folder in the
release package and open index.html with any web browser. Note that the docs are not
installed, so you need to keep your own copy.

