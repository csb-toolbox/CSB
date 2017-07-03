## Manual Builds

The build and test processes are implemented in CSB itself. Use ``csb.build`` 
to build, test and package the entire project:


    $ python csb/build.py -o <output directory>

This will produce a standard installable release package. If Epydoc is installed, 
our HTML API docs will be compiled and included as well.


## Build and Test Automation

Under a Continuous Integration model (CI), CSB should be built and tested after each 
commit). CI builds can be configured using a simple 
Travis CI job which runs the ``csb/build.py`` program at regular intervals 
using arbitrary Python interpreter versions (Python 2.7 and 3.6 by default). To 
automate this task, we provide a Travis CI job configuration, which performs the following:

* Checkout the project from the central repository at GitHub
* Run ``csb/build.py`` with the desired Python interpreter and collect all build 
  artifacts and API docs
* Collect and analyze the output build log (success, failure, crash)
* Notify the specified operators (send email)

## Release Management

**Note**: When testing a release package, make sure it does not interfere 
with other releases that might be installed on the same machine. This is 
achieved by installing the package in a virtual environment. Here is an 
example with anaconda:

    $ conda create -n test python=2.7
    $ pip install csb-x.y.z
    $ csb-test

Follow these steps to build and deploy a new release:

1. Tag and build the release
	* Run the entire test suite on as many platforms as possible (e.g. Linux, Mac, Windows) 
	  and make sure there are no failures
	* Run the entire test suite with different python versions on at least one platform (e.g. Linux)
	* Bump the version number stored in ``csb.__version__`` and commit
	* Add a new repository tag matching the new version: ``R-x.y.z``; commit and push
	* Run a new build
2. Push to GitHub
	* Login with your GitHub credentials
	* Go to [Releases](https://github.com/csb-toolbox/CSB/releases) and click "Draft a new release"
	* Fill in the Name (CSB x.y.z) and Release Notes (e.g. change log, fixed issues, etc)
	* Browse and upload the release package (tar.gz) 
3. Push to PyPi
	* Login with your PyPi credentials, Gmail or OpenID
	* Go to [PyPi](https://pypi.python.org/pypi?%3Aaction=submit_form)
	* Browse the PKG-INFO file extracted from the release package and click "Add"
	* Click on "files" next to the new version and add the release package (tar.gz) with type=source
	* Extract the API docs from the release package and archive them in a zip file
	* Browse the zip file and click "Upload Documentation"
