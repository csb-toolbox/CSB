## Manual Builds

The build and test processes are implemented in CSB itself. Use ``csb.build`` 
to build, test and package the entire project:


    $ python csb/build.py -o <output directory>

This will produce a standard installable release package. If Epydoc is installed, 
our HTML API docs will be compiled and included as well.

## Build and Test Automation

Under a Continuous Integration model (CI), CSB should be built and tested nightly 
(or multiple times per day, or even after each commit). Given the current workload, 
nightly builds are sufficient. Nightly builds can be configured using a simple 
scheduled (cron) job that runs the ``csb/build.py`` program at regular intervals 
using arbitrary Python interpreter versions (Python 2.6 and 3.2 by default). To 
automate this task, we provide a shell script 
([build.sh](build.sh)), which performs the following steps:

* Checkout the project from the central repository at GitHub
* Run ``csb/build.py`` with the desired Python interpreter and collect all build 
  artifacts and API docs
* Collect and analyze the output build log (success, failure, crash)
* Notify the specified operators (send email)

#### Configuration

* Install and configure a local SMTP server, e.g. postfix: 

    $ sudo apt-get install postfix

Make sure it works:

    
    $ python -c "import smtplib; smtplib.SMTP('localhost').sendmail('your@email.here', 'your@email.here', 'test')"

If you can't install or use a local SMTP, just register a new gmail account, e.g. _csb.build@gmail.com_.

* Open [build.sh](build.sh) and set the build directory ($ROOT)
* Configure also $BOTMAIL (the email address of the build bot, e.g. _csb.build@gmail.com_), 
  $SMTP (smtp host, e.g. smtp.gmail.com) and $OPERATORS (list of recipients)
* If $SMTP requires authentication, configure also $BOTPWD; otherwise, leave it blank
* Save your changes

#### Nightly Builds

* Configure scheduled tasks for automated builds with Python 2 and Python 3:

     0 0 * * * bash -lc '/path/to/build.sh'         # Python 2 (default)
    30 0 * * * bash -lc '/path/to/build.sh python3' # Python 3

The optional parameter of build.sh is the python interpreter to use.  You could include 
additional tasks for arbitrary platforms (e.g. python2.6, python3.3), assuming that such 
python executables exist on your machine.

#### Build Artifacts

The build script will produce the following artifacts:

* $ROOT/pythonx.y/CSB/log --- build log
* $ROOT/pythonx.y/CSB/build.zip --- wraps the final release tarball (csb-x.y.z.tar.gz)
* $ROOT/pythonx.y/CSB/csb-x.y.z.tar.gz --- release tarball 
* $ROOT/pythonx.y/CSB/build/docs/api --- API docs

You can expose them using a local web server for easy access:

    Alias /csb /path/to/BUILD/python2/CSB
    Alias /csb3 /path/to/BUILD/python3/CSB

and then bookmark the following links:

    localhost/csb/log
    localhost/csb3/log

    localhost/csb/build.zip
    localhost/csb3/build.zip

    localhost/csb/build/docs/api/index.html

## Release Management

**Note**: When testing a release package, make sure it does not interfere 
with other releases that might be installed on the same machine. This is 
achieved by manipulating the $PYTHONPATH:

    
    $ tar -xzvf csb-x.y.z.tar.gz
    $ mv csb-x.y.z release   # dots in directory names cause problems with python
    $ export PYTHONPATH=release
    $ python release/csb/test/app.py

Follow these steps to build and deploy a new release:

1. Tag and build the release
	* Run the entire test suite on as many platforms as possible (e.g. Linux, Mac, Windows) 
	  and make sure there are no failures
	* Run the entire test suite with different python versions on at least one platform (e.g. Linux)
	* Bump the version number stored in ``csb.__version__`` and commit
	* Add a new repository tag matching the new version: ``R-x.y.z``; commit
	* Trigger a new build with your build.sh
2. Push to CodePlex
	* Login with your CodePlex credentials
	* Go to [Downloads](csb.codeplex.com/releases) and click "New Release"
	* Fill in Name (csb x.y.z) and Release Notes (e.g. change log, fixed issues, etc)
	* Browse the release package (tar.gz) and set Type=Source Code, Status=Stable, Show=True
3. Push to PyPi
	* Login with your PyPi credentials, Gmail or OpenID
	* Go to [PyPi](https://pypi.python.org/pypi?%3Aaction=submit_form)
	* Browse the PKG-INFO file extracted from the release package and click "Add"
	* Click on "files" next to the new version and add the release package (tar.gz) with type=source
	* Extract the API docs from the release package and archive them in a zip file
	* Browse the zip file and click "Upload Documentation"
