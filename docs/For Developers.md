## Getting the Source

You will need to install a ``git`` client first.

#### On the command line
Once you have git, checkout the entire project:

    $ git clone https://github.com/csb-toolbox/CSB.git

This will create a project called "CSB" in the current directory. 
For committing any changes, you will also need to provide your 
GitHub username and password or setup SSH authentication.

Next, install all dependencies and the package in editable mode:

    $ pip install --editable .[dev]

#### Using PyCharm

Open the ``CSB`` folder in PyCharm. Don't commit any ``.idea`` files.

## Building

Use ``csb.build`` to build, test and package the entire project:


    $ python csb/build.py -o <output directory>
    
Full details are provided [here](Build Automation and Deployment).

## Testing

You should run the complete test suite before each commit. All CSB 
tests are executed with the ``csb.test.Console``. To execute all tests:


    $ csb-test "csb.test.cases.*"

or just:


    $ csb-test         

To execute test suites for specific modules only:


    $ csb-test "csb.test.cases.module1.*" "csb.test.cases.module2"...

Run with "-h" for help. For more details on our test framework, 
including guidelines for writing unit test, please refer to the 
[API documentation](http://pythonhosted.org/csb/), package 
[csb.test](http://pythonhosted.org/csb/csb.test-module.html).

You should also test the syntax of your docstrings before each commit:


    $ epydoc --html -v -o /tmp --name CSB --fail-on-error --fail-on-warning --fail-on-docstring-warning csb

For convenience, you can run [test.sh](test.sh) before pushing code to 
the repository to make sure that all test cases and docstring checks pass.

## Methodology

We have embraced a Test-driven development (TDD) approach. That is, all 
new features and bug fixes must be accompanied by corresponding test 
cases in a subpackage of ``csb.test.cases``. For more details on this 
opic please read the [API documentation](http://pythonhosted.org/csb/), 
package [csb.test](http://pythonhosted.org/csb/csb.test-module.html).

## Style Guide

#### Compatibility

* make sure your code is cross-platform (Windows, Linux, Mac)
* make sure your code is cross-interpreter compatible (python 2.7+:  
2.7, 3.5, 3.6, ...)
To ensure cross-interpreter compatibility, please write Python 3 
compatible code. In the rare event that a given construct in Python 3 
is not backwards compatible with Python 2.7, please ensure portability 
by following these guidelines:

* you shouldn't feel the need to use ``print`` in CSB (even in apps); 
  but if really necessary, import and use the print function (``print()``)
* use ``csb.io.Pickle``, ``csb.io.MemoryStream`` and ``csb.io.urllib`` instead of importing the standard library modules (details in [csb.io](http://pythonhosted.org/csb/csb.io-module.html))
* use ``csb.core.string`` for string ``isinstance`` checks and ``csb.core.metaclass()`` to define metaclass inheritance (details in [csb.core](http://pythonhosted.org/csb/csb.core-module.html))

#### Commit Messages

You should always provide a message for each commit. Commit comments should follow the general Mercurial guidelines:

* the first line should be a complete, independent sentence, which summarizes all changes
* next lines should give as much detail about affected modules as possible
* if resolving an issue, include the issue number
Here is a sample commit comment:

``
Fixed issue #00010: SomeObject crashes with an exception for structures with no title

csb.bio.structure:
 - added a new Structure.xxxx property

csb.bio.io.wwpdb:
 - SomeObject will now check for empty titles

csb.test.cases.bio.io.wwpdb:
 - added regression tests	
``

#### Coding Style, Commenting and Formatting

Format your code using 4-spaces indentation (not tabs). In Eclipse just use {"Ctrl+Shift+F"}.
Your code should follow the PEP8 style guide strictly. Here is a summary:

* ``ClassName``
* ``method_name or function_name``; the only exception are test methods: ``testMethodName``
* ``variable_name``
* ``property_name``
* ``CONSTANT_NAME``
Additional requirements:
* ``_internal_method, _internal_field`` -- unencapsulated classes are not allowed (except if defining a simple struct); mark private methods and fields with a single or double underscore
* ``@property, @property.setter`` -- getters and setters are not allowed; use @property to expose a private field
* ``str.format()`` -- avoid using the deprecated "%" operator
* ``super(ChildClass, self).method()`` -- do not use the deprecated ``ParentClass.method()``
* always define all instance variables in the constructor of a class
* always derive classes from ``object`` explicitly
* define abstract classes with ``ABCMeta, @abstractmethod, @abstractproperty``
* organize imports at the top of each module; avoid private imports
* provide an Epydoc docstring for each public class, method, property
* provide an Epydoc docstring for each module
* prefer short, clean, obvious code and design patterns over writing comments
* always prefix static fields and methods with the corresponding class name; do not access static fields/methods with ``self`` or hide static fields/methods within derived classes

Here is a code template:

```python
"""
Module short summary.

Module multi-line description...
"""

import module1
import module2

from module3 import ObjectOne, ObjectTwo
from module4 import ObjectThree


class ClassName(object):
    """
    Class description...
    ... ends here.

    @param param_name: description
    @type param_name: L{SomeObject} or built-in (str, int, etc)
    """

    __metaclass__ = ABCMeta

    CONSTANT_NAME = value

    def __init__(self, param_name):
        self._private_field = None

    @property
    def property_name(self):
        """
        Description...
        """
        return self._private_field
    @property_name.setter
    def property_name(self, value):
        # validate value here...
        self._property_name = value

    def method_name(self, param_name)
        """
        Description...
        
        @param param_name: description
        @type param_name: L{SomeObject} or built-in (str, int, etc)

        @return: description of @rtype
        @rtype: L{SomeObject} or built-in (str, int, etc)
        
        @raise ExceptionName: condition
        """
        (body)

    @abstractmethod
    def abstract_method_name(self, param_name):
        """
        (docstring)
        """
        pass

    @staticmethod
    def static_method_name(param_name):
        """
        (docstring)
        """
        (body)
```
