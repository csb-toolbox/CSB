## Introduction

CSB is a granular, object-oriented bioinformatics library, which provides 
abstractions you can use to build your own bioinformatics applications. 
We put a strong emphasis on ease of use and convenience, but maybe our 
definition of ease of use is a bit different from what is frequently 
expected from a python / bioinformatics library. We strongly believe that 
an API is easy to use when it is designed with the following properties in mind:

* high architectual and coding standards
* reusability and extensibility
* clean, abstracted and well-encapsulated interfaces
* strictly consistent naming conventions in accord with the industry standards
* use of classic design patterns
* preference of good/obvious design over extensive code comments and documentation

In other words, if you are an experienced developer and you don't like sloppy 
programming, you should feel right at home.

Here we present a brief tutorial on some of the most fundamental parts of 
CSB. It will help you to get started and explore the scope of the library. 
After going through these examples you will get enough confidence in order 
to start experimenting with the APIs and reading the comprehensive 
[API docs](http://pythonhosted.org/csb/), packaged with every release.

## Overview

The library is composed of a set of highly branched python packages 
(namespaces). Some of the packages are meant to be directly used by the 
clients (core library), while others are utility modules and take part 
in the development of the library:

* Core class library
* Application framework -- executable console applications ("protocols"), 
which consume objects from the core library. The framework ensures that each 
CSB application is also reusable and can be instantiated as a regular python 
object without any ugly side effects (sys.exit() and friends). See ``csb.apps`` 
in our [API docs](http://pythonhosted.org/csb/) for more details.
* Test framework -- built on top of the standard unittest as a thin wrapping 
layer. Provides some sugar like transparent management of test data files, 
and modular test execution. ``csb.test`` will give you all the details.

The core library is roughly composed of:
* bioinformatics API: ``csb.bio``, which includes stuff like ``csb.bio.io``, 
``csb.bio.structure``, ``csb.bio.sequence``, ``csb.bio.hmm``
* statistics API: ``csb.statistics``, ``csb.numeric``
* utilities - ``csb.io``, ``csb.core`

## Tutorials

* [Basic Tutorial](Basic Tutorial.md)
* [Structures and PDB I/O](Structures and PDB IO.md)
* [Sequences, Alignments and FASTA I/O](Sequence Alignments and FASTA IO.md)
* [Profile HMMs and HHpred I/O](Profile HMMs and HHpred IO.md)
* [Probability Distributions](Probability Distributions.md)
* [Monte Carlo Sampling](Monte Carlo Sampling.md)

## Installation

Detailed instructions are provided [here](Installation.md).

## Development

Documentation for developers can be found [here](For Developers.md).

## More Documentation

Be sure to check our tutorials and ultimately our comprehensive 
[API docs](http://pythonhosted.org/csb/), bundled with each release package.

