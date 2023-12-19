import sys

from setuptools import setup, find_packages
from csb.build import ROOT 
from io import open

try:
    __doc__ = open('README.rst', encoding="utf-8").read()
except IOError:
    __doc__ = ""


NAME = ROOT
AUTHOR = "Michael Habeck et al."
EMAIL = "ivan.kalev@gmail.com"
URL = "http://github.com/csb-toolbox"
SUMMARY = "Computational Structural Biology Toolbox"
DESCRIPTION = __doc__
LICENSE = 'MIT'

REQUIREMENTS = open("requirements.txt", encoding="utf-8").readlines()
DEV_REQUIREMENTS = []

v = {}
exec(open(ROOT + "/__init__.py", encoding="utf-8").read(), v)
VERSION = v["Version"]()


def build():

    return setup(
        name=NAME,
        packages=find_packages(),
        include_package_data=True,
        version=VERSION.short,
        author=AUTHOR,
        author_email=EMAIL,
        url=URL,
        description=SUMMARY,
        long_description=DESCRIPTION,
        license=LICENSE,
        install_requires=REQUIREMENTS,
        tests_require=DEV_REQUIREMENTS,
        extras_require={
            'dev': DEV_REQUIREMENTS
        },
        test_suite="csb.test.cases",
        entry_points={
            'console_scripts': [
                'csb-test = csb.test.app:main',
                'csb-bfit = csb.apps.bfit:main',
                'csb-bfite = csb.apps.bfite:main',
                'csb-csfrag = csb.apps.csfrag:main',
                'csb-hhfrag = csb.apps.hhfrag:main',
                'csb-buildhmm = csb.apps.buildhmm:main',
                'csb-hhsearch = csb.apps.hhsearch:main',
                'csb-precision = csb.apps.precision:main',
                'csb-promix = csb.apps.promix:main',
                'csb-embd = csb.apps.embd:main'
            ]
        },
        classifiers=(
            'Development Status :: 5 - Production/Stable',
            'Intended Audience :: Developers',
            'Intended Audience :: Science/Research',
            'License :: OSI Approved :: MIT License',
            'Operating System :: OS Independent',
            'Programming Language :: Python',
            'Programming Language :: Python :: 3.6',
            'Programming Language :: Python :: 3.7',
            'Programming Language :: Python :: 3.8',
            'Programming Language :: Python :: 3.9',
            'Programming Language :: Python :: 3.10',
            'Programming Language :: Python :: 3.11',
            'Programming Language :: Python :: 3.12',
            'Topic :: Scientific/Engineering',
            'Topic :: Scientific/Engineering :: Bio-Informatics',
            'Topic :: Scientific/Engineering :: Mathematics',
            'Topic :: Scientific/Engineering :: Physics',
            'Topic :: Software Development :: Libraries'
        )
    )


if __name__ == '__main__':
    build()
