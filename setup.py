from __future__ import print_function
import setuptools
from distutils.core import setup
from distutils.extension import Extension
import numpy as np
from distutils.ccompiler import new_compiler
import os
import sys
import tempfile
def main():
    setup(name="crossmatchK2",
          version="0.0.1",
          description="Internal crossmatch functions for K2 and EDR3 Binaries",
          author="Sam Christian",
		  url='https://github.com/Sam-2727/crossmatchK2',
		  packages =['crossmatch'],
		  zip_safe=False,
		  classifiers=[
			  "Intended Audience :: Science/Research",
			  "License :: OSI Approved :: MIT License",
			  "Programming Language :: Python :: 3 :: Only",
		  ],
          author_email="samchristian@mit.edu",
		  licence="MIT")

if __name__ == "__main__":
    main()