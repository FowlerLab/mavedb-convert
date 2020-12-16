[![Build Status](https://travis-ci.com/VariantEffect/mavedbconvert.svg?branch=main)](https://travis-ci.com/VariantEffect/mavedbconvert)
[![Coverage Status](https://coveralls.io/repos/github/VariantEffect/mavedbconvert/badge.svg?branch=main)](https://coveralls.io/github/VariantEffect/mavedbconvert?branch=main)
[![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black)

# mavedbconvert
A command line tool for converting Multiplex Assay of Variant Effect datasets into a MaveDB-ready format.

## Important note
This tool may generate input files that are invalid under the forthcoming MAVE-HGVS format,
and therefore will not be suitable for upload to MaveDB.

In particular, variants should be in sorted order and variant strings containing ambiguous symbols such as 
"?" "N" or "Xaa" will not be accepted.

# Installation
Download the mavedbconvert source and navigate to that directory.
We recommend creating a [virtual environment](https://docs.python.org/3/library/venv.html) before proceeding with the installation.

Install the package using pip:

    pip3 install .

## Python 3.9 compatibility
The tests are currently failing under Python 3.9 due to multiple HDF5 files being open.
This is being worked on in an updated version.
If you encounter this error using mavedbconvert for normal use, please open a GitHub issue.