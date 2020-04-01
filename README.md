[![Build Status](https://travis-ci.com/VariantEffect/mavedbconvert.svg?branch=master)](https://travis-ci.com/VariantEffect/mavedbconvert)
[![Coverage Status](https://coveralls.io/repos/github/VariantEffect/mavedbconvert/badge.svg?branch=master)](https://coveralls.io/github/VariantEffect/mavedbconvert?branch=master)
[![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black)

# mavedbconvert
A command line tool for converting Multiplex Assay of Variant Effect datasets into a MaveDB-ready format.

# Installation
Download the mavedbconvert source and navigate to that directory.
We recommend creating a [virtual environment](https://docs.python.org/3/library/venv.html) before proceeding with the installation.

Install the package using pip:

    pip3 install .

## Troubleshooting
If you are a OSX user, you may experience header related issues when installing pysam. The current workaround 
is to install pysam v0.13 manually before installing the requirements:

    pip3 install pysam==0.13

This is the latest version known to compile without errors.

Although pysam is not required for mavedbconvert directly, it is installed by some of our dependencies.
Until it is removed or made optional by those libraries, mavedbconvert will unfortunately not be installable on Windows.
