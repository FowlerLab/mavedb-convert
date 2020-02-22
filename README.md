# mavedb-convert
A command line tool for converting alternate file formats into a MaveDB compliant format.

# Installation
Download the `mavedb-convert` source and navigate to that directory.
We recommend creating a [virtual environment](https://docs.python.org/3/library/venv.html) before proceeding with the installation.

Install dependencies using the requirements file and then install the package:

    pip3 install -r requirements/install.txt
    pip3 install .

Additional requirements needed for running the unit tests and doing package development are in `reuirements/dev.txt`

## Troubleshooting
If you are a OSX user, you may experience header related issues when installing `pysam`. The current workaround 
is to install pysam version `0.13` manually before installing the requirements:

    pip install pysam==0.13

This is the latest version known to compile without errors.

Although `pysam` is not required for `mavedb-convert` directly, it is installed by some of our dependencies. Until it is removed or made optional by those libraries, `mavedb-convert` will unfortunately not be installable on Windows.
