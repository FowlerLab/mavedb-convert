# mavedb-convert
A command line tool for converting alternate file formats into a MaveDB compliant format.

# Attention!
If you are a OSX user, you may experience header related issues when installing `pysam`. The current workaround 
is to install pysam version `0.13` manually:

`pip install pysam==0.13`

This is the latest version known to compile without errors.

Although `pysam` is not required for `mavedb-convert` directly, it is installed by some of our dependencies. Until it is removed or made optional by those libraries, `mavedb-convert` will unfortunately not be installable on Windows.
