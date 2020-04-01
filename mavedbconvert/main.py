"""
Converts an output file from Enrich, Enrich2, or EMPIRIC to a csv file that
can be uploaded to MaveDB.

EMPIRIC users: Please supply an excel/tab separated file with a single row
at the top of the file containing your column names. Required columns are
'Position' and 'Amino Acid'. If you want nucleotide level HGVS events to be
inferred then you must also provide a 'Codon' column. All other columns must be
numeric where NaN values must be encoded using 'NaN', 'Na', 'None', 'N/A',
'undefined', 'Null', or a blank cell.

All outputs are in 1-based coordinates.

Usage:
  mavedbconvert enrich2 <src> [--dst=D] [--wtseq=W] [--offset=O] [--hgvs-column=A] [--input-type=T] [--skip-header=H] [--skip-footer=H] [--non-coding]
  mavedbconvert enrich <src> [--dst=D] [--wtseq=W] [--offset=O]  [--score-column=C] [--input-type=T] [--sheet-name=S] [--skip-header=H] [--skip-footer=H]
  mavedbconvert empiric <src> [--dst=D] [--wtseq=W] [--offset=O] [--zero-based] [--score-column=C] [--input-type=T] [--sheet-name=S] [--skip-header=H] [--skip-footer=H]
  mavedbconvert -h | --help
  mavedbconvert --version
  

Options:
  -h --help         Show this screen.

  --version         Show version.

  <program>         The program the input was generated from.
                    Currently supports enrich, enrich2 and empiric.

  <src>             Path to input file to convert to MaveDB format.

  -d --dst=D        Directory to save the output file to. An attempt will be
                    made to create the directory tree and check write access.
                    If input is a H5 file and a directory is not supplied, a
                    subdirectory with the same name as the input file will be
                    created. In all other cases, saves to the output to the
                    input files's directory if not supplied. [default: None]

  --wtseq=W         A wild-type DNA sequence or fasta file containing the
                    sequence. Required when inputs are from Enrich, Enrich2 or
                    EMPIRIC. Length should be a multiple of 3 for
                    Enrich/EMPIRIC sources or when source is Enrich2
                    and the non-coding flag is absent. [default: None]

  --offset=O        Value to subtract from reported Enrich2 nucleotide positions.
                    For Enrich/EMPIRIC this value should be a multiple of 3 so
                    it can be subtracted from the reported codon positions.
                    [default: 0]
                    
  --input-type=T    Type of input. Can be either 'counts' or 'scores'.
                    Ignored for Enrich2 HD5 files. [default: scores]

  --score-column=C  Column to use as scores. Ignored for Enrich2.
                    [default: None]
  
  --hgvs-column=A   Column containing HGVS variants. Enrich2 TSV only.
                    [default: hgvs]

  --skip-header=H   Integer representing the number of rows to skip at the
                    beginning of an Excel file. [default: 0]

  --skip-footer=F   Integer representing the number of rows to skip at the
                    end of an Excel file. [default: 0]

  --sheet-name=S    The sheet name in an Excel file to parse. [default: None]
  
  --zero-based      Set if the coordinates in an EMPIRIC input file contains
                    one-based coordinates. Ignored for Enrich and Enrich2
                    [default: False]
                    
  --non-coding      Set Enrich2 input file specifies non-coding HGVS syntax.
                    [default: False]
"""
import sys
import docopt
import logging

from . import enrich, enrich2, empiric, constants, LOGGER, parsers


logger = logging.getLogger(LOGGER)


def parse_args(docopt_args=None):
    if docopt_args is None:
        docopt_args = docopt.docopt(__doc__, version="0.4.0-alpha")
    return parsers.parse_docopt(docopt_args)


def main():
    try:
        program, kwargs = parse_args()
        if program == "enrich":
            enrich.Enrich(**kwargs).convert()
        elif program == "enrich2":
            enrich2.Enrich2(**kwargs).convert()
        elif program == "empiric":
            empiric.Empiric(**kwargs).convert()
        else:
            logger.error(
                "Supported programs are {}".format(
                    ", ".join(constants.supported_programs)
                )
            )
            sys.exit()
    except Exception as e:
        logger.exception("A critical error has occurred during conversion.")
        sys.exit(getattr(e, "errno", 0))


if __name__ == "__main__":
    main()
