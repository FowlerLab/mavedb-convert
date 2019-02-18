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
  mavedb-convert enrich2 <src> [--dst=D] [--wtseq=W] [--offset=O] [--hgvs_column=A] [--non_coding] [--input_type=T] [--skip_header=H] [--skip_footer=H]
  mavedb-convert enrich <src> [--dst=D] [--wtseq=W] [--score_column=C] [--input_type=T] [--sheet_name=S] [--skip_header=H] [--skip_footer=H]
  mavedb-convert empiric <src> [--dst=D] [--wtseq=W] [--one_based] [--score_column=C] [--input_type=T] [--sheet_name=S] [--skip_header=H] [--skip_footer=H]
  mavedb-convert -h | --help
  mavedb-convert --version
  

Options:
  -h --help         Show this screen.

  --version         Show version.

  <program>         The program the input was generated from.
                    Currently supports Enrich, enrich2 and EMPIRIC.

  <src>             Path to input file to convert to MaveDB format.

  -d --dst=D        Directory to save the output file to. An attempt will be
                    made to create the directory tree and check write access.
                    If input is a H5 file and a directory is not supplied, a
                    subdirectory with the same name as the input file will be
                    created. In all other cases, saves to the output to the
                    input files's directory if not supplied. [default: None]

  --wtseq=W         A wildtype DNA sequence or fasta file containing the
                    sequence. Required when inputs are from Enrich, Enrich2 or
                    EMPIRIC. Not used for other input sources. [default: None]

  --one_based       Set if the coordinates in an EMPIRIC input file contains
                    one-based coordinates. Ignored for Enrich and Enrich2
                    [default: False]
                    
  --non_coding      Set Enrich2 input file specifies non-coding HGVS syntax.
                    [default: False]

  --offset=O        Number of bases at the beginning of `wtseq` to ignore
                    before beginning translation. [default: 0]

  --score_column=C  Column to use as scores. Ignored for Enrich2. [default: None]
  
  --hgvs_column=A   Column containing HGVS variants. Enrich2 TSV only. [default: hgvs]

  --input_type=T    Type of input. Can be either 'counts' or 'scores'.
                    Ignored for Enrich2 HD5 files. [default: scores]

  --skip_header=H   Integer representing the number of rows to skip at the
                    beginning of an Excel file. [default: 0]

  --skip_footer=F   Integer representing the number of rows to skip at the
                    end of an Excel file. [default: 0]

  --sheet_name=S    The sheet name in an Excel file to parse. [default: None]
"""
import os
import sys
import docopt
import logging

from . import enrich, enrich2, empiric, constants, LOGGER, fasta


logger = logging.getLogger(LOGGER)


def parse_args(docopt_args=None):
    if docopt_args is None:
        docopt_args = docopt.docopt(__doc__, version='0.2.1-alpha')
    kwargs = {}
    program = None
    for k, v in docopt_args.items():
        if k == '<src>':
            path = os.path.normpath(os.path.expanduser(v))
            try:
                open(path, 'rt').close()
            except FileNotFoundError as e:
                logger.error(
                    "Could not find <src> file '{}'".format(path))
                sys.exit(e.errno)
            except IsADirectoryError as e:
                logger.error(
                    "<src> must be a file not a directory")
                sys.exit(e.errno)
            except PermissionError as e:
                logger.error(
                    "Permission denied for <src> file '{}'".format(path))
                sys.exit(e.errno)
            except IOError as e:
                logger.error(
                    "Unable to open <src> file '{}'".format(path))
                sys.exit(e.errno)
            kwargs[k[1:-1]] = path
        elif k == '--dst':
            if str(v).capitalize() != 'None':
                path = os.path.normpath(os.path.expanduser(v))
                try:
                    if not os.path.isdir(path):
                        os.makedirs(path, exist_ok=True)
                    os.access(path, mode=os.W_OK)
                except FileNotFoundError as e:
                    logger.error(
                        "Could not create directory {}. "
                        "Please ensure it is a valid path.".format(path)
                    )
                    sys.exit(e.errno)
                except PermissionError as e:
                    logger.error(
                        "Permission denied when creating {}.".format(path)
                    )
                    sys.exit(e.errno)
            else:
                path = None
            kwargs[k[2:]] = path
        elif k in constants.supported_programs and v:
            program = k
        elif k in ('--wtseq', '--offset', '--one_based', '--sheet_name',
                   '--score_column', '--hgvs_column', '--input_type',
                   '--skip_header', '--skip_footer', '--non_coding'):
            if isinstance(v, str):
                if v == 'None':
                    kwargs[k[2:]] = None
                else:
                    if k in ('--skip_footer', '--skip_header'):
                        try:
                            kwargs[k[2:]] = int(v)
                        except ValueError:
                            logger.error("{} must be an integer.".format(k))
                            sys.exit()
                    else:
                        kwargs[k[2:]] = v.strip()
            else:
                kwargs[k[2:]] = v
        else:
            continue

    # ---- Validate input types ----- #
    if kwargs['input_type'] not in constants.types:
        logger.error(
            "Supported dataset input types are {}".format(
                ' or '.join(constants.types)))
        sys.exit()
    
    if program not in constants.supported_programs:
        logger.error(
            "Supported programs are {}".format(
                ', '.join(constants.supported_programs)))
        sys.exit()

    # ---- Validate score column ----- #
    if program in ('enrich', 'empiric',):
        defines_column = kwargs['score_column'] and \
                         kwargs['score_column'] is not None
        input_type = kwargs['input_type']
        if input_type == constants.score_type and not defines_column:
            logger.error("A scores column name must be specified.")
            sys.exit()
    
    # ---- Validate offset and wt_seq ----- #
    if program in ('enrich2', ):
        try:
            offset = int(kwargs['offset'])
            # is_coding = kwargs['is_coding']
            # if is_coding and abs(offset) % 3 != 0:
            #     logger.error("Coding offset must be a multiple of three.")
            #     sys.exit()
            kwargs['offset'] = offset
        except (TypeError, ValueError):
            logger.error("Offset must be an integer")
            sys.exit()
            
    if program in ('enrich', 'empiric', 'enrich2'):
        wt_seq = kwargs['wtseq']
        if wt_seq is None or not wt_seq.strip():
            logger.error(
                "A wildtype sequence must be supplied with "
                "Enrich, Enrich2 and Empiric inputs.")
            sys.exit()

        if os.path.isfile(os.path.normpath(os.path.expanduser(wt_seq))):
            wt_seq = fasta.parse_fasta(
                os.path.normpath(os.path.expanduser(wt_seq)))
            kwargs['wtseq'] = wt_seq
        
        if not constants.dna_re.fullmatch(wt_seq):
            logger.error(
                "'{}' is not a valid wildtype DNA sequence.".format(wt_seq))
            sys.exit()
        
        is_mult_of_three = len(wt_seq) % 3 == 0
        if program in ('enrich', 'empiric'):
            if not is_mult_of_three:
                logger.error("Wild-type sequence must be a multiple of three.")
                sys.exit()
        else:
            is_coding = not kwargs['non_coding']
            if is_coding and not is_mult_of_three:
                logger.error("Coding wild-type sequence must be a multiple of three.")
                sys.exit()

    return program, kwargs
    

def main():
    try:
        program, kwargs = parse_args()
        kwargs['wt_sequence'] = kwargs.pop('wtseq')
        kwargs['skip_header_rows'] = kwargs.pop('skip_header')
        kwargs['skip_footer_rows'] = kwargs.pop('skip_footer')
        kwargs['is_coding'] = not kwargs.pop('non_coding')
        if program == 'enrich':
            kwargs.pop('one_based')
            kwargs.pop('offset')
            kwargs.pop('hgvs_column')
            
            enrich.Enrich(**kwargs).convert()
        elif program == 'enrich2':
            kwargs.pop('one_based')
            kwargs.pop('score_column')
            
            enrich2.Enrich2(**kwargs).convert()
        elif program == 'empiric':
            kwargs.pop('offset')
            kwargs.pop('hgvs_column')
            
            empiric.Empiric(**kwargs).convert()
        else:
            logger.error("Supported programs are {}".format(
                ', '.join(constants.supported_programs)))
            sys.exit()
    except Exception:
        logger.exception("An error occured during conversion.")
        sys.exit()
    
        
if __name__ == '__main__':
    main()
