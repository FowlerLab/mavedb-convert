import os
import logging
from fqfa.fasta.fasta import parse_fasta_records
from fqfa.validator.validator import dna_bases_validator

from . import LOGGER, constants, exceptions


logger = logging.getLogger(LOGGER)


def parse_boolean(value):
    return str(value) == "True"


def parse_numeric(value, name=None, dtype=int):
    try:
        return dtype(value)
    except (ValueError, TypeError):
        raise ValueError("Arg '{}' must be type {}".format(name, dtype))


def parse_string(value):
    if not value:
        return None
    if not str(value).strip():
        return None
    if str(value).strip().lower() == "none":
        return None
    return str(value).strip()


def parse_src(src):
    src = parse_string(src)
    if not src:
        raise ValueError("<src> argument is required.")
    path = os.path.normpath(os.path.expanduser(src))
    try:
        open(path, "rt").close()
    except FileNotFoundError as e:
        logger.error("Could not find <src> file '{}'".format(path))
        raise e
    except IsADirectoryError as e:
        logger.error("<src> must be a file not a directory")
        raise e
    except PermissionError as e:
        logger.error("Permission denied for <src> file '{}'".format(path))
        raise e
    except IOError as e:
        logger.error("Unable to open <src> file '{}'".format(path))
        raise e
    return path


def parse_dst(dst):
    dst = parse_string(dst)
    if dst:
        path = os.path.normpath(os.path.expanduser(dst))
        try:
            if not os.path.isdir(path):
                os.makedirs(path, exist_ok=True)
            os.access(path, mode=os.W_OK)
        except PermissionError as e:
            logger.error("Permission denied when creating {}.".format(path))
            raise e
        return path
    return None


def parse_program(program):
    if isinstance(program, dict):
        if program.get("enrich", False):
            program = parse_program("enrich")
        elif program.get("enrich2", False):
            program = parse_program("enrich2")
        elif program.get("empiric", False):
            program = parse_program("empiric")
        else:
            raise ValueError("<program> is required.")
    if program not in constants.supported_programs:
        raise ValueError("{} is not a recognised format.".format(program))
    return program


def parse_wt_sequence(wtseq, coding=True):
    if os.path.isfile(os.path.normpath(os.path.expanduser(wtseq))):
        with open(os.path.normpath(os.path.expanduser(wtseq))) as fh:
            _, wtseq = next(parse_fasta_records(fh))

    if not dna_bases_validator(wtseq.upper()):
        raise exceptions.InvalidWildTypeSequence(
            "Wild-type sequence contains invalid characters."
        )

    if coding and len(wtseq) % 3 != 0:
        raise exceptions.SequenceFrameError(
            f"Enrich2 wild-type sequence for a coding dataset "
            "must be a multiple of three. Found length {len(wtseq)}."
        )

    return wtseq.upper()


def parse_input_type(value):
    value = parse_string(value)
    if value not in (constants.score_type, constants.count_type):
        raise ValueError(
            "'{}' is not a recognised input type. "
            "Choices are '{}' or '{}'.".format(
                value, constants.score_type, constants.count_type
            )
        )
    return value


def parse_score_column(value, input_type, program):
    value = parse_string(value)
    defines_column = value is not None
    if program in ("enrich", "empiric"):
        if input_type == constants.score_type and not defines_column:
            raise ValueError(
                "A scores column name must be "
                "specified if input type is '{}'.".format(constants.score_type)
            )
    return value


def parse_offset(offset, coding=True):
    offset = parse_numeric(offset, name="offset", dtype=int)
    mult_of_three = abs(offset) % 3 == 0
    if coding and not mult_of_three:
        raise ValueError("Offset for a coding dataset must be a multiple of three.")

    return offset


def parse_docopt(docopt_args):
    parsed_kwargs = {}

    program = parse_program(docopt_args)
    parsed_kwargs["src"] = parse_src(docopt_args.get("<src>", None))
    parsed_kwargs["dst"] = parse_dst(docopt_args.get("--dst", None))

    # Parse booleans
    parsed_kwargs["is_coding"] = not parse_boolean(
        docopt_args.get("--non-coding", False)
    )
    parsed_kwargs["one_based"] = not parse_boolean(
        docopt_args.get("--zero-based", False)
    )

    # Parse WT and Offset fields
    parsed_kwargs["wt_sequence"] = parse_wt_sequence(
        docopt_args.get("--wtseq", None),
        coding=parsed_kwargs["is_coding"],
    )
    parsed_kwargs["offset"] = parse_offset(
        docopt_args.get("--offset", 0),
        coding=parsed_kwargs["is_coding"],
    )

    # Parse Input related fields
    parsed_kwargs["input_type"] = parse_input_type(
        docopt_args.get("--input-type", None)
    )
    parsed_kwargs["score_column"] = parse_score_column(
        docopt_args.get("--score-column", None),
        program=program,
        input_type=parsed_kwargs["input_type"],
    )
    parsed_kwargs["hgvs_column"] = parse_string(docopt_args.get("--hgvs-column", None))

    # Parse Excel related fields
    parsed_kwargs["sheet_name"] = parse_string(docopt_args.get("--sheet_name", None))
    parsed_kwargs["skip_header_rows"] = parse_numeric(
        docopt_args.get("--skip-header", 0), name="skip_header", dtype=int
    )
    parsed_kwargs["skip_footer_rows"] = parse_numeric(
        docopt_args.get("--skip-footer", 0), name="skip_footer", dtype=int
    )
    return program, parsed_kwargs
