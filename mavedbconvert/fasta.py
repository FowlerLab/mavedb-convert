import os
import gzip
import bz2


__all__ = ["split_fasta_path", "parse_fasta"]


from . import constants


def split_fasta_path(fname):
    """
    Check that *fname* exists and has a valid FASTQ_ file extension. Valid
    file extensions are ``.fastq`` or ``.fq``, optionally followed by ``.gz``
    or ``.bz2`` if the file is compressed.

    Returns a tuple containing the directory, the file base name with no
    extension, the FASTQ_ file extension used, and the compression format
    (``"gz"``, ``"bz2"``, or ``None``).

    Raises an ``IOError`` if the file doesn't exist. Returns ``None`` if the
    file extension is not recognized.

    Parameters
    ----------
    fname : `str`
        Path to the fastq file.

    Returns
    -------
    `tuple`
        Filename split into ``head``, ``base``, ``extension``,
        ``compression extension``
    """
    if os.path.isfile(fname):
        compression = None
        head, tail = os.path.split(fname)
        base, ext = os.path.splitext(tail)
        if ext.lower() == ".bz2":
            compression = "bz2"
            base, ext = os.path.splitext(base)
        elif ext.lower() == ".gz":
            compression = "gz"
            base, ext = os.path.splitext(base)
        if ext.lower() in (".fa", ".fasta"):
            return head, base, ext, compression
        else:
            raise IOError(
                "Warning: unexpected file extension for '{fname}'".format(fname=fname)
            )
    else:
        raise IOError("file '{fname}' doesn't exist".format(fname=fname))


def parse_fasta(file):
    """Reads the first sequence from a fasta file."""
    _, _, ext, compression = split_fasta_path(file)
    if compression is None and ext in (".fa", ".fasta"):  # raw FASTQ
        open_func = open
    elif compression == "bz2":
        open_func = bz2.open
    elif compression == "gz":
        open_func = gzip.open
    else:
        raise IOError("File extension must be 'gz', 'bz2', 'fa' or 'fasta'.")

    sequence = ""
    sequence_count = 0
    in_sequence = False
    handle = open_func(file, "rt")
    lines = handle.readlines()
    i = 0
    while i < len(lines):
        if i == 0 and not lines[i].strip().startswith(">"):
            raise IOError("'{}' is not a valid FASTA file.".format(file))

        if lines[i].strip().startswith(">"):
            in_sequence = True
            sequence_count += 1
            i += 1
            if sequence_count > 1:
                raise ValueError("Fasta file must contain a single sequence.")

        while i < len(lines) and in_sequence:
            if lines[i].strip().startswith(">"):
                in_sequence = False
            elif not lines[i].strip():
                i += 1
                continue
            elif not constants.dna_re.fullmatch(lines[i].strip().upper()):
                raise ValueError(
                    "Invalid nucleotide characters in line '{}. "
                    "Supported characters are ATCGatcg.".format(lines[i].strip())
                )
            else:
                sequence += lines[i].strip()
                i += 1

    handle.close()
    return sequence.upper()
