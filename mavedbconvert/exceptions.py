class VariantMismatchError(Exception):
    """
    Throw exception when a variant defines a SNP with a reference base
    that does not match that in its wild-type sequence/codon list.
    """
    pass


class HGVSValidationError(Exception):
    """
    Throw exception when a variant defines a SNP with a reference base
    that does not match that in its wild-type sequence/codon list.
    """
    pass


class WildTypeIndexError(Exception):
    """
    Throw exception when a variant defines a SNP with a position that goes
    out of bounds relative to its wild-type sequence.
    """
    pass


class VariantNotOneBased(Exception):
    """
    Throw exception when a variant is not 1-based.
    """
    pass


class VariantNotZeroBased(Exception):
    """
    Throw exception when a variant is not 0-based.
    """
    pass


class HGVSMatchError(Exception):
    """
    Throw exception when a variant could not be pattern matched.
    """
    pass


class OutOfFrameError(Exception):
    """
    Throw exception when a DNA sequence is out of frame (length is not
    a mulitple of 3)
    """


class HGVSParseError(Exception):
    """
    Throw exception when a a multi variant defines events that span
    multiple codons.
    """