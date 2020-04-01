import re

from hgvsp import dna

from hgvs.dataproviders import uta
from hgvs.parser import Parser


supported_programs = ("enrich", "enrich2", "empiric")
extra_na = (
    "None",
    "none",
    "NONE",
    "undefined",
    "Undefined",
    "UNDEFINED",
    "na",
    "Na",
    "N/a",
    "Null",
    "",
    " ",
)
null_value_re = re.compile(r"\s+|nan|na|none|undefined|n/a|null")
surrounding_brackets_re = re.compile(r"\((.*)\)")
dna_re = re.compile(r"[ATCGatcg]+", flags=re.IGNORECASE)

# HGVS
dummy_ref = "NM_000000000.0"
hgvs_parser = Parser()
hdp = None


def get_hdp(*args, **kwargs):
    """Keep a single HDP network connection instance."""
    global hdp
    if hdp is None:
        if "pooling" not in kwargs:
            kwargs["pooling"] = True
        hdp = uta.connect(*args, **kwargs)
    return hdp


# HGVSP constants
hgvsp_nt_pos = "position"
hgvsp_pro_pos = "position"
hgvsp_nt_ref = "ref"
hgvsp_nt_alt = "alt"
hgvsp_pro_ref = "pre"
hgvsp_pro_alt = "post"
hgvsp_silent = "silent"


# Enrich2 constants
enrich2_synonymous = "_sy"
enrich2_wildtype = "_wt"
special_variants = ("_wt", "_sy")
synonymous_table = "synonymous"
variants_table = "variants"


# MaveDB constants
nt_variant_col = "hgvs_nt"
pro_variant_col = "hgvs_pro"
variant_columns = (nt_variant_col, pro_variant_col)
score_type = "scores"
count_type = "counts"
types = (score_type, count_type)
mavedb_score_column = "score"


# Conversions between single- and three-letter amino acid codes
AA_CODES = {
    "Ala": "A",
    "A": "Ala",
    "Arg": "R",
    "R": "Arg",
    "Asn": "N",
    "N": "Asn",
    "Asp": "D",
    "D": "Asp",
    "Cys": "C",
    "C": "Cys",
    "Glu": "E",
    "E": "Glu",
    "Gln": "Q",
    "Q": "Gln",
    "Gly": "G",
    "G": "Gly",
    "His": "H",
    "H": "His",
    "Ile": "I",
    "I": "Ile",
    "Leu": "L",
    "L": "Leu",
    "Lys": "K",
    "K": "Lys",
    "Met": "M",
    "M": "Met",
    "Phe": "F",
    "F": "Phe",
    "Pro": "P",
    "P": "Pro",
    "Ser": "S",
    "S": "Ser",
    "Thr": "T",
    "T": "Thr",
    "Trp": "W",
    "W": "Trp",
    "Tyr": "Y",
    "Y": "Tyr",
    "Val": "V",
    "V": "Val",
    "Ter": "*",
    "*": "Ter",
    "???": "?",
    "?": "???",
    "X": "Xaa",
    "Xaa": "X",
}
