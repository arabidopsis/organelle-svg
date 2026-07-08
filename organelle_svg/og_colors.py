from __future__ import annotations

import re
from typing import Literal
from typing import TypeAlias
from typing import TypedDict


def get_colors_for_genome(genome: str) -> list[Colour]:
    return Plastome if genome != "mitochondrion" else Chondriome


def mregex(s: str) -> re.Pattern[str]:
    return re.compile(s, re.I)


DOTANY = mregex(".*")


WHITE = 255, 255, 255
RED = 255, 0, 0
GREEN = 0, 255, 0
BLUE = 0, 0, 255

ATP = 151, 190, 13
PSA = 0, 102, 44
PSB = 50, 137, 37
RBCL = 31, 161, 45
PET = 121, 156, 19
TRN = 22, 41, 131
ORF = 87, 185, 168
NDH = 255, 236, 0
CLP = 233, 93, 15
RRN = 226, 0, 26
RPO = 189, 18, 32
RPS = 219, 170, 115
RPL = 158, 119, 66
YCF = 255, 250, 208
ORI = 255, 128, 128
SDH = 52, 211, 77
COB = 200, 250, 40
COX = 255, 180, 255

VIOLET = 171, 37, 157

vvlgrey = 240, 240, 240  # taken from /circos-colors.html

ColourTuple: TypeAlias = tuple[int, int, int] | tuple[int, int, int, float]


class Colour(TypedDict):
    type: re.Pattern[str]  # gene etc.
    pattern: re.Pattern[str]  # e.g. ^psa.* for photosystem I genes
    color: ColourTuple
    fullname: str
    drawflag: bool


PatternKeys: TypeAlias = Literal["type", "pattern"]

Plastome: list[Colour] = [
    {
        "type": mregex("gene"),
        "pattern": mregex("^psa.*"),
        "color": PSA,
        "fullname": "photosystem I",
        "drawflag": True,
    },
    {
        "type": mregex("gene"),
        "pattern": mregex("^psb.*"),
        "color": PSB,
        "fullname": "photosystem II",
        "drawflag": True,
    },
    {
        "type": mregex("gene"),
        "pattern": mregex("^pet.*"),
        "color": PET,
        "fullname": "cytochrome b/f complex",
        "drawflag": True,
    },
    {
        "type": mregex("gene"),
        "pattern": mregex("^atp.*"),
        "color": ATP,
        "fullname": "ATP synthase",
        "drawflag": True,
    },
    {
        "type": mregex("gene"),
        "pattern": mregex("^ndh.*"),
        "color": NDH,
        "fullname": "NADH dehydrogenase",
        "drawflag": True,
    },
    {
        "type": mregex("gene"),
        "pattern": mregex("^rbc[lL].*"),
        "color": RBCL,
        "fullname": "RubisCO large subunit",
        "drawflag": True,
    },
    {
        "type": mregex("gene"),
        "pattern": mregex("^rpo.*"),
        "color": RPO,
        "fullname": "RNA polymerase",
        "drawflag": True,
    },
    {
        "type": mregex("gene"),
        "pattern": mregex("^rps.*"),
        "color": RPS,
        "fullname": "ribosomal proteins (SSU)",
        "drawflag": True,
    },
    {
        "type": mregex("gene"),
        "pattern": mregex("^rpl.*"),
        "color": RPL,
        "fullname": "ribosomal proteins (LSU)",
        "drawflag": True,
    },
    {
        "type": mregex("gene"),
        "pattern": mregex("^(clp|mat).*"),
        "color": CLP,
        "fullname": "clpP, matK",
        "drawflag": True,
    },
    {
        "type": mregex("gene"),
        "pattern": mregex("^ycf.*"),
        "color": YCF,
        "fullname": "hypothetical chloroplast reading frames (ycf)",
        "drawflag": True,
    },
    {
        "type": mregex("CDS|gene"),
        "pattern": mregex("^orf.*"),
        "color": ORF,
        "fullname": "ORFs",
        "drawflag": True,
    },
    {
        "type": mregex("tRNA"),
        "pattern": mregex("trn.*"),
        "color": TRN,
        "fullname": "transfer RNAs",
        "drawflag": True,
    },
    {
        "type": mregex("rRNA"),
        "pattern": DOTANY,
        "color": RRN,
        "fullname": "ribosomal RNAs",
        "drawflag": True,
    },
    {
        "type": mregex("rep_origin"),
        "pattern": mregex("^ori.*"),
        "color": ORI,
        "fullname": "origin of replication",
        "drawflag": True,
    },
    {
        "type": mregex("intron"),
        "pattern": DOTANY,
        # "color": WHITE,
        "color": vvlgrey,
        "fullname": "introns",
        "drawflag": True,
    },
    {
        "type": mregex("other"),
        "pattern": DOTANY,
        "color": VIOLET,
        "fullname": "other genes",
        "drawflag": True,
    },
    {
        "type": mregex("_operon_"),
        "pattern": DOTANY,
        "color": RED,
        "fullname": "polycistronic transcripts",
        "drawflag": False,
    },
]


Chondriome: list[Colour] = [
    {
        "type": mregex("gene"),
        "pattern": mregex("^(nad|nd).*"),
        "color": NDH,
        "fullname": "complex I (NADH dehydrogenase)",
        "drawflag": True,
    },
    {
        "type": mregex("gene"),
        "pattern": mregex("^sdh.*"),
        "color": SDH,
        "fullname": "complex II (succinate dehydrogenase)",
        "drawflag": True,
    },
    {
        "type": mregex("gene"),
        "pattern": mregex("^cob.*"),
        "color": COB,
        "fullname": "complex III (ubichinol cytochrome c reductase)",
        "drawflag": True,
    },
    {
        "type": mregex("gene"),
        "pattern": mregex("^cox.*"),
        "color": COX,
        "fullname": "complex IV (cytochrome c oxidase)",
        "drawflag": True,
    },
    {
        "type": mregex("gene"),
        "pattern": mregex("^atp.*"),
        "color": ATP,
        "fullname": "ATP synthase",
        "drawflag": True,
    },
    {
        "type": mregex("gene"),
        "pattern": mregex("^ccb.*"),
        "color": PSB,
        "fullname": "cytochrome c biogenesis",
        "drawflag": True,
    },
    {
        "type": mregex("gene"),
        "pattern": mregex("^rpo.*"),
        "color": RPO,
        "fullname": "RNA polymerase",
        "drawflag": True,
    },
    {
        "type": mregex("gene"),
        "pattern": mregex("^rps.*"),
        "color": RPS,
        "fullname": "ribosomal proteins (SSU)",
        "drawflag": True,
    },
    {
        "type": mregex("gene"),
        "pattern": mregex("^rpl.*"),
        "color": RPL,
        "fullname": "ribosomal proteins (LSU)",
        "drawflag": True,
    },
    {
        "type": mregex("gene"),
        "pattern": mregex("^(clp|mat).*"),
        "color": CLP,
        "fullname": "maturases",
        "drawflag": True,
    },
    {
        "type": mregex("gene"),
        "pattern": mregex("^orf.*"),
        "color": ORF,
        "fullname": "ORFs",
        "drawflag": True,
    },
    {
        "type": mregex("tRNA"),
        "pattern": DOTANY,
        "color": TRN,
        "fullname": "transfer RNAs",
        "drawflag": True,
    },
    {
        "type": mregex("rRNA"),
        "pattern": DOTANY,
        "color": RRN,
        "fullname": "ribosomal RNAs",
        "drawflag": True,
    },
    {
        "type": mregex("rep_origin"),
        "pattern": mregex("^ori.*"),
        "color": ORI,
        "fullname": "origin of replication",
        "drawflag": True,
    },
    {
        "type": mregex("intron"),
        "pattern": DOTANY,
        # "color": WHITE,
        "color": vvlgrey,
        "fullname": "introns",
        "drawflag": True,
    },
    {
        "type": mregex("_operon_"),
        "pattern": DOTANY,
        "color": RED,
        "fullname": "polycistronic transcripts",
        "drawflag": False,
    },
]
