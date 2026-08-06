from __future__ import annotations

import re
from dataclasses import dataclass
from typing import Literal
from typing import TypeAlias


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


@dataclass(kw_only=True)
class GeneColor:
    type: re.Pattern[str]  # gene etc.
    pattern: re.Pattern[str]  # e.g. ^psa.* for photosystem I genes
    color_tuple: ColourTuple
    fullname: str
    drawflag: bool = True

    @property
    def color_str(self) -> str:
        c = self.color_tuple
        if len(c) == 3:
            r, g, b = c
            return f"#{r:2X}{g:2X}{b:2X}".replace(" ", "0")
        return f"rgba({','.join(str(s) for s in c)})"


PatternKeys: TypeAlias = Literal["type", "pattern"]

Plastome: list[GeneColor] = [
    GeneColor(
        type=mregex("gene"),
        pattern=mregex("^psa.*"),
        color_tuple=PSA,
        fullname="photosystem I",
    ),
    GeneColor(
        type=mregex("gene"),
        pattern=mregex("^psb.*"),
        color_tuple=PSB,
        fullname="photosystem II",
    ),
    GeneColor(
        type=mregex("gene"),
        pattern=mregex("^pet.*"),
        color_tuple=PET,
        fullname="cytochrome b/f complex",
    ),
    GeneColor(
        type=mregex("gene"),
        pattern=mregex("^atp.*"),
        color_tuple=ATP,
        fullname="ATP synthase",
    ),
    GeneColor(
        type=mregex("gene"),
        pattern=mregex("^ndh.*"),
        color_tuple=NDH,
        fullname="NADH dehydrogenase",
    ),
    GeneColor(
        type=mregex("gene"),
        pattern=mregex("^rbc[lL].*"),
        color_tuple=RBCL,
        fullname="RubisCO large subunit",
    ),
    GeneColor(
        type=mregex("gene"),
        pattern=mregex("^rpo.*"),
        color_tuple=RPO,
        fullname="RNA polymerase",
    ),
    GeneColor(
        type=mregex("gene"),
        pattern=mregex("^rps.*"),
        color_tuple=RPS,
        fullname="ribosomal proteins (SSU)",
    ),
    GeneColor(
        type=mregex("gene"),
        pattern=mregex("^rpl.*"),
        color_tuple=RPL,
        fullname="ribosomal proteins (LSU)",
    ),
    GeneColor(
        type=mregex("gene"),
        pattern=mregex("^(clp|mat).*"),
        color_tuple=CLP,
        fullname="clpP, matK",
    ),
    GeneColor(
        type=mregex("gene"),
        pattern=mregex("^ycf.*"),
        color_tuple=YCF,
        fullname="hypothetical chloroplast reading frames (ycf)",
    ),
    GeneColor(
        type=mregex("CDS|gene"),
        pattern=mregex("^orf.*"),
        color_tuple=ORF,
        fullname="ORFs",
    ),
    GeneColor(
        type=mregex("tRNA"),
        pattern=mregex("trn.*"),
        color_tuple=TRN,
        fullname="transfer RNAs",
    ),
    GeneColor(
        type=mregex("rRNA"),
        pattern=DOTANY,
        color_tuple=RRN,
        fullname="ribosomal RNAs",
    ),
    GeneColor(
        type=mregex("rep_origin"),
        pattern=mregex("^ori.*"),
        color_tuple=ORI,
        fullname="origin of replication",
    ),
    GeneColor(
        type=mregex("intron"),
        pattern=DOTANY,
        # "color": WHITE,
        color_tuple=vvlgrey,
        fullname="introns",
    ),
    GeneColor(
        type=mregex("other"),
        pattern=DOTANY,
        color_tuple=VIOLET,
        fullname="other genes",
    ),
    GeneColor(
        type=mregex("_operon_"),
        pattern=DOTANY,
        color_tuple=RED,
        fullname="polycistronic transcripts",
        drawflag=False,
    ),
]


Chondriome: list[GeneColor] = [
    GeneColor(
        type=mregex("gene"),
        pattern=mregex("^(nad|nd).*"),
        color_tuple=NDH,
        fullname="complex I (NADH dehydrogenase)",
    ),
    GeneColor(
        type=mregex("gene"),
        pattern=mregex("^sdh.*"),
        color_tuple=SDH,
        fullname="complex II (succinate dehydrogenase)",
    ),
    GeneColor(
        type=mregex("gene"),
        pattern=mregex("^cob.*"),
        color_tuple=COB,
        fullname="complex III (ubichinol cytochrome c reductase)",
    ),
    GeneColor(
        type=mregex("gene"),
        pattern=mregex("^cox.*"),
        color_tuple=COX,
        fullname="complex IV (cytochrome c oxidase)",
    ),
    GeneColor(
        type=mregex("gene"),
        pattern=mregex("^atp.*"),
        color_tuple=ATP,
        fullname="ATP synthase",
    ),
    GeneColor(
        type=mregex("gene"),
        pattern=mregex("^ccb.*"),
        color_tuple=PSB,
        fullname="cytochrome c biogenesis",
    ),
    GeneColor(
        type=mregex("gene"),
        pattern=mregex("^rpo.*"),
        color_tuple=RPO,
        fullname="RNA polymerase",
    ),
    GeneColor(
        type=mregex("gene"),
        pattern=mregex("^rps.*"),
        color_tuple=RPS,
        fullname="ribosomal proteins (SSU)",
    ),
    GeneColor(
        type=mregex("gene"),
        pattern=mregex("^rpl.*"),
        color_tuple=RPL,
        fullname="ribosomal proteins (LSU)",
    ),
    GeneColor(
        type=mregex("gene"),
        pattern=mregex("^(clp|mat).*"),
        color_tuple=CLP,
        fullname="maturases",
    ),
    GeneColor(
        type=mregex("gene"),
        pattern=mregex("^orf.*"),
        color_tuple=ORF,
        fullname="ORFs",
    ),
    GeneColor(
        type=mregex("tRNA"),
        pattern=DOTANY,
        color_tuple=TRN,
        fullname="transfer RNAs",
    ),
    GeneColor(
        type=mregex("rRNA"),
        pattern=DOTANY,
        color_tuple=RRN,
        fullname="ribosomal RNAs",
    ),
    GeneColor(
        type=mregex("rep_origin"),
        pattern=mregex("^ori.*"),
        color_tuple=ORI,
        fullname="origin of replication",
    ),
    GeneColor(
        type=mregex("intron"),
        pattern=DOTANY,
        # "color": WHITE,
        color_tuple=vvlgrey,
        fullname="introns",
    ),
    GeneColor(
        type=mregex("_operon_"),
        pattern=DOTANY,
        color_tuple=RED,
        fullname="polycistronic transcripts",
        drawflag=False,
    ),
]

COLOR_DICT = {"mitochondrion": Chondriome, "plastid": Plastome, "chloroplast": Plastome}


def get_colors_for_genome(genome_type: str) -> list[GeneColor]:
    return COLOR_DICT.get(genome_type, Plastome)
