# fonts to use for SVG diagrams
# Note that it is the *viewer* of the svg (or their machine) which
# will select which font to use.
from __future__ import annotations

import re

SVG_DEFAULT_FONTS = "Roboto, Ubuntu, Helvetica, Arial, sans-serif"
SOURCE = re.compile("^.*(chloroplast|plastid|mitochondrion)(.*)$")
# types of features to include in the SVG
TYPES = {"CDS", "tRNA", "rRNA", "mRNA"}
# names in the feature.qualifiers to try for labeling features in the SVG
TRY_NAMES = ["Name", "gene", "ID", "Parent", "product", "db_xref", "locus_tag"]
