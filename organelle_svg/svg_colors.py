from __future__ import annotations

from typing import TYPE_CHECKING


from .og_colors import get_colors_for_genome


if TYPE_CHECKING:
    from typing import Callable
    from .og_colors import PatternKeys

HISTOGRAM_COLORS = {"depth": "rgba(250,0,0,.6)", "coverage": "rgba(0,0,255,.6)"}
STRAND_COLORS = {1: "rgb(239,59,44)", -1: "rgb(66,146,198)"}
# background color for SFF in normal plot
BG_COLOR_SFF = "#aaaaff"


def colorer(
    genome_type: str, *, by: PatternKeys = "pattern", default: str = "white"
) -> Callable[[str], str]:
    colours = get_colors_for_genome(genome_type)
    C = [(getattr(d, by), d.color_str) for d in colours if d.drawflag]

    def color(gene: str) -> str:

        for m, color in C:
            if m.match(gene):
                return color

        return default

    return color
