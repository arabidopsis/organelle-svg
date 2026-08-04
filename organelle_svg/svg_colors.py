from __future__ import annotations

from typing import TYPE_CHECKING


from .og_colors import get_colors_for_genome


if TYPE_CHECKING:
    from typing import Callable
    from .og_colors import PatternKeys, ColourTuple

HISTOGRAM_COLORS = {"depth": "rgba(250,0,0,.6)", "coverage": "rgba(0,0,255,.6)"}
STRAND_COLORS = {1: "rgb(239,59,44)", -1: "rgb(66,146,198)"}
# background color for SFF in normal plot
BG_COLOR_SFF = "#aaaaff"


def colorer(genome: str, *, by: PatternKeys = "pattern") -> Callable[[str], str]:
    colours = get_colors_for_genome(genome)
    C = [(d[by], d["color"]) for d in colours if d["drawflag"]]

    def color(gene: str) -> str:
        c: ColourTuple = 0, 0, 0

        for m, color in C:
            if m.match(gene):
                c = color
                break
        a = "" if len(c) == 3 else "a"
        return f"rgb{a}({','.join(str(s) for s in c)})"

    return color
