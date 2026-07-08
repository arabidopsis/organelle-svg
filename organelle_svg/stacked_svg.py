from __future__ import annotations

from typing import override
from typing import TYPE_CHECKING

from intervaltree import IntervalTree

from .ogdraw_svg import BaseDraw
from .ogdraw_svg import create_band
from .svg import add_center_text
from .svg_colors import colorer
from .svg_utils import circle
from .svg_utils import group
from .svg_utils import styles_to_classes
from .svg_utils import ticks

if TYPE_CHECKING:
    from typing import Any
    from typing import Unpack
    from typing import Iterator
    from typing import Sequence
    from xml.etree.ElementTree import Element
    from Bio.SeqRecord import SeqRecord
    from Bio.SeqFeature import SeqFeature
    from .ogdraw_svg import BaseDrawArgs


class StackedDraw(BaseDraw):
    class_prefix = "stacks"
    band_loc = 0.6
    tick_loc = 0.8
    dr = 0.02

    def __init__(
        self,
        name: str,
        recs: Sequence[SeqRecord],
        type: str = "chloroplast",
        **kwargs: Unpack[BaseDrawArgs],
    ):
        super().__init__(**kwargs)
        N = max(len(r) for r in recs)

        def get_angle(pos: float) -> float:
            return (N - pos) * 360 / N

        self.get_angle = get_angle
        self.recs = recs
        self.N = N
        self.stacks: dict[int | None, set[int]] = {}
        self.name = name
        self.colorf = colorer(type)

    def ir_pos(self) -> float:
        return self.ir_pos_with(any(self.has_seq(r) for r in self.recs))

    def ticks(self) -> None:
        # loc = self.r(self.band_loc)
        tloc = self.r(self.tick_loc)
        dr = self.r(self.dr)
        g = group(klass="ticks")

        ticks(
            self.N,
            tloc,
            g,
            self.r,
            get_angle=self.get_angle,
            stroke="grey",
            grid=(tloc - dr / 2, self.r(0.3)),
            grid_opacity=0.4,
            font_size=self.r(0.03),
            opacity=0.5,
            # text_fill="#ffffff",
        )
        g.append(circle(0, 0, r=tloc, fill=None, stroke="grey", stroke_width=1))
        self.g.append(g)

    def layers(self) -> Iterator[tuple[SeqFeature, int]]:
        itree: dict[int | None, IntervalTree] = {
            1: IntervalTree(),
            -1: IntervalTree(),
        }
        itree[0] = itree[-1]
        itree[None] = itree[1]

        stacks: dict[int | None, set[int]] = {1: set(), -1: set()}
        stacks[0] = stacks[-1]
        stacks[None] = stacks[1]
        for rec in self.recs:
            for feat in rec.features:
                f = feat.location
                if f is None:
                    continue
                s, e = f.start, f.end
                # assert isinstance(s, ExactPosition) and isinstance(e, ExactPosition)
                it = itree[f.strand]
                nn = stacks[f.strand]
                oo = it[s:e]
                if oo:
                    # if (f.start, f.end) in {(o.begin, o.end) for o in oo}:
                    #     continue
                    s = {o.data for o in oo}
                    n = max(s) + 1
                    for m in range(n):
                        if m not in s:
                            n = m
                            break
                else:
                    n = 0
                nn.add(n)
                it[f.start : f.end] = n
                n = n if f.strand == 1 else -n
                yield feat, n

        self.stacks = stacks

    def bands(self) -> None:
        loc = self.r(self.band_loc)
        get_angle, r, colorf = self.get_angle, self.r, self.colorf
        dr = r(self.dr)

        g = group(klass="band")

        for feat, n in self.layers():
            if feat.location is None:
                continue
            gene = self.get_gene(feat)
            g.append(
                create_band(
                    loc + n * dr,
                    feat.location,
                    color=colorf(gene),
                    dr=dr,
                    get_angle=get_angle,
                    stroke="grey",
                    stroke_width=1,
                ),
            )

        g.append(circle(0, 0, r=loc, fill=None, stroke="grey", stroke_width=3))
        self.g.append(g)

    def get_gene(self, feat: SeqFeature) -> str:
        if "gene" in feat.qualifiers:
            return feat.qualifiers["gene"][0]
        return feat.type

    def postscript(self) -> None:
        ngenomes = len({r.id for r in self.recs})
        # ["Chloë"]
        text: list[str | tuple[str, Any]] = [
            self.name,
            f"genomes={ngenomes}",
        ]
        w = 2 * self.radius
        self.svg.append(self.g)
        add_center_text(
            self.svg,
            text,
            font_size=self.r(0.08),
            y_offset=0,  # self.r(0.02),
            # color="#ffffff",
            inline=True,
            width=w,
            height=w,
            opacity=0.6,
        )
        if self.styles_to_classes:
            styles_to_classes(self.svg, self.class_prefix)

    @override
    def draw(self, **attrib: Any) -> Element:
        self.ticks()
        self.bands()
        self.postscript()
        return self.svg
