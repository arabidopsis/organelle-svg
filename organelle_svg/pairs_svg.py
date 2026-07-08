from __future__ import annotations

from typing import override
from typing import TYPE_CHECKING

from .bands import get_links
from .ogdraw_svg import BaseDraw
from .ogdraw_svg import get_gene_names
from .ogdraw_svg import show_band
from .svg_colors import STRAND_COLORS
from .svg_utils import arc
from .svg_utils import circle
from .svg_utils import group
from .svg_utils import maker
from .svg_utils import middle
from .svg_utils import Overlap
from .svg_utils import path
from .svg_utils import ribbon_path
from .svg_utils import styles_to_classes
from .svg_utils import text_arc
from .svg_utils import text_horz
from .svg_utils import text_len
from .svg_utils import text_perp
from .svg_utils import ticks

if TYPE_CHECKING:
    from typing import Any
    from typing import Callable
    from typing import Unpack
    from xml.etree.ElementTree import Element
    from Bio.SeqRecord import SeqRecord
    from .ogdraw_svg import BaseDrawArgs


class Chromosome:
    strand_colors = STRAND_COLORS
    radial_title: bool = True

    def __init__(
        self,
        rec: SeqRecord,
        r: Callable[[float], float],
        start: int,
        end: int,
        radius: float,
        title: str | None = None,
        units: int = 1000,
        genome: str | None = None,
    ):
        self.rec = rec
        self.start = start
        self.end = end
        self.r = r
        self.units = units
        assert genome is not None
        self.genome = genome
        self.title_ = title or self.rec.id or ""
        self.radius = radius

        def get_angle(base: float) -> float:
            return base * (end - start) / len(rec) + start

        self.get_angle = get_angle

    def ticks(
        self,
        g: Element,
        loc: float,
        tick_info: dict[str, Any] | None = None,
        **attrib: Any,
    ) -> None:
        attrib = {**dict(stroke="grey", opacity=1, fill="grey"), **attrib}
        r = self.r
        ticks(
            len(self.rec),
            loc,
            g,
            get_angle=self.get_angle,
            units=self.units,
            grid=(r(0.25), r(0.75)),
            font_size=r(0.03),
            r=r,
            tick_info=tick_info,
            grid_opacity=0.2,
            **attrib,
        )

    def names(self, g: Element, pos: float, fs: float, **attrib: Any) -> None:
        r = self.r
        sw = r(1 / self.radius)
        attrib = {**dict(stroke_width=sw, opacity=1, stroke="grey"), **attrib}

        def dotext(data: list[Overlap], strand: Any, fill: str) -> None:
            attr = {"fill": fill, "klass": "name", **attrib}
            txt = text_perp(
                [d for d in data if d.props.strand == strand],
                r0=pos,
                fs=fs,
                r=maker(r(0.2)),
                outside=True,
                radius=self.radius,
                offset=0.1,
                dp=0.05,
                **attr,
            )
            g.append(txt)

        # use fs * r to get more squeezed text
        efs = fs * 0.7
        res = get_gene_names(
            self.rec,
            r0=pos,
            fs=efs,
            get_angle=self.get_angle,
            strand_key=lambda f: 1,
        )
        if 1 not in res:
            return
        data = res[1]
        dotext(data, strand=1, fill=self.strand_colors[1])
        dotext(data, strand=-1, fill=self.strand_colors[-1])

    def bands(
        self,
        g: Element,
        pos: float,
        bw: float = 0.5,
        **attrib: Any,
    ) -> None:
        # rec = self.rec
        r = self.r
        g.append(
            arc(
                pos,
                self.start,
                self.end,
                fill=None,
                stroke="grey",
                stroke_width=r(1 / self.radius),
            ),
        )
        show_band(
            g,
            self.rec,
            pos,
            r,
            self.genome,
            get_angle=self.get_angle,
            dr=r(bw),
            **attrib,
        )

    def title(
        self,
        g: Element,
        fs: float,
        offset: float = 0,
        extra: str = "",
        **attrib: Any,
    ) -> None:
        r = self.r
        r0 = r(0.9)
        fs *= 2
        angle = middle(self.start, self.end)
        if not self.radial_title:
            txt = text_horz(
                [Overlap(angle, self.title_, delta_r=offset)],
                r0=r0,
                fs=fs,
                r=r,
                klass="title",
                **attrib,
            )
        else:
            angle -= text_len(self.title_, fs) / fs
            txt = text_arc(self.title_, angle, r0 + offset / 1.5, font_size=fs)
        g.append(txt)


def do_links(
    chr1: Chromosome,
    chr2: Chromosome,
    r0: float,
    g: Element,
    r: Callable[[float], float],
    **attrib: Any,
) -> None:
    rec1, rec2 = chr1.rec, chr2.rec

    # attrib.pop("fill", None)
    g2 = group(klass="links")
    for m, l1, l2, _ in get_links(rec1, rec2):
        if m == "matched":
            attr = {
                "data-info": f"{l1.gene}[{l1.start + 1}..{l1.end}][{l2.start + 1}..{l2.end}]",
                **attrib,
            }

            p = path(
                d=ribbon_path(
                    r0,
                    chr1.get_angle(l1.start),
                    chr1.get_angle(l1.end),
                    r0,
                    chr2.get_angle(l2.end),
                    chr2.get_angle(l2.start),
                    bezier_radius=r(0.20),
                ),
                **attr,
            )
            g2.append(p)
        else:
            if l1.name == rec1.id:
                get_angle = chr1.get_angle
            else:
                get_angle = chr2.get_angle
            dr = r(0.025)
            attr = {
                "data-info": f"{l1.gene}[{l1.start + 1}..{l1.end}]",
                "fill": "black",
                **attrib,
            }
            a = arc(
                r(0.50),
                get_angle(l1.start),
                get_angle(l1.end),
                dr=dr if l1.strand > 0 else -dr,
                **attr,
            )
            g.append(a)

        g.append(g2)


class PairsDraw(BaseDraw):
    class_prefix = "pairs"
    ribbon_attrs = dict(stroke="grey", fill="pink", opacity=0.1)

    def __init__(
        self,
        rec0: SeqRecord,
        rec1: SeqRecord,
        title: str | None = None,
        **kwargs: Unpack[BaseDrawArgs],
    ):
        super().__init__(**kwargs)

        self.genome = self.genome or self.genome_info_class().genome(rec0)

        chr1 = Chromosome(rec0, self.r, -80, 80, radius=self.radius, genome=self.genome)
        chr2 = Chromosome(
            rec1,
            self.r,
            260,
            100,
            radius=self.radius,
            title=title,
            genome=self.genome,
        )

        self.chr1 = chr1
        self.chr2 = chr2

    def get_tick_pos(self) -> float:
        return self.r(0.8)

    def ir_pos(self) -> float:
        return self.ir_pos_with(
            self.has_seq(self.chr1.rec) or self.has_seq(self.chr2.rec),
        )

    def ticks(self, **attrib: Any) -> None:
        g, r = self.g, self.r
        tick_pos = self.get_tick_pos()
        attrib = self.merge_attr(
            attrib,
            stroke="black",
            stroke_width=r(2 / self.radius),
            opacity=1,
        )
        attrib.pop("fill", None)
        g.append(circle(0, 0, tick_pos, fill=None, **attrib))

    def chromosomes(self, **attrib: Any) -> None:
        nattrib = self.merge_attr(attrib)
        battrib = self.merge_attr(
            attrib,
            stroke="black",
            stroke_width=self.r(1 / self.radius),
        )
        tattrib = self.merge_attr(attrib, tick_info={10: dict(height=12 / 100)})
        tick_pos = self.get_tick_pos()
        g, r, fs = self.g, self.r, self.fs
        extra = ""
        for c in [self.chr1, self.chr2]:
            # c.ticks(g, tick_pos)
            c.ticks(g, tick_pos, **tattrib)
            c.names(g, tick_pos + r(0.01), fs, **nattrib)
            c.bands(g, r(0.75), bw=0.05, **battrib)
            c.title(
                g,
                fs,
                offset=r(0.08),
                extra=extra if c == self.chr2 else "",
                **nattrib,
            )

    def links(self, **attrib: Any) -> None:
        chr1, chr2 = self.chr1, self.chr2
        g, r = self.g, self.r
        attr = self.merge_attr(attrib, **self.ribbon_attrs)
        do_links(chr1, chr2, r(0.75), g, r, **attr)
        attrib.pop("fill", None)
        # g.append(circle(0, 0, r(0.5), fill=None, stroke="black"))
        g.append(arc(r(0.5), -80, 80, fill=None, **attrib))
        g.append(arc(r(0.5), 260, 100, fill=None, **attrib))

    def postscript(self, **attrib: Any) -> None:
        # put bling on top
        self.svg.append(self.g)
        self.add_bling(
            self.chr1.rec,
            rec1=self.chr2.rec,
            **attrib,
        )

    @override
    def draw(self, **attrib: Any) -> Element:
        self.ticks(**attrib)
        self.chromosomes(**attrib)
        self.links(**attrib)
        self.postscript(**attrib)
        if self.styles_to_classes:
            styles_to_classes(self.svg, self.class_prefix)
        return self.svg
