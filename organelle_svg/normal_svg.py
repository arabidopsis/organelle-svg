from __future__ import annotations

from typing import override
from typing import TYPE_CHECKING

from .band_utils import iter_features
from .ogdraw_svg import BaseDraw
from .ogdraw_svg import default_namer
from .ogdraw_svg import DepthDraw
from .ogdraw_svg import get_gene_names
from .ogdraw_svg import get_gene_names_from_features
from .ogdraw_svg import rotate
from .ogdraw_svg import show_band
from .svg_colors import BG_COLOR_SFF
from .svg_colors import STRAND_COLORS
from .svg_utils import annular_path
from .svg_utils import circle
from .svg_utils import group
from .svg_utils import maker
from .svg_utils import path
from .svg_utils import styles_to_classes
from .svg_utils import text_perp

if TYPE_CHECKING:
    from typing import Any
    from typing import Unpack
    from Bio.SeqFeature import SeqFeature
    from Bio.SeqRecord import SeqRecord
    from xml.etree.ElementTree import Element
    from .ogdraw_svg import BaseDrawArgs
    from .svg import TEXT_TYPE
    from .svg_utils import Overlap


class NormalDraw(BaseDraw):
    histogram_colors = DepthDraw.histogram_colors
    strand_colors = STRAND_COLORS
    class_prefix = "normal"
    bg_rec1 = dict(fill=BG_COLOR_SFF, opacity=str(0.2))
    # override
    gc_pos = 0.35
    gc_width = 0.05

    def __init__(
        self,
        rec0: SeqRecord,
        rec1: SeqRecord,
        **kwargs: Unpack[BaseDrawArgs],
    ):
        super().__init__(**kwargs)
        n0 = len(rec0)
        n1 = len(rec1)
        N = max(n0, n1)
        ir_swapped = False
        if self.rotate_image:
            ir = self.get_IR(rec1) if n1 >= n0 else self.get_IR(rec0)
            rot, ir_swapped = rotate(ir, N) if ir else (0, False)
        else:
            rot = 0

        def get_angle(p: float) -> float:
            return (N - p - rot) * 360 / N

        def get_tick_angle(p: float) -> float:
            return (N - p) * 360 / N

        self.genome = self.genome or self.genome_info_class().genome(rec0)
        assert self.genome is not None
        self.rec0 = rec0
        self.rec1 = rec1
        self.get_angle = get_angle
        self.get_tick_angle = get_tick_angle
        self.sw = self.r(1 / self.radius)
        self.rotate_image_angle = rot
        self.ir_swapped = ir_swapped

    def ir_pos(self) -> float:
        return self.ir_pos_with(self.has_seq(self.rec0) or self.has_seq(self.rec1))

    def ticks(self, **attrib: Any) -> None:
        self.base_ticks(
            max(len(self.rec0), len(self.rec1)),
            self.get_tick_angle,
            **attrib,
        )

    def gene_names(self, **attrib: Any) -> None:
        r, get_angle, fs, sw = self.r, self.get_angle, self.fs, self.sw

        try_name = default_namer(self.rec0.name)

        def namer(feat: SeqFeature) -> str:
            return try_name(feat.qualifiers)

        g = group(klass="name")

        txt = r(0.75)
        datad = get_gene_names(
            self.rec0,
            txt,
            fs,
            get_angle,
            subfeatures=self.subfeatures,
            namer=namer,
            strand_key=lambda f: 1,  # all the strands share the same circle space
        )
        if 1 not in datad:  # from strand_key
            return

        data = datad[1]

        inner = circle(0, 0, txt, fill=None, stroke=self.stroke, stroke_width=sw)
        g.append(inner)

        attrib = self.merge_attr(attrib, opacity=1, klass="name", stroke=self.stroke)

        def dotext(
            strand: Any,
            data: list[Overlap],
            fill: str,
            outside: bool = True,
            pos: float = txt,
        ) -> None:
            data = [d for d in data if d.props.strand == strand]
            if not data:
                return
            txt = text_perp(
                data,
                pos + (r(0.0005) if outside else -r(0.0005)),
                fs,
                maker(r(0.2)),
                outside=outside,
                radius=self.radius,
                offset=0.05,  # relative values
                dp=0.02,  # relative values
                fill=fill,
                stroke_width=sw,
                **attrib,
            )
            g.append(txt)

        dotext(1, data, fill=self.strand_colors[1])
        dotext(-1, data, fill=self.strand_colors[-1])
        # other strand posibilities
        dotext(None, data, fill="black")
        dotext(0, data, fill="black")

        # add any missing names from rec1

        done = {namer(feat) for _, feat in iter_features(self.rec0)}
        rec1 = [feat for _, feat in iter_features(self.rec1) if namer(feat) not in done]

        datad = get_gene_names_from_features(
            rec1,
            txt,
            fs,
            get_angle,
            target_length=len(self.rec1.seq or ""),
            subfeatures=self.subfeatures,
            namer=namer,
            strand_key=lambda f: 1,
        )
        if 1 in datad:  # from strand_key above
            _, pos1, dr = self.get_band_pos()
            g.append(
                circle(0, 0, pos1 - dr, fill=None, stroke=self.stroke, stroke_width=sw),
            )
            ov: list[Overlap] = datad[1]
            dotext(1, ov, fill="black", outside=False, pos=pos1 - dr)
            dotext(-1, ov, fill="black", outside=False, pos=pos1 - dr)

        self.g.append(g)

    def get_band_pos(self) -> tuple[float, float, float]:
        return self.r(0.70), self.r(0.58), self.r(0.05)

    def bands(self, **attrib: Any) -> None:
        r, g, get_angle, sw = self.r, self.g, self.get_angle, self.sw

        attrib = self.merge_attr(attrib, **{"stroke": "black", "stroke_width": sw})

        pos0, pos1, dr = self.get_band_pos()
        klass = "band"

        def dobands(rec: SeqRecord, pos: float) -> None:
            g.append(
                circle(
                    0,
                    0,
                    pos,
                    klass=klass,
                    fill=None,
                    stroke=self.stroke,
                    stroke_width=sw,
                ),
            )
            assert self.genome is not None
            show_band(
                g,
                rec,
                pos,
                r,
                genome=self.genome,
                get_angle=get_angle,
                dr=dr,
                **attrib,
            )

        if self.bg_rec1:
            g.append(
                path(
                    d=annular_path(0, 0, pos1 - dr, 2 * dr),
                    klass=klass,
                    **self.bg_rec1,
                ),
            )

        dobands(self.rec0, pos0)
        dobands(self.rec1, pos1)

    def histograms(self, **attrib: Any) -> None:
        get_angle = self.get_angle
        pos, _ = self.get_gc_pos()
        for rec in [self.rec0, self.rec1]:
            if self.has_seq(rec):
                self.do_gc_histogram(rec, get_angle, **attrib)
                self.g.append(
                    circle(
                        0,
                        0,
                        pos,
                        fill=None,
                        stroke=self.stroke,
                        stroke_width=self.sw,
                        klass="gc",
                    ),
                )
                break
        for rec in [self.rec1, self.rec0]:
            if self.is_sff(rec):
                self.do_coverage(rec, get_angle, self.histogram_colors, **attrib)
                break

    def get_rec1_str(self, rec1: SeqRecord) -> str | tuple[str, dict[str, str]]:
        ret = super().get_rec1_str(rec1)
        if self.bg_rec1 and isinstance(ret, str):
            return ret, {"fill": self.bg_rec1["fill"], "font-weight": "bold"}
        return ret

    def doirscan(self, **attrib: Any) -> None:
        # by convention rec1 is an SFF file
        ir = self.base_doirscan(self.rec1, self.get_angle, **attrib)
        if not ir:
            self.base_doirscan(self.rec0, self.get_angle, **attrib)

    def postscript(self, **attrib: Any) -> None:
        # put bling ontop
        self.svg.append(self.g)
        self.add_bling(self.rec0, rec1=self.rec1)

    @override
    def draw(self, **attrib: Any) -> Element:
        self.ticks(**attrib)
        self.gene_names(**attrib)
        self.bands(**attrib)
        self.histograms(**attrib)
        self.doirscan(**attrib)
        self.postscript(**attrib)
        if self.styles_to_classes:
            styles_to_classes(self.svg, self.class_prefix)
        return self.svg


class NCBINormalDraw(NormalDraw):
    def center_text_list(self, rec: SeqRecord) -> tuple[TEXT_TYPE, str]:
        a, genome = super().center_text_list(rec)
        a.append(("(chloe annotation is outer ring)", dict(font_size=30)))
        return a, genome
