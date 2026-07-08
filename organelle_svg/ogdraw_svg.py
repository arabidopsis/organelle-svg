from __future__ import annotations

import abc
from collections import defaultdict
from math import pi
from typing import Any
from typing import override
from typing import TYPE_CHECKING
from typing import TypedDict

from .band_utils import has_seq as has_seq_data
from .band_utils import iter_features
from .band_utils import Overlapping
from .histograms import coverage_histogram
from .histograms import depth_histogram
from .histograms import gc_histogram
from .irscan import IRScan
from .irscan import IRScanResult
from .svg import add_center_text
from .svg import add_legend
from .svg import savesvg
from .svg import tobytes
from .svg import tostring
from .svg_colors import colorer
from .svg_colors import GC_HIST_BACKGROUND
from .svg_colors import GenomeInfo
from .svg_colors import HISTOGRAM_COLORS
from .svg_utils import arc
from .svg_utils import circle
from .svg_utils import fix_text_overlap
from .svg_utils import group
from .svg_utils import klass
from .svg_utils import maker
from .svg_utils import middle
from .svg_utils import Overlap
from .svg_utils import radial_line
from .svg_utils import rect
from .svg_utils import styles_to_classes
from .svg_utils import svge
from .svg_utils import text_horz
from .svg_utils import text_len
from .svg_utils import text_perp
from .svg_utils import ticks

if TYPE_CHECKING:
    from typing import Callable
    from typing import Unpack
    from typing import Sequence
    from typing import Iterator
    from xml.etree.ElementTree import Element
    from Bio.SeqFeature import SeqFeature
    from Bio.SeqRecord import SeqRecord
    from Bio.SeqFeature import CompoundLocation
    from Bio.SeqFeature import SimpleLocation
    from .svg import TEXT_TYPE
    from .band_utils import Overlap as IOverlap


def create_band(
    r0: float,
    part: SimpleLocation | CompoundLocation,
    color: str,
    dr: float,
    get_angle: Callable[[float | int], float],
    gene: str | None = None,
    **attrib: Any,
) -> Element:
    if gene:
        c = ""
        if hasattr(part, "comment"):
            c = getattr(part, "comment") or ""
            if c:
                c = " " + str(c)
            else:
                c = ""
        fmt = f"{gene}:[{part.start + 1}..{part.end}]{c}"  # type: ignore
        attr = {"data-info": fmt, **attrib}
    else:
        attr = attrib
    start, end = part.start, part.end
    angle1 = get_angle(start)  # type: ignore
    angle2 = get_angle(end)  # type: ignore
    w = dr if part.strand is None or part.strand > 0 else -dr
    a = arc(r0, angle2, angle1, dr=w, fill=color, **attr)
    return a


def show_band(
    g2: Element,
    rec: SeqRecord,
    r0: float,
    r: Callable[[float], float],
    genome: str,
    get_angle: Callable[[int], float],
    dr: float = 5,
    show_span: bool = False,
    **attrib: Any,
) -> None:
    g = group(klass=klass(attrib, "band"))

    # N = len(rec)

    getcolor = colorer(genome)
    getcolort = colorer(genome, by="type")

    icolor = getcolort("intron")

    def isintron(loc: SimpleLocation) -> bool:
        if hasattr(loc, "intron"):
            return getattr(loc, "intron")
        return False

    namer = default_namer(rec.name)

    # use Proteome colors
    attrib.pop("fill", None)

    def create_glyph(gene: str, part: SimpleLocation, color: str) -> Element:
        c = ""
        if hasattr(part, "comment"):
            c = getattr(part, "comment") or ""
            if c:
                c = " " + str(c)
            else:
                c = ""
        fmt = f"{gene}:[{part.start + 1}..{part.end}]{c}"  # type: ignore
        attr = {"data-info": fmt, **attrib}
        start, end = part.start, part.end
        angle1 = get_angle(start)  # type: ignore
        angle2 = get_angle(end)  # type: ignore
        w = dr if part.strand and part.strand > 0 else -dr
        a = arc(r0, angle2, angle1, dr=w, fill=color, **attr)
        return a

    # do any introns first below rest

    for _, feat in iter_features(rec):
        if feat.location is None:
            continue
        if feat.type == "intron" and feat.location is not None:
            g.append(create_glyph("intron", feat.location, icolor))  # type: ignore
        else:
            for part in feat.location.parts:
                if isintron(part):
                    g.append(create_glyph("intron", part, icolor))

    for _, feat in iter_features(rec):
        if feat.location is None:
            continue
        if feat.type == "intron":
            continue
        q = feat.qualifiers
        gene = namer(q)

        color = getcolor(gene)

        for part in feat.location.parts:
            if isintron(part):
                # already drawn
                continue
            a = create_glyph(gene, part, color)
            g.append(a)

    if show_span:
        dxx = r(0.003)
        for _, feat in iter_features(rec):
            if feat.location is None:
                continue
            if feat.location.strand:  # not 0 or None
                start, end = feat.location.start, feat.location.end
                angle1 = get_angle(start)  # type: ignore
                angle2 = get_angle(end)  # type: ignore
                w = dr if feat.location.strand > 0 else -dr
                a = arc(r0 + w, angle2, angle1, stroke="black", stroke_width=dxx)
                g.append(a)
    g2.append(g)


def draw_inverted_repeat(
    ir: IRScanResult,
    r0: float,
    r: Callable[[float], float],
    fs: float,
    get_angle: Callable[[int], float],
    scan_width: float,
    swapped: bool = False,
    **attrib: Any,
) -> Element:
    g = group(klass="irscan")

    sw = scan_width
    attrib.setdefault("stroke_width", sw)
    attrib.setdefault("stroke", "black")
    attrib.setdefault("font_weight", "bold")

    def txt(s: int, e: int, txt: str) -> None:
        sa = get_angle(s)
        ea = get_angle(e)
        a = arc(r0, sa, ea, **attrib)

        attr = {
            "stroke_linecap": "square",
            **attrib,
            "stroke_width": attrib["stroke_width"] / 2,
        }

        sbar = radial_line(sa, r0, r0 + r(0.01), **attr)
        ebar = radial_line(ea, r0, r0 + r(0.01), **attr)
        g.extend([a, sbar, ebar])
        g.append(
            text_horz(
                [Overlap(middle(sa, ea), txt, delta_a=0, delta_r=20)],
                r0,
                fs=fs * 1.5,
                r=r,
                **attrib,
            ),
        )

    txt(ir.IRA_start, ir.IRA_end, "IRB" if swapped else "IRA")
    txt(ir.IRB_start, ir.IRB_end, "IRA" if swapped else "IRB")
    return g


def show_overlap(
    g: Element,
    rec: SeqRecord,
    r0: float,
    r: Callable[[float], float],
    get_angle: Callable[[int], float],
    dr: float = 5,
    **attrib: Any,
) -> None:
    g2 = group(klass=klass(attrib, "overlap"))
    attrib.setdefault("fill", "red")
    spread = 0.1  # degrees
    overlaps = Overlapping(rec)

    dnamer = default_namer(rec.name)

    def isintron(part: SimpleLocation) -> bool:
        return hasattr(part, "intron") and getattr(part, "intron")

    def namer(feat: SeqFeature, part: SimpleLocation) -> str:
        if isintron(part):
            return "intron"
        return dnamer(feat.qualifiers)

    def draw_overlap(o: IOverlap) -> None:
        ir = o.overlap
        n1, n2 = namer(o.feat1, o.part1), namer(o.feat2, o.part2)
        s = f"{n1}-{n2} overlap [{ir.start + 1}..{ir.end}]"
        attr = {"data-info": s, **attrib}
        start, end = ir.start, ir.end
        angle1 = get_angle(start)
        angle2 = get_angle(end)
        strd = int(o.strand)
        w = dr if strd and strd > 0 else -dr
        a = arc(r0 + w, angle2 - spread, angle1 + spread, dr=w, **attr)
        g2.append(a)

    for o in overlaps.overlaps():
        # ignore overlaps with introns
        if isintron(o.part1) or isintron(o.part2):
            continue
        draw_overlap(o)
    if len(g2) > 0:
        g.append(g2)


def default_namer(default: str) -> Callable[[dict[str, Any]], str]:
    from .config import TRY_NAMES

    nc = 0

    def try_name(q: dict[str, Any], /) -> str:
        nonlocal nc
        # gene from sff, Parent from .gff3
        for n in TRY_NAMES:
            if n in q:
                return q[n][0]
        nc += 1
        return f"{default}/{nc}"

    return try_name


def get_gene_names(
    rec: SeqRecord,
    r0: float,
    fs: float,
    get_angle: Callable[[float], float],
    subfeatures: bool = True,
    strand_key: Callable[[Any], Any] | None = None,
    ntries: int = 4,
    namer: Callable[[SeqFeature], str] | None = None,
) -> dict[Any, list[Overlap]]:
    def try_name(q: dict[str, Any], /) -> str:
        return "unknown"

    def mynamer(feat: SeqFeature) -> str:
        return try_name(feat.qualifiers)

    if namer is None:
        try_name = default_namer(rec.name)  # noqa: F811
        namer = mynamer

    return get_gene_names_from_features(
        [feat for _, feat in iter_features(rec)],
        r0,
        fs,
        get_angle,
        target_length=len(rec.seq if rec.seq else ""),
        subfeatures=subfeatures,
        strand_key=strand_key,
        ntries=ntries,
        namer=namer,
    )


def genome_middle(target_length: int, pos1: int, pos2: int) -> int:
    pos1 = genome_wrap(target_length, pos1)
    pos2 = genome_wrap(target_length, pos2)
    if abs(pos1 - pos2) < target_length / 2:
        return (pos1 + pos2) // 2
    return genome_wrap(target_length, (pos1 + pos2 - target_length) // 2)


def strandkey(loc: SimpleLocation | CompoundLocation) -> Any:
    return loc.strand


def get_gene_names_from_features(  # noqa: C901
    features: Sequence[SeqFeature],
    r0: float,
    fs: float,
    get_angle: Callable[[float], float],
    target_length: int,
    subfeatures: bool = True,
    strand_key: Callable[[SimpleLocation | CompoundLocation], Any] | None = None,
    ntries: int = 4,
    namer: Callable[[SeqFeature], str] | None = None,
) -> dict[Any, list[Overlap]]:
    if strand_key is None:
        strand_key = strandkey

    if namer is None:
        try_name = default_namer("unknown")
        namer = lambda feat: try_name(feat.qualifiers)  # noqa: E731

    def position(
        location: SimpleLocation,
        feat: SeqFeature,
        sf: bool,
    ) -> Iterator[tuple[int, SimpleLocation]]:
        start, end = location.start, location.end
        if not isinstance(start, int) or not isinstance(end, int):
            return
        # find a position for this name
        # if we a naming all features then
        # the middle is good
        if sf:
            yield genome_middle(target_length, start, end), location
            return
        # single part middle is good
        if len(feat.location.parts) == 1:  # type: ignore
            yield genome_middle(target_length, start, end), location
            return
        # if we are a small angle apart
        # middle is good
        a = abs(get_angle(start) - get_angle(end))
        if min(a, 360 - a) < 5:
            yield genome_middle(target_length, start, end), location
            return

        # not using subfeatures but this is just too
        # spread out
        if feat.location is not None:
            for loc in feat.location.parts:
                yield genome_middle(target_length, loc.start, loc.end), loc  # type: ignore

    data = defaultdict(list)

    for feat in features:
        if not feat.location:
            continue
        gene = namer(feat)
        sf = subfeatures
        if subfeatures:
            parts = feat.location.parts
        else:
            # ok this feature is all over the place
            # so break it up anyway
            if feat.location.strand in {None, 0}:
                parts = feat.location.parts
                sf = True  # pretend we are subfeature
            else:
                parts = [feat.location]
        for fno, location in enumerate(parts):
            # don't name introns
            if hasattr(location, "intron") and getattr(location, "intron"):
                if not subfeatures:
                    continue
            locs = list(position(location, feat, sf))  # type: ignore
            for fno1, (pos, loc) in enumerate(locs):
                angle = get_angle(pos)

                if subfeatures and len(parts) > 1:
                    name = f"{gene}/{fno + 1}"
                elif len(locs) > 1:  # no sub features but multiloc
                    name = f"{gene}/{fno1 + 1}"
                else:
                    name = gene
                if location:
                    data[strand_key(location)].append(Overlap(angle, name, 0, 0, loc))

    rfs = fs * 180 / (r0 * pi)

    def text_width(txt: str, _fs: float) -> float:
        return text_len(txt, fs)

    ret = {}
    for k in data:
        sdata = sorted(data[k], key=lambda o: o.angle)
        sdata = fix_text_overlap(sdata, rfs, ntries=ntries, text_width=text_width)
        ret[k] = sdata

    return ret


def merge_attr(
    attrib: dict[str, Any],
    arg_attr: dict[str, Any],
    attrib_prefix: str = "",
    **local_defaults: Any,
) -> dict[str, Any]:
    def get_target(d: dict[str, Any]) -> dict[str, Any]:
        if not attrib_prefix:
            return d
        dl = {}
        if attrib_prefix in d:
            dl.update(d[attrib_prefix])
        # look for 'target_fill'
        tt = attrib_prefix + "_"
        n = len(tt)
        for k in d:
            if k.startswith(tt):
                dl[k[n:]] = d[k]
        return dl

    if not attrib_prefix:
        return {**attrib, **local_defaults, **arg_attr}

    return {**attrib, **get_target(local_defaults), **get_target(arg_attr)}


def genome_wrap(genome_length: int, position: int) -> int:
    if 0 <= position < genome_length:
        return position

    while position >= genome_length:
        position -= genome_length

    while position < 0:
        position += genome_length

    return position


def rotate(ir: IRScanResult, target_length: int) -> tuple[int, bool]:
    n = ir.IRA_end - ir.IRA_start
    ira = genome_wrap(target_length, ir.IRA_start)
    irb = genome_wrap(target_length, ir.IRB_start)
    ir_swapped = False
    if ira > irb:
        ir_swapped = True
        ira, irb = irb, ira
    ssc = irb - (ira + n)
    # assert ssc >= 0, ssc
    if ssc < ira:
        return target_length - (irb + n), ir_swapped
    return target_length - (ira + n), not ir_swapped


class AttrMerger:
    def __init__(self, **attrib: Any):
        self.attrib = attrib

    def merge_attr(
        self,
        arg_attr: dict[str, Any],
        **local_defaults: Any,
    ) -> dict[str, Any]:
        attrib_prefix = local_defaults.pop("attrib_prefix", "")
        return merge_attr(
            self.attrib,
            arg_attr,
            attrib_prefix=attrib_prefix,
            **local_defaults,
        )


class BaseDrawArgs(TypedDict, total=False):
    radius: int
    bg: bool | str
    genome: str | None
    subfeatures: bool
    scale: float
    opacity: float
    irscan: bool | IRScan | None
    rotate_image: bool
    attrib: dict[str, Any]


class BaseDraw(AttrMerger, abc.ABC):
    styles_to_classes: bool = True
    gc_pos: float = 0.38
    gc_width: float = 0.075
    gc_background: str = GC_HIST_BACKGROUND
    show_id: bool = False
    class_prefix: str = "cls"
    add_legend: bool = True
    stroke: str = "grey"
    ir_swapped: bool = False
    genome_info_class: type[GenomeInfo] = GenomeInfo

    def __init__(
        self,
        radius: int = 1000,
        bg: bool | str = False,
        genome: str | None = None,
        subfeatures: bool = False,
        scale: float = 1.0,
        opacity: float = 1.0,
        irscan: bool | IRScan | None = True,
        rotate_image: bool = False,
        attrib: dict[str, Any] | None = None,
    ):
        super().__init__(**(attrib or {}))
        self._id = 0
        if irscan is True:
            irscan = IRScan()
        elif irscan is False:
            irscan = None

        self.radius = radius
        self.bg = bg
        self.subfeatures = subfeatures
        self.scale = scale
        self.opacity = opacity
        self.irscan = irscan
        self.genome = genome
        self.rotate_image = rotate_image
        self.rotate_image_angle = 0

        r = maker(radius * scale)
        # svg = create_svg(radius)
        svg = svge(2 * radius)
        if bg:
            svg.append(rect(0, 0, 2 * radius, fill=bg, opacity=opacity, klass="bg"))

        # effective font width in degrees....

        self.fs = scale * pi * radius / 180.0
        self.svg = svg
        self.r = r
        self.sw = r(1 / self.radius)
        self.irscan_width = r(6 / self.radius)
        self.reinit()

    def getid(self) -> str:
        self._id += 1
        return f"id{self._id}"

    def get_gc_pos(self) -> tuple[float, float]:
        r = self.r
        pos = r(self.gc_pos)
        gcw = r(self.gc_width)
        return pos, gcw

    def reinit(self) -> None:
        radius = self.radius
        T = f"translate({radius} {radius})"
        g = group(transform=T)
        self.g = g

    def save(self, svgname: str, pretty_print: bool = False) -> None:
        savesvg(self.svg, svgname, pretty_print=pretty_print)

    def get_rec1_str(
        self,
        rec1: SeqRecord,
    ) -> str | tuple[str, dict[str, str]]:
        return f"({rec1.id or ''})"

    def center_text_list(self, rec: SeqRecord) -> tuple[TEXT_TYPE, str]:
        return self.genome_info_class().center_text_list(rec)

    def add_bling(
        self,
        rec: SeqRecord,
        rec1: SeqRecord | None = None,
        **attrib: Any,
    ) -> None:
        genome = self.add_center_text(rec, rec1, **attrib)
        if self.add_legend:
            self.add_legend_colors(genome)

    def add_center_text(
        self,
        rec: SeqRecord,
        rec1: SeqRecord | None = None,
        **attrib: Any,
    ) -> str:
        width = 2 * self.radius
        svg = self.svg
        genome = self.genome

        def patch(center_text: TEXT_TYPE) -> None:
            if len(center_text) >= 2:
                v = (
                    str(center_text[1]),
                    dict(font_style="italic", font_weight="bold"),
                )
                center_text[0], center_text[1] = (
                    v,
                    str(center_text[0]) + " genome",
                )

        center_text, g = self.center_text_list(rec)
        if genome is None:
            genome = g

        patch(center_text)

        if self.show_id:
            center_text.append(f"({rec.id})")

        if rec1 and rec1.id != rec.id:
            s = self.get_rec1_str(rec1)
            if s:
                center_text.extend(["vs.", s])

        r = width / (2 * self.radius)

        if self.rotate_image_angle != 0:
            center_text.append(
                (
                    f"rotated {self.rotate_image_angle:,} bp",
                    {"font-size": str(25 * r), "color": "rgba(0,0,0,.6)"},
                ),
            )

        attrib = {
            "opacity": 1.0,
            "font_size": 45 * r,
            "color": "rgba(0,0,0,1.0)",
            **attrib,
        }

        add_center_text(
            svg,
            center_text,
            width=width,
            height=width,
            inline=True,
            y_offset=20 * r,
            **attrib,
        )
        return genome

    def add_legend_colors(self, genome: str) -> None:
        width = 2 * self.radius
        r = 1.0
        add_legend(
            self.svg,
            self.genome_info_class().legend_text(genome),
            font_size=25 * r,
            y_offset=40 * r,
            width=width,
            height=width,
            inline=True,
            box_size=20 * r,
        )

    def get_IR(self, rec: SeqRecord) -> IRScanResult | None:
        def fl(res: list[SeqFeature], watson: bool) -> Iterator[SimpleLocation]:
            for r in res:
                if r.location is not None:
                    for loc in r.location.parts:
                        if loc.strand is not None:
                            if watson:
                                if loc.strand > 0:
                                    yield loc
                            else:
                                if loc.strand < 0:
                                    yield loc

        res: list[SeqFeature] = [r for r in rec.features if r.type == "repeat_region"]
        if res:
            fwd = list(fl(res, watson=True))
            rev = list(fl(res, watson=False))
            if not fwd or not rev:
                # maybe from a Genbank file
                return None
            lfwd, lrev = fwd[0], rev[0]
            return IRScanResult(
                int(lfwd.start),  # type: ignore
                int(lfwd.end),  # type: ignore
                int(lrev.start),  # type: ignore
                int(lrev.end),  # type: ignore
            )
        return None

    def has_seq(self, rec: SeqRecord) -> bool:
        return has_seq_data(rec)

    def is_sff(self, rec: SeqRecord) -> bool:
        return getattr(rec, "is_sffseq", False)

    def get_tick_pos(self) -> float:
        r = self.r
        return r(1) - r(0.15)

    def base_ticks(
        self,
        nbp: int,
        get_angle: Callable[[float], float],
        **attrib: Any,
    ) -> None:
        r, g = self.r, self.g
        loc = self.get_tick_pos()
        stroke = self.stroke

        attrib = self.merge_attr(
            attrib,
            opacity=1,
            fill="grey",
            grid_opacity=0.2,
            stroke=stroke,
        )
        g2 = ticks(
            nbp,
            loc,
            g,
            units=1000 if nbp > 100_000 else 500,
            grid=(r(0.35), r(0.75)),
            font_size=r(0.03),
            r=r,
            get_angle=get_angle,
            **attrib,
        )
        g2.append(
            circle(0, 0, loc, fill=None, stroke=stroke, stroke_width=0.5 * self.sw),
        )

    def ir_pos_with(self, has_seq: bool) -> float:
        pos, w = self.get_gc_pos()
        return pos + (w if has_seq else self.irscan_width)

    def ir_pos(self) -> float:
        # position of inverted repeat
        raise NotImplementedError()

    def base_doirscan(
        self,
        rec: SeqRecord,
        get_angle: Callable[[int], float],
        **attrib: Any,
    ) -> IRScanResult | None:
        g, r, fs = self.g, self.r, self.fs
        attrib = self.merge_attr(attrib)
        has_seq = self.has_seq(rec)

        ir = self.get_IR(rec)
        if ir is None and self.irscan and has_seq:
            ir = self.irscan.doirscan(str(rec.id), rec.seq)
        if ir is None:
            return None

        g.append(
            draw_inverted_repeat(
                ir,
                r0=self.ir_pos(),
                r=r,
                fs=fs,
                get_angle=get_angle,
                scan_width=self.irscan_width,
                swapped=self.ir_swapped,
                **attrib,
            ),
        )
        return ir

    def do_coverage(
        self,
        rec: SeqRecord,
        get_angle: Callable[[int], float],
        histogram_colors: dict[str, str],
        **attrib: Any,
    ) -> None:
        if not self.is_sff(rec):
            return
        sffrec = rec  # cast(SFFSeqRecord, rec)
        sw = 1000 / self.radius
        pos, gcw = self.get_gc_pos()
        height = gcw / 2

        r = self.r
        g = group(klass="sff")
        cattr = dict(fill=None, stroke_width=sw, stroke="grey")
        g.append(circle(0, 0, pos, **cattr))
        attrib = self.merge_attr(
            attrib,
            background=self.gc_background,  # , opacity=self.opacity
        )
        attrib.pop("fill", None)
        depth_histogram(
            sffrec,
            g,
            r0=pos - height,
            r=maker(height),
            get_angle=get_angle,
            fill=histogram_colors["depth"],
            **attrib,
            # opacity=self.opacity,
        )
        g.append(circle(0, 0, r(pos - height), **cattr))
        coverage_histogram(
            sffrec,
            g,
            r0=pos - 2 * height,
            r=maker(height),
            get_angle=get_angle,
            fill=histogram_colors["coverage"],
            **attrib,
            # opacity=self.opacity,
        )
        self.g.append(g)

    def do_gc_histogram(
        self,
        rec: SeqRecord,
        get_angle: Callable[[int], float],
        **attrib: Any,
    ) -> None:
        if not self.has_seq(rec):
            return
        attrib = self.merge_attr(
            attrib,
            background=self.gc_background,
            fill="grey",
            opacity=self.opacity,
        )
        g = self.g
        pos, gcw = self.get_gc_pos()
        gc_histogram(
            rec,
            g,
            r0=pos,
            r=maker(gcw),
            get_angle=get_angle,
            **attrib,
        )

    @abc.abstractmethod
    def draw(self, **attrib: Any) -> Element:
        raise NotImplementedError()

    def to_string(self, pretty_print: bool = False, **attrib: Any) -> str:
        return tostring(
            self.draw(**attrib),
            pretty_print=pretty_print,
        )

    def to_bytes(
        self,
        pretty_print: bool = False,
        **attrib: Any,
    ) -> bytes:
        return tobytes(
            self.draw(**attrib),
            pretty_print=pretty_print,
        )


class OGDraw(BaseDraw):
    show_span = False
    class_prefix = "ogdraw"

    def __init__(
        self,
        rec: SeqRecord,
        show_ticks: bool = True,
        **kwargs: Unpack[BaseDrawArgs],
    ):
        super().__init__(**kwargs)
        N = len(rec)
        ir_swapped = False
        if self.rotate_image:
            ir = self.get_IR(rec)
            rot, ir_swapped = rotate(ir, N) if ir else (0, False)
        else:
            rot = 0

        def get_angle(pos: float) -> float:
            return (N - pos - rot) * 360 / N

        def get_tick_angle(pos: float) -> float:
            return (N - pos) * 360 / N

        self.genome = self.genome or self.genome_info_class().genome(rec)
        self.rec = rec
        self.get_angle = get_angle
        self.get_tick_angle = get_tick_angle
        self.show_ticks = show_ticks
        self.rotate_image_angle = rot
        self.ir_swapped = ir_swapped

    def draw_rec(self, **attrib: Any) -> None:
        g, r, get_angle, rec, fs = self.g, self.r, self.get_angle, self.rec, self.fs

        try_name = default_namer(rec.name)

        def namer(feat: SeqFeature) -> str:
            g = try_name(feat.qualifiers)
            if feat.location:
                for loc in feat.location.parts:
                    if hasattr(loc, "comment") and getattr(loc, "comment"):
                        return g + "*"  # f" ({loc.comment})"
            return g

        r0 = r(1) - r(0.2)
        band_width = r(0.05)
        sw = r(1 / self.radius)

        data = get_gene_names(
            rec,
            r0,
            fs,
            get_angle,
            subfeatures=self.subfeatures,
            namer=namer,
        )
        # offset and dp are relative
        tattr = self.merge_attr(
            attrib,
            stroke_width=sw,
            offset=30 / 1000,
            dp=7 / 1000,
            opacity=1.0,
            klass="name",
        )

        show_overlap(g, rec, r0 - r(0.1), r, get_angle, dr=r(0.05), **attrib)

        if 1 in data:
            g.append(
                text_perp(
                    data[1],
                    r0 - band_width,
                    fs,
                    r,
                    radius=self.radius,
                    outside=True,
                    **tattr,
                ),
            )
        if -1 in data:
            g.append(
                text_perp(
                    data[-1],
                    r0 - 3 * band_width,
                    fs,
                    r,
                    radius=self.radius,
                    outside=False,
                    **tattr,
                ),
            )

        g.append(
            circle(
                0,
                0,
                r0 - r(0.1),
                stroke="black",
                klass="band-path",
                fill=None,
                stroke_width=3 * sw,
            ),
        )

        battr = self.merge_attr(
            attrib,
            stroke_width=sw,
            stroke="black",
            offset=30 / 1000,
            dp=7 / 1000,
            opacity=1.0,
            klass="band",
        )

        show_band(
            g,
            rec,
            r0 - r(0.1),
            r,
            self.genome or "<unknown>",
            get_angle,
            dr=band_width,
            show_span=self.show_span,
            **battr,
        )

    def ir_pos(self) -> float:
        return self.ir_pos_with(self.has_seq(self.rec))

    def doirscan(self, **attrib: Any) -> None:
        self.base_doirscan(self.rec, self.get_angle, **attrib)

    def postscript(self, **attrib: Any) -> None:
        self.doirscan(**attrib)
        # put bling ontop
        self.svg.append(self.g)
        self.add_bling(self.rec, **attrib)

    def extra(self, **attrib: Any) -> None:
        pass

    def ticks(self, **attrib: Any) -> None:
        self.base_ticks(len(self.rec), self.get_tick_angle, **attrib)

    @override
    def draw(self, **attrib: Any) -> Element:
        if self.show_ticks:
            self.ticks(**attrib)
        self.draw_rec(**attrib)
        self.extra(**attrib)
        self.postscript(**attrib)
        if self.styles_to_classes:
            styles_to_classes(self.svg, self.class_prefix)
        return self.svg


class GCOGDraw(OGDraw):
    @override
    def extra(self, **attrib: Any) -> None:
        self.do_gc_histogram(self.rec, self.get_angle, **attrib)


class DepthDraw(GCOGDraw):
    histogram_colors = HISTOGRAM_COLORS

    @override
    def extra(self, **attrib: Any) -> None:
        super().extra(**attrib)
        self.do_coverage(self.rec, self.get_angle, self.histogram_colors)
