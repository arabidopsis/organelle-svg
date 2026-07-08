from __future__ import annotations

import gzip
from html import escape as htmlescape
from io import StringIO
from pathlib import Path
from typing import Any
from typing import cast
from typing import TYPE_CHECKING
from typing import TypeAlias
from xml.etree.ElementTree import Element
from xml.etree.ElementTree import indent
from xml.etree.ElementTree import parse as etree_parse
from xml.etree.ElementTree import tostring as etree_tostring

from .svg_utils import FontFamily
from .svg_utils import RGB_Type
from .svg_utils import toattrs
from .svg_utils import tostyle
from .svg_utils import tuple_color

if TYPE_CHECKING:
    from typing import IO


def getsvg(svgname: str) -> Element:
    def read(fp: IO[bytes] | gzip.GzipFile) -> Element:
        x = etree_parse(fp, parser=None)
        return x.getroot()

    if svgname.endswith(".svgz"):
        with gzip.open(svgname, "rb") as fp:
            return read(fp)
    else:
        with open(svgname, "rb") as fp2:
            return read(fp2)


def tostring(svg: Element, pretty_print: bool = False) -> str:
    if pretty_print:
        indent(svg)
    return etree_tostring(svg, encoding="unicode")


def tobytes(svg: Element, pretty_print: bool = False) -> bytes:
    if pretty_print:
        indent(svg)
    return etree_tostring(svg)


def savesvg(
    svg: bytes | Element,
    svgname: str | Path,
    pretty_print: bool = False,
) -> None:
    def write(fp: IO[bytes] | gzip.GzipFile) -> None:
        if isinstance(svg, bytes):
            fp.write(svg)
        else:
            assert isinstance(svg, Element)
            fp.write(tobytes(svg, pretty_print=pretty_print))

    svgname = Path(svgname)
    if svgname.name.endswith(".svgz"):
        with gzip.open(svgname, "wb") as fp:
            write(fp)
    else:
        with open(svgname, "wb") as fp2:
            write(fp2)


def add_legend(
    svg: Element,
    lines: list[tuple[str, Any]],
    font_size: float = 30,
    y_offset: float = 10,
    sep: int = 20,
    left: int = 10,
    opacity: float = 0.5,
    width: float = 2000,
    height: float = 2000,
    font_family: str | None = None,
    box_size: float | None = None,
    inline: bool = False,
    escape: bool = True,
    **attrib: Any,
) -> Element:
    font_family = font_family or FontFamily
    # nlines = len(lines)
    lines = list(reversed(lines))
    size = font_size // 1.5 if box_size is None else box_size
    g = create_element('<g class="legend"></g>')
    svg.append(g)
    rstyle = {"stroke-width": 2, "stroke": "black"}
    tstyle = {
        "opacity": opacity,
        "font-family": font_family,
        "font-size": f"{font_size}px",
    }
    tstyle = {**attrib, **tstyle}
    if inline:
        t1 = rstyle
        t2 = tstyle
    else:
        t1 = {}
        t2 = {}
    for i, (line, color) in enumerate(lines):
        if escape:
            line = htmlescape(line)
        y_pos = height - y_offset - i * font_size - sep
        if color:
            if isinstance(color, tuple):
                tx = {
                    **t1,
                    **tuple_color(
                        cast(RGB_Type, color[:3]),
                        color[4] if len(color) > 3 else None,
                    ),
                }
            else:
                tx = {**t1, "fill": color}
            txs = toattrs(tx)
            e = create_element(
                f"""<rect x="{left}" y="{y_pos - size}px"
                width="{size}" height="{size}"
                class="legend-box" {txs} />""",
            )
            g.append(e)
        txs = toattrs(t2) if t2 else ""
        e = create_element(
            f"""<text x="{left + 2 * size}" y="{y_pos}"
            {txs}
            text-anchor="left" class="legend-text">{line}</text>""",
        )
        g.append(e)
    if not inline:
        add_styles(
            svg,
            [("legend-text", tostyle(tstyle)), ("legend-box", tostyle(rstyle))],
        )
    return svg


TEXT_TYPE: TypeAlias = list[str | tuple[str, Any]]


def add_center_text(
    svg: Element,
    text: str | TEXT_TYPE,
    font_size: float = 60,
    y_offset: float = 0,
    width: float = 2000,
    height: float = 2000,
    sep: float = 20,
    opacity: float = 0.6,
    color: str = "#000000",
    font_family: str | None = None,
    inline: bool = False,
    escape: bool = True,
    **attrib: Any,
) -> Element:
    font_family = font_family or FontFamily
    if isinstance(text, str):
        lines: TEXT_TYPE = [t.strip() for t in text.splitlines() if t.strip()]
    else:
        lines = text
    nlines = len(lines)
    mid = nlines / 2
    w2 = int(width / 2)
    y_start = height / 2 + y_offset - mid * font_size

    style = {
        "opacity": opacity,
        "fill": color,
        "font-family": font_family,
        "font-size": f"{font_size}px",
    }
    style = {**attrib, **style}

    g = create_element('<g class="center-text"></g>')
    svg.append(g)
    y_pos = y_start
    for line in lines:
        t2 = style if inline else {}
        fs = font_size
        if not isinstance(line, str):
            line, s = line
            if s:
                fs = float(s.get("font-size", font_size))
                t2 = {**t2, **s}

        if escape:
            line = htmlescape(line)
        t2s = toattrs(t2) if t2 else ""
        e = create_element(
            f"""<text x="{w2}px" y="{int(y_pos)}px" class="center-text-line"
                text-anchor="middle" {t2s}
        >{line}</text>""",
        )
        g.append(e)
        y_pos += fs + sep

    if not inline:
        add_styles(svg, [("center-text-line", tostyle(style))])
    return svg


def viewbox(svg: Element) -> tuple[int, int]:
    if "width" in svg.attrib and "height" in svg.attrib:
        w = svg.attrib["width"]
        h = svg.attrib["height"]
        del svg.attrib["width"]
        del svg.attrib["height"]
        w = w.replace("px", "")
        h = h.replace("px", "")
        svg.attrib["viewBox"] = f"0 0 {w} {h}"
        return int(w), int(h)
    return 0, 0


STYLE_T = """.{cls} {{ {style} }}"""


def add_styles(svg: Element, styles: list[tuple[str, str]]) -> None:
    ll = [STYLE_T.format(cls=cls, style=style) for cls, style in styles]

    sep = "\n"
    style = f"""<style>{sep}{sep.join(ll)}{sep}</style>"""
    estyle = create_element(style)
    svg.insert(0, estyle)
    # svg.getchildren()[0].addprevious(style)


def create_element(txt: str) -> Element:
    return etree_parse(StringIO(txt), parser=None).getroot()


def get_dim(svg: Element) -> tuple[int, int]:
    w, h = svg.attrib["viewBox"].split()[-2:]
    return int(w), int(h)
