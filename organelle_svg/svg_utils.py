from __future__ import annotations

import re
from collections import defaultdict
from math import cos
from math import pi
from math import sin
from math import tan
from typing import TYPE_CHECKING
from typing import TypeAlias
from uuid import uuid4
from xml.etree.ElementTree import Element

from .band_utils import intrange
from .band_utils import intrangeset
from .config import SVG_DEFAULT_FONTS


if TYPE_CHECKING:
    from typing import Any
    from typing import Callable
    from typing import Iterable
    from .histograms import Hist

RGBA = re.compile(
    r"^rgba\(\s*([0-9]+)\s*,\s*([0-9]+)\s*,\s*([0-9]+)\s*,\s*([01]?\.[0-9]+)\)$",
)

# see M.2 https://www.w3.org/TR/SVG11/attindex.html#PresentationAttributes
VALID_STYLES = {
    "alignment-baseline",
    "baseline-shift",
    "clip-path",
    "clip-rule",
    "clip",
    "color-interpolation-filters",
    "color-interpolation",
    "color-profile",
    "color-rendering",
    "color",
    "cursor",
    "direction",
    "display",
    "dominant-baseline",
    "enable-background",
    "fill-opacity",
    "fill-rule",
    "fill",
    "filter",
    "flood-color",
    "flood-opacity",
    "font-family",
    "font-size-adjust",
    "font-size",
    "font-stretch",
    "font-style",
    "font-variant",
    "font-weight",
    "glyph-orientation-horizontal",
    "glyph-orientation-vertical",
    "image-rendering",
    "kerning",
    "letter-spacing",
    "lighting-color",
    "marker-end",
    "marker-mid",
    "marker-start",
    "mask",
    "opacity",
    "overflow",
    "pointer-events",
    "shape-rendering",
    "stop-color",
    "stop-opacity",
    "stroke-dasharray",
    "stroke-dashoffset",
    "stroke-linecap",
    "stroke-linejoin",
    "stroke-miterlimit",
    "stroke-opacity",
    "stroke-width",
    "stroke",
    "text-anchor",
    "text-decoration",
    "text-rendering",
    "unicode-bidi",
    "visibility",
    "word-spacing",
    "writing-mode",
}


FontFamily = SVG_DEFAULT_FONTS


def rgba(c: str, m: re.Match[str] | None = None) -> dict[str, str]:
    m = RGBA.match(c) if m is None else m
    if not m:
        return dict(fill=c)
    r, g, b, a = m.group(1, 2, 3, 4)
    color = (int(r), int(g), int(b))
    return tuple_color(color, a=float(a))


RGB_Type: TypeAlias = tuple[int, int, int]


def tuple_color(c: RGB_Type, a: float | None = None) -> dict[str, str]:
    r, g, b = c
    fill = f"#{r:2X}{g:2X}{b:2X}".replace(" ", "0")
    if a is None:
        return dict(fill=fill)

    return dict(fill=fill, opacity=str(a))


def kebab_case(s: str) -> str:
    if s == "klass":
        return "class"
    if "_" in s:
        s = "-".join(s.split("_"))
    return s


def tostr(v: Any) -> str | dict[str, str]:
    if v is None:
        return "none"
    if isinstance(v, float):
        s = f"{v:.4f}"
    else:
        s = str(v)
        # allow for rgba(r,g,b,a) type colors
        m = RGBA.match(s)
        if m:
            return rgba(s, m)
    return s


def add_units(k: str, v: str) -> str:
    if k == "font-size":  # font-size without units as a class doesn't work
        if not v.endswith(("px", "em", "rem")):
            v = f"{v}px"
    return v


def attr(d: dict[str, Any]) -> dict[str, str]:
    ret: dict[str, str] = {}
    style: dict[str, Any] = {}
    for k, v2 in d.items():
        if k == "klass" and v2 is None:
            continue
        if isinstance(v2, dict):
            a = attr(v2)
            ret.update(a)
            continue
        k, v = kebab_case(k), tostr(v2)
        if k in VALID_STYLES:
            tgt = style
        else:
            tgt = ret
        if isinstance(v, dict):
            tgt.update(v)
        else:
            tgt[k] = add_units(k, v)

    if style:
        s = ";".join(f"{k}:{v}" for k, v in sorted(style.items())) + ";"
        if "style" in ret:
            r = ret["style"]
            if not r.endswith(";"):
                r += ";"
            s = r + s

        ret["style"] = s
    return ret


def attr2(**d: Any) -> dict[str, str]:
    return attr(d)


def tostyle(d: dict[str, Any]) -> str:
    style = None
    ret = attr(d)
    if "style" in ret:
        style = ret.pop("style")
    sret = ";".join(f"{k}:{str(v)}" for k, v in sorted(ret.items()))
    if sret:
        sret += ";"
    if style is not None:
        sret += style
    return sret


def toattrs(d: dict[str, Any]) -> str:
    return " ".join(f'{k}="{str(v)}"' for k, v in sorted(attr(d).items()))


radiansconversion = pi / 180.0


def frompolar(a: float, r: float) -> tuple[float, float]:
    return r * cos(a * radiansconversion), r * sin(a * radiansconversion)


def aaf(
    x: float,
    y: float,
    r: float,
    startangle: float,
    endangle: float,
) -> tuple[float, float, float, float]:
    xs = x + r * cos(startangle * radiansconversion)
    ys = y + r * sin(startangle * radiansconversion)
    xe = x + r * cos(endangle * radiansconversion)
    ye = y + r * sin(endangle * radiansconversion)
    return xs, ys, xe, ye


def relative_arc(r: float, a0: float, a1: float, x: float = 0, y: float = 0) -> str:
    large_arc_flag = 1 if abs(a0 - a1) > 180 else 0
    d = 1 if a1 > a0 else 0
    xe0, ye0 = frompolar(a1, r)
    xe0 += x
    ye0 += y
    return "A {:.4f} {:.4f} 0 {:d} {:d} {:.4f} {:.4f}".format(
        r,
        r,
        large_arc_flag,
        d,
        xe0,
        ye0,
    )


def line_arc(r: float, a0: float, a1: float, x: float = 0, y: float = 0) -> str:
    xe0, ye0 = frompolar(a0, r)
    xe0 += x
    ye0 += y
    return f"M {xe0:.4f} {ye0:.4f} " + relative_arc(r, a0, a1, x, y)


def pie(
    x: float,
    y: float,
    r: float,
    startangle: float,
    endangle: float,
    dr: float,
) -> str:
    r1 = r + dr
    xs0, ys0, xe0, ye0 = aaf(x, y, r, startangle, endangle)
    xs1, ys1, xe1, ye1 = aaf(x, y, r1, startangle, endangle)
    # If we want to plot angles larger than 180 degrees we need this
    large_arc_flag = 1 if abs(endangle - startangle) > 180 else 0
    d = 1 if endangle > startangle else 0
    ret = [
        f"M {xs0:.4f} {ys0:.4f}",
        "A {:.4f} {:.4f} 0 {:d} {:d} {:.4f} {:.4f}".format(
            r,
            r,
            large_arc_flag,
            d,
            xe0,
            ye0,
        ),
        f"L {xe1:.4f} {ye1:.4f}",
        "A {:.4f} {:.4f} 0 {:d} {:d} {:.4f} {:.4f}".format(
            r1,
            r1,
            large_arc_flag,
            abs(d - 1),
            xs1,
            ys1,
        ),
        f"L {xs0:.4f} {ys0:.4f} Z",
    ]
    return " ".join(ret)


def middle(a0: float, a1: float) -> float:
    return (a0 + a1 + 360) / 2 - 360 if abs(a0 - a1) > 180 else (a0 + a1) / 2


# https://css-tricks.com/svg-path-syntax-illustrated-guide/


def ribbon_path(
    r0: float,
    s0: float,
    e0: float,
    r1: float,
    s1: float,
    e1: float,
    bezier_radius: float = 0,
    crest: float = 0,
) -> str:
    xs0, ys0, xe0, ye0 = aaf(0, 0, r0, s0, e0)
    xs1, ys1, xe1, ye1 = aaf(0, 0, r1, s1, e1)
    laf0 = 1 if abs(e0 - s0) > 180 else 0
    laf1 = 1 if abs(e1 - s1) > 180 else 0
    d0 = 1 if e0 > s0 else 0
    d1 = 1 if e1 > s1 else 0
    #     m0 = (e0 + s1 + 360) / 2 - 360 if abs(e0 - s1) > 180 else (e0 + s1) / 2
    #     m1 = (e1 + s0 + 360) / 2 - 360 if abs(e1 - s0) > 180 else (e1 + s0) / 2

    xm0, ym0 = frompolar(e0, bezier_radius)
    xm1, ym1 = frompolar(s1, bezier_radius)
    #     crestx, cresty = frompolar(m0, crest)
    #     xm0, ym0 = crestx + (xe0 - crestx)*bezier_radius/r0, cresty + (ye0 - cresty)*bezier_radius/r0
    #     xm1, ym1 = crestx + (xs1 - crestx)*bezier_radius/r1, cresty + (ys1 - cresty)*bezier_radius/r1

    ret = [
        f"M {xs0:.4f} {ys0:.4f}",
        f"A {r0:.4f} {r0:.4f} 0 {laf0:d} {d0:d} {xe0:.4f} {ye0:.4f}",
        "C {:.4f} {:.4f} {:.4f} {:.4f} {:.4f} {:.4f}".format(
            xm0,
            ym0,
            xm1,
            ym1,
            xs1,
            ys1,
        ),
    ]
    xm0, ym0 = frompolar(e1, bezier_radius)
    xm1, ym1 = frompolar(s0, bezier_radius)
    #     crestx, cresty = frompolar(m1, crest)
    #     xm0, ym0 = crestx + (xe1 - crestx)*bezier_radius/r0, cresty + (ye1 - cresty)*bezier_radius/r0
    #     xm1, ym1 = crestx + (xs0 - crestx)*bezier_radius/r1, cresty + (ys0 - cresty)*bezier_radius/r1
    ret.extend(
        [
            "A {:.4f} {:.4f} 0 {:d} {:d} {:.4f} {:.4f}".format(
                r1,
                r1,
                laf1,
                d1,
                xe1,
                ye1,
            ),
            "C {:.4f} {:.4f} {:.4f} {:.4f} {:.4f} {:.4f} Z".format(
                xm0,
                ym0,
                xm1,
                ym1,
                xs0,
                ys0,
            ),
        ],
    )
    return " ".join(ret)


def bezier_path(
    r0: float,
    a0: float,
    r1: float,
    a1: float,
    bezier_radius: float = 0,
    crest: float = 0,
) -> str:
    xs0, ys0 = frompolar(a0, r0)
    xs1, ys1 = frompolar(a1, r1)
    # laf = 1 if abs(a0 - a1) > 180 else 0
    # d0 = 1 if a0 > a1 else 0

    #     m = (a0 + a1 + 360) / 2 - 360 if abs(a0 - a1) > 180 else (a0 + a1) / 2

    xm0, ym0 = frompolar(a0, bezier_radius)
    xm1, ym1 = frompolar(a1, bezier_radius)
    #     crestx, cresty = frompolar(m, crest)
    #     xm0, ym0 = crestx + (xs0 - crestx)*bezier_radius/r0, cresty + (ys0 - cresty)*bezier_radius/r0
    #     xm1, ym1 = crestx + (xs1 - crestx)*bezier_radius/r1, cresty + (ys1 - cresty)*bezier_radius/r1
    ret = [
        f"M {xs0:.4f} {ys0:.4f}",
        "C {:.4f} {:.4f} {:.4f} {:.4f} {:.4f} {:.4f}".format(
            xm0,
            ym0,
            xm1,
            ym1,
            xs1,
            ys1,
        ),
    ]
    return " ".join(ret)


def annular(x: float, y: float, r: float, dr: float) -> str:
    r1 = r + dr
    xs0, ys0, xe0, ye0 = aaf(x, y, r, 0, 359.999)
    xs1, ys1, xe1, ye1 = aaf(x, y, r1, 0, 359.999)

    ret = [
        f"M {xs0:.4f} {ys0:.4f}",
        f"A {r:.4f} {r:.4f} 0 1 1 {xe0:.4f} {ye0:.4f}z",
        f"M {xe1:.4f} {ye1:.4f}",
        f"A {r1:.4f} {r1:.4f} 0 1 0 {xs1:.4f} {ys1:.4f}Z",
    ]
    return " ".join(ret)


def annular_path(x: float, y: float, r: float, dr: float) -> str:
    # simpler without sin cosine comp
    o = r + dr
    return f"""M {x:.4f},{y - o:.4f} A{o:.4f},{o:.4f} 0 1,1 {x - 0.001:.4f},{y - o:.4f} z
               M {x:.4f},{y + r:.4f} A{r:.4f},{r:.4f} 0 1,0 {x - 0.001:.4f},{y + r:.4f} z"""


def svge(width: int, height: int | None = None, **attrib: Any) -> Element:
    if height is None:
        height = width
    return Element(
        "svg",
        attrib=attr2(
            viewBox=f"0 0 {width} {height}",
            xmlns="http://www.w3.org/2000/svg",
            **attrib,
        ),
    )


def rect(
    x: float,
    y: float,
    width: int,
    height: int | None = None,
    **attrib: Any,
) -> Element:
    if height is None:
        height = width
    return Element(
        "rect",
        attrib=attr2(x=x, y=y, width=width, height=height, **attrib),
    )


def circle(cx: float, cy: float, r: float, **attrib: Any) -> Element:
    return Element("circle", attrib=attr2(cx=cx, cy=cy, r=r, **attrib))


def line(
    x1: float,
    y1: float,
    x2: float,
    y2: float,
    **attrib: Any,
) -> Element:
    return Element("line", attrib=attr2(x1=x1, x2=x2, y1=y1, y2=y2, **attrib))


def radial_line(a: float, r0: float, r1: float, **attrib: Any) -> Element:
    x1, y1 = frompolar(a, r0)
    x2, y2 = frompolar(a, r1)
    return line(x1, y1, x2, y2, **attrib)


def text(txt: str, **attrib: Any) -> Element:
    attrib.setdefault("font_family", FontFamily)
    e = Element("text", attrib=attr(attrib))
    e.text = txt
    return e


def group(**attrib: Any) -> Element:
    return Element("g", attrib=attr(attrib))


def path(d: str, **attrib: Any) -> Element:
    return Element("path", attrib=attr2(d=d, **attrib))


def arc(
    r: float,
    sa: float,
    ea: float,
    dr: float = 0,
    **attrib: Any,
) -> Element:
    return pos_arc(r, sa, ea, dr=dr, x=0, y=0, **attrib)


def pos_arc(
    r: float,
    sa: float,
    ea: float,
    dr: float = 0,
    x: float = 0,
    y: float = 0,
    **attrib: Any,
) -> Element:
    if dr == 0:
        attrib["d"] = line_arc(r, sa, ea, x, y)
        attrib["fill"] = None
    else:
        attrib["d"] = pie(x, y, r, sa, ea, dr=dr)
    return Element("path", attrib=attr(attrib))


def text_arc(
    texts: str,
    sa: float,
    r0: float,
    pthid: str | None = None,
    maxext: float = 359.99,
    **attrib: Any,
) -> Element:
    # sa is the start angle....
    # can't seem to center this text
    # pthid must be a unique id
    g, txtf = radial_text(r0, sa=0, pthid=pthid, maxext=maxext)
    txtf(texts, sa, attrib)
    return g


def text_path(txt: str, pthid: str, **attrib: Any) -> Element:
    tp = Element("textPath", attrib={"href": "#" + pthid})
    tp.text = txt
    te = Element("text", attrib=attr(attrib))
    te.append(tp)
    return te


# e.g.
# g, txt = radial_text(self.r(0.8), font_size=30)
# for i in range(0, 360, 10):
#     txt(f'text at {i}deg', i)
# self.g.append(g)
def radial_text(
    r0: float,
    sa: float = 0,
    pthid: str | None = None,
    offset: float = 0,
    maxext: float = 359.99,
    **attrib: Any,
) -> tuple[Element, Callable[[str, float, dict[str, Any]], None]]:
    rpthid = f"id{uuid4().hex}" if pthid is None else pthid

    attrib.setdefault("font_family", FontFamily)
    g = group(**attrib)
    pth = arc(r0, sa, sa + maxext, id=rpthid)
    g.append(pth)

    def tfunc(txt: str, sa: float, attrib: dict[str, Any]) -> None:
        attrib.setdefault("transform", f"rotate({sa + offset})")
        g.append(text_path(txt, rpthid, **attrib))

    return g, tfunc


def maker(radius: float) -> Callable[[float], float]:
    def r(v: float) -> float:
        return v * radius

    return r


def klass(attrib: dict[str, Any], default: str | None = None) -> str | None:
    return attrib.pop("klass") if "klass" in attrib else default


def hist(
    data: Iterable[Hist],  # actually Hist(start,end,value)
    pos: float,
    r: Callable[[float], float],
    get_angle: Callable[[int], float],
    offset: bool = True,
    background: bool = False,
    inverted: bool = False,
    rule: Callable[[Hist], tuple[float, dict[str, Any] | None]] | None = None,
    **attrib: Any,
) -> Element:
    def myrule(d: Hist) -> tuple[float, dict[str, Any] | None]:
        return d.value, None

    if rule is None:
        rule = myrule

    g = group(klass=klass(attrib, None))
    if background:
        g.append(path(d=annular_path(0, 0, pos, r(1)), fill=background, opacity=0.6))
    rr = pos + (r(0.5) if offset else 0)
    if offset:
        g.append(circle(0, 0, rr, fill=None, stroke="white", stroke_width=r(0.02)))
    top = r(1)
    for d in data:
        v, a = rule(d)
        s, e = get_angle(d.start), get_angle(d.end)
        attrs = ({**attrib, **a}) if a else attrib
        v = r(v) / 2 if offset else r(v)
        if inverted:
            a1 = arc(rr + top, s, e, -v, **attrs)
        else:
            a1 = arc(rr, s, e, v, **attrs)
        g.append(a1)
    return g


def lines(
    data: Iterable[tuple[int, int, float]],
    inner_r: float,
    r: Callable[[float], float],
    get_angle: Callable[[float], float] = lambda a: a,
    inverted: bool = False,
    background: str | None = None,
    **attrib: Any,
) -> Element:
    # data is [(angle1, angle2, radius)] where angle is interpreted by get_angle
    stroke = attrib.get("stroke", "white")
    sw = attrib.get("stroke_width", r(0.2))

    one = r(1)

    g = group(klass=klass(attrib, None))
    if background:
        g.append(
            path(
                d=annular_path(0, 0, inner_r, one),
                fill=background,
                opacity=0.1,
                stroke_width=0,
            ),
        )

    p = inner_r
    end = inner_r + one
    dp = r(0.1)
    # draw annular circles from inner_r to inner_r + one
    while p < end:
        g.append(circle(0, 0, p, fill=None, stroke=stroke, stroke_width=sw))
        p += dp
    llines = []
    for s, e, radius in data:
        a1 = get_angle(s)
        a2 = get_angle(e)
        a = middle(a1, a2)
        if inverted:
            x, y = frompolar(a, inner_r + one - r(radius))
        else:
            x, y = frompolar(a, inner_r + r(radius))
        llines.append((x, y))
    llines.append(llines[0])
    pth = (
        [f"M{x},{y}" for x, y in llines[:1]]
        + [f"L{x},{y}" for x, y in llines[1:]]
        + ["M{},{}".format(*frompolar(360, inner_r + (one if inverted else 0)))]
    )

    pth.append(relative_arc(inner_r + (one if inverted else 0), 360, 0) + " Z")
    pths = " ".join(pth)
    ad = {"stroke_width": sw, **attrib}
    g.append(path(d=pths, **ad))
    return g


def text_len(txt: str, fs: float) -> float:
    # Just guessing
    n = len(txt)
    ni = txt.count("i") + txt.count("l") + txt.count("j")
    nm = txt.count("m")
    dx = fs / 2
    return (n - ni - nm) * dx + ni * (dx / 2) + nm * (dx * 1.4)


def extract_attrs(attrib: dict[str, Any], prefix: str) -> dict[str, Any]:
    if not prefix.endswith("_"):
        prefix += "_"
    n = len(prefix)
    return {k[n:]: v for k, v in attrib.items() if k.startswith(prefix)}


def exclude_attrs(attrib: dict[str, Any], *prefixes: str) -> dict[str, Any]:
    prefixes = tuple({p if p.endswith("_") else p + "_" for p in prefixes})
    return {k: v for k, v in attrib.items() if not k.startswith(prefixes)}


def text_perp(
    data: Iterable[Overlap],
    r0: float,
    fs: float,
    r: Callable[[float], float],
    outside: bool = True,
    offset: float = 0.0,
    dp: float = 0.02,
    radius: float = 1000,
    **attrib: Any,
) -> Element:
    """Draw polar text at radius r0 of fontsize fs.

    :param data: a list of [Overlap(angle, text, delta_angle, delta_r)]
    """
    g = group(klass=klass(attrib, None))
    tattrib = {k: v for k, v in attrib.items() if k not in {"stroke", "stroke_width"}}
    if "font_family" not in tattrib:
        tattrib["font_family"] = FontFamily
    sw = r(3 / radius)
    if "stroke_width" in attrib:
        sw = attrib.pop("stroke_width")
    if "fill" in attrib:
        attrib.pop("fill")
    if "stroke" not in attrib:
        attrib["stroke"] = "grey"
    if outside:
        dp = abs(r(dp))
    else:
        dp = -abs(r(dp))
    offset = r(offset)

    for overlap in data:
        # txt = f'{a}'

        r_angle = angle_wrap(overlap.angle + overlap.delta_a)

        # tlen = text_len(overlap.text, fs)

        if abs(r_angle) > 90 and abs(r_angle) < 270:
            rot = 180
        else:
            rot = 0
        dx, dy = fs / 2, fs / 3

        if outside:
            text_anchor = "start" if rot == 0 else "end"
            delta_r = abs(overlap.delta_r)
            r_0 = r0 + delta_r + offset
        else:
            text_anchor = "end" if rot == 0 else "start"
            dx = -dx
            delta_r = -abs(overlap.delta_r)
            r_0 = r0 + delta_r - offset

        t1 = text(
            overlap.text,
            x=0,
            y=0,
            dx=dx if rot == 0 else -dx,
            dy=dy,
            font_size=fs,
            text_anchor=text_anchor,
            **tattrib,
        )
        tr1 = f"rotate({r_angle:.4f}) translate({r_0:.4f} 0) rotate({rot})"
        t1.attrib["transform"] = tr1

        g.append(t1)
        delta_a = angle_wrap(overlap.delta_a)
        x, y = frompolar(delta_a, r_0 + dx)
        x -= r0

        aa = angle_wrap(overlap.angle)
        tr = f"rotate({aa:.4f}) translate({r0:.4f} 0)"
        if delta_a == 0:
            p = line(
                dp / 2,
                0,
                x,
                y,
                stroke_width=sw,
                fill=None,
                transform=tr,
                **attrib,
            )
        else:
            dy = tan(delta_a * pi / 180.0) * dp
            p = path(
                d=f"M0,0 L{dp:.4f},0 L{x - dp:.4f},{y - dy:.4f} L{x:.4f},{y:.4f}",
                stroke_width=sw,
                fill=None,
                transform=tr,
                **attrib,
            )
        g.append(p)
    return g


def angle_wrap(a: float) -> float:
    while a >= 360:
        a -= 360
    while a < 0:
        a += 360
    return a


def text_horz(
    data: list[Overlap],
    r0: float,
    fs: float,
    r: Callable[[float], float],
    dp: float = 0.02,
    **attrib: Any,
) -> Element:
    """Draw horizonal text at radius r0 of fontsize fs.

    :param data: a list of [Overlap(angle, text, delta_angle, delta_r)]
    """
    tattrib = {k: v for k, v in attrib.items() if k not in {"stroke", "stroke_width"}}
    attrib = {
        k: v for k, v in attrib.items() if k not in {"font_family", "font_weight"}
    }
    if "font_family" not in tattrib:
        tattrib["font_family"] = FontFamily
    if "fill" in attrib:
        attrib.pop("fill")
    if "stroke" not in attrib:
        attrib["stroke"] = "grey"
    dp = r(dp)
    g2 = group(klass=klass(attrib, None))

    for overlap in data:
        delta_r, delta_a, angle = (
            overlap.delta_r,
            angle_wrap(overlap.delta_a),
            angle_wrap(overlap.angle),
        )

        rot = angle + 90
        tr = f"rotate({rot + delta_a:.4f}) translate(0 -{r0 + delta_r:.4f})"
        if 0 < angle < 180:
            tr += " rotate(180)"
        t_ = text(
            overlap.text,
            x=0,
            y=0,
            dy=fs / 3,
            font_size=fs,
            text_anchor="middle",
            transform=tr,
            **tattrib,
        )
        g2.append(t_)

        # no need for ticks if we are close
        if abs(delta_r) < dp and abs(delta_a) < dp:
            continue
        dpp = dp if delta_r > 0 else -dp
        tr = f"rotate({angle:.4f}) translate({r0:.4f}, 0)"
        if delta_a == 0:
            p = line(0, 0, delta_r - dpp, 0, transform=tr, fill=None, **attrib)
        else:
            x, y = frompolar(delta_a, r0 + delta_r - dpp)
            x -= r0
            dy = tan(delta_a * pi / 180.0) * dp
            x1 = max(x - dpp, dpp) if delta_r > 0 else min(x - dpp, dpp)
            p = path(
                d=f"M{dpp / 2:.4f},0 L{dpp:.4f},0 L{x1:.4f},{y - dy:.4f} L{x:.4f},{y:.4f}",
                fill=None,
                transform=tr,
                **attrib,
            )
        g2.append(p)
    return g2


class Overlap:
    def __init__(
        self,
        angle: float,
        txt: str,
        delta_a: float = 0,
        delta_r: float = 0,
        props: Any = None,
    ):
        self.angle = angle
        self.text = txt
        self.delta_a = delta_a
        self.delta_r = delta_r
        self.props = props

    def __repr__(self) -> str:
        return f"Overlap({self.angle:3f}deg,{self.text}"


def fix_text_overlap(
    data: list[Overlap],
    fs: float,
    ntries: int = 3,
    delta_r_fudge: float = 1.5,
    delta_a_fudge: float = 1.0,
    text_width: Callable[[str, float], float] = text_len,
    outside: bool = True,
    text_height: Callable[[str, float], float] = lambda txt, fs: fs,
    level: int = 0,
) -> list[Overlap]:
    # currently fs is the angular width of the font i.e.
    # the number of *degrees* that it takes up at the
    # radius we a specifiying
    # assumes data is sorted by angle
    scale = 10.0
    irs = intrangeset()

    def irf(angle: float, txt: str) -> intrange:
        w = text_height(txt, fs)
        return intrange(int(scale * (angle - w / 2)), int(scale * (angle + w / 2)))

    # FIXME: doesn't deal with circularity of genome
    ndata = []
    unplaced = []
    for idx, o in enumerate(data):
        # warning angle should never change it is
        # the "correct" position.
        a = o.angle + o.delta_a
        i = irf(a, o.text)
        if irs.overlaps(i):
            # try for some space near us
            rng = i.end - i.start
            # try [-1,1,-2,2,-3,3] etc
            # for ii in [iii for j in range(1, ntries + 1) for iii in [-j, j]]:
            # only downstream
            for iii in range(1, ntries + 1):
                s = i.start + rng * iii
                i2 = intrange(s, s + rng)
                if not irs.overlaps(i2):
                    o.delta_a = angle_wrap(
                        o.delta_a + iii * text_height(o.text, fs) * delta_a_fudge,
                    )
                    i = i2
                    break
            else:
                unplaced.append((idx, o))
                continue
        irs.add(i)
        ndata.append(o)

    if unplaced:
        nndata = []
        for idx, o in unplaced:
            mx: float | None = None
            ii: int | None = None
            # try to find a spot with the shortest name
            for j in range(max(0, idx - ntries), min(len(data), idx + ntries)):
                oo = data[j]
                tw = text_width(oo.text, fs)
                if mx is None or tw < mx:
                    mx = tw
                    ii = j
            if ii is not None and mx is not None:
                oo = data[ii]
                dr = mx * delta_r_fudge
                o.delta_r += dr if outside else -dr
                o.delta_a = angle_wrap(
                    (oo.angle + oo.delta_a) - o.angle,
                )  # so o.angle + o.delta_a == oo.angle + oo.delta_a

            else:
                # what to do? do we every get here?
                dr = text_width(o.text, fs)
                o.delta_r += dr if outside else -dr

            nndata.append(o)
            # FIXME: effective font-size changes if moving radially

        return ndata + fix_text_overlap(
            nndata,
            fs,
            ntries=ntries,
            delta_a_fudge=delta_a_fudge,
            delta_r_fudge=delta_r_fudge,
            text_width=text_width,
            text_height=text_height,
            outside=outside,
            level=level + 1,
        )
    return ndata


def ticks(
    nbases: int,
    loc: float,
    g: Element,
    r: Callable[[float], float],
    get_angle: Callable[[float], float],
    font_size: float = 0.0,
    grid: tuple[float, float] | None = None,
    units: int = 1000,
    grid_opacity: float = 1.0,
    tick_info: dict[str, Any] | None = None,
    **attrib: Any,
) -> Element:
    g2 = group(klass=klass(attrib, "ticks"))

    if get_angle is None:
        get_angle = lambda b: b * 360 / nbases  # noqa: E731
    if not font_size:
        font_size = r(0.04)
    dx = 0
    tattrib = exclude_attrs(attrib, "text_", "grid_")
    txtattrib = {**tattrib, **extract_attrs(attrib, "text_")}
    txtattrib.setdefault("font_family", FontFamily)
    gattrib = {
        **tattrib,
        **extract_attrs(attrib, "grid_"),
    }

    # end = loc

    gattrib["opacity"] = grid_opacity
    tick_info = tick_info or {}

    def grid_line(w: float) -> None:
        if grid:
            g2.append(
                line(
                    x1=grid[0],
                    y1=0,
                    x2=grid[1],
                    y2=0,
                    stroke_width=r(w),
                    transform=rotate,
                    **gattrib,
                ),
            )

    def tick(w: float, h: float) -> None:
        g2.append(
            line(
                x1=loc + dx,
                y1=0,
                x2=loc + r(h),
                y2=0,
                stroke_width=r(w),
                transform=rotate,
                **tattrib,
            ),
        )

    def dolabel(base: float, angle: float, pos: float) -> None:
        txt = f"{base // 1000}kB"
        if 90 < angle < 270:
            rot = 180
        else:
            rot = 0

        text_anchor = "start" if rot == 0 else "end"

        tr = rotate + f"translate({loc + r(pos):.4f} 0) rotate({rot})"
        fs = font_size
        g2.append(
            text(
                txt,
                x=0,
                y=0,
                dx=fs / 4 if rot == 0 else -fs / 4,
                dy=fs / 3,
                font_size=fs,
                text_anchor=text_anchor,
                transform=tr,
                stroke_width=0,
                **txtattrib,
            ),
        )

    def get_defs(t_: Any, w: float, h: float) -> tuple[float, float]:
        if t_ not in tick_info:
            return w, h
        ti = tick_info[t_]
        w = ti.get("width", w)
        h = ti.get("height", h)
        return w, h

    for b in range(0, nbases, units):
        a = get_angle(b)
        rotate = f"rotate({a:.4f})"

        if b % (10 * units) == 0:
            w, h = get_defs(10, 3 / 1000, 0.03)
            if abs(b - nbases) > units:
                dolabel(b, a, h)
            tick(w, h)
            grid_line(3 / 1000)
        elif b % (5 * units) == 0:
            w, h = get_defs(5, 3 / 1000, 0.02)
            tick(w, h)
            grid_line(2 / 1000)
        elif b % units == 0:
            w, h = get_defs(1, 2 / 1000, 0.01)
            tick(w, h)
            grid_line(1 / 1000)

    g.append(g2)
    return g2


def xpath(svg: Element, path: str) -> list[Element]:
    return svg.findall(path)


def styles_to_classes(svg: Element, prefix: str = "cls") -> None:
    """Gather all styles up and create classes."""
    from .svg import add_styles

    styled = defaultdict(list)
    for p in xpath(svg, ".//*[@style]"):
        style = p.attrib["style"]
        style = style.strip()
        if not style.endswith(";"):
            style += ";"
        styled[style].append(p)

    styles = []
    cnames = {}
    for i, style in enumerate(styled):
        cls = f"{prefix}-{i}"
        cnames[style] = cls
        styles.append((cls, style))

    # remove style and replace with classNames
    for style in styled:
        for e in styled[style]:
            if "style" in e.attrib:
                del e.attrib["style"]
            if "class" in e.attrib:
                e.attrib["class"] += " " + cnames[style]
            else:
                e.attrib["class"] = cnames[style]
    add_styles(svg, styles)
