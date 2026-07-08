from __future__ import annotations

import math
from typing import TYPE_CHECKING

from .band_utils import iter_features
from .svg_utils import hist

if TYPE_CHECKING:
    from typing import Any
    from typing import Callable
    from typing import Iterator
    from xml.etree.ElementTree import Element
    from Bio.SeqRecord import SeqRecord
    from Bio.SeqFeature import SeqFeature
    from Bio.SeqFeature import SimpleLocation


def merge_attr(attr: dict[str, Any], **defaults: Any) -> dict[str, Any]:
    return {**defaults, **attr}


class Hist:
    def __init__(
        self,
        start: int,
        end: int,
        value: float,
        data: SeqFeature | None = None,
    ):
        self.data = data
        self.start = start
        self.end = end
        self.value = value


def gc_hist(rec: SeqRecord, span: int = 100) -> Iterator[Hist]:
    if rec.seq is None:
        return
    seq = rec.seq.lower()
    n = len(seq)
    # span = n // span
    for i in range(0, n, span):
        r = seq[i : i + span]
        start = i  # + 1
        end = i + len(r)
        v = r.count("g") + r.count("c")
        v = v / len(r)
        yield Hist(start=start, end=end, value=v)


def depth_hist(sff_rec: SeqRecord) -> Iterator[Hist]:
    def gattr(part: SimpleLocation) -> float:
        if not hasattr(part, "depth"):
            return -0.1
        v = getattr(part, "depth")
        if math.isnan(v):
            return -0.1
        return v

    yield from sff_hist(sff_rec, gattr)


def coverage_hist(sff_rec: SeqRecord) -> Iterator[Hist]:
    def gattr(part: SimpleLocation) -> float:
        v = getattr(part, "coverage", math.nan)
        if math.isnan(v):
            return -0.1
        return v / 100.0

    yield from sff_hist(sff_rec, gattr)


def sff_hist(
    sff_rec: SeqRecord,
    gattr: Callable[[SimpleLocation], float],
) -> Iterator[Hist]:
    # if not isinstance(sff_rec, SFFSeqRecord):
    #     raise ValueError("not a SFF SeqRecord")
    s: int
    e: int
    for _, feat in iter_features(sff_rec):
        if feat.location is None:
            continue
        for part in feat.location.parts:
            v = gattr(part)  # type: ignore
            s, e = part.start, part.end  # type: ignore
            yield Hist(start=s, end=e, value=v, data=feat)


def feat_to_data(d: Hist) -> dict[str, str] | None:
    if d.data is None:
        return None
    feat: SeqFeature = d.data
    for k in ("gene", "Name"):
        if k in feat.qualifiers:
            v = feat.qualifiers[k][0]
            return {"data-info": f"{v}={d.value:.3f}"}

    return None


def error_rule(
    error: str,
) -> Callable[[Hist], tuple[float, dict[str, str] | None]]:
    errord = dict(fill=error)

    def rule(d: Hist) -> tuple[float, dict[str, str] | None]:
        v = d.value
        if v < 0:
            return 1.0, errord
        return v, None  # feat_to_data(d)

    return rule


def depth_histogram(
    rec: SeqRecord,
    g: Element,
    r0: float,
    r: Callable[[float], float],
    get_angle: Callable[[int], float],
    **attrib: Any,
) -> None:
    error = attrib.pop("error", "red")
    attrib.setdefault("inverted", False)
    attrib.setdefault("offset", False)

    attrib = merge_attr(attrib, klass="depth")
    g.append(
        hist(
            depth_hist(rec),
            r0,
            r=r,
            get_angle=get_angle,
            rule=error_rule(error),
            **attrib,
        ),
    )


def coverage_histogram(
    rec: SeqRecord,
    g: Element,
    r0: float,
    r: Callable[[float], float],
    get_angle: Callable[[int], float],
    **attrib: Any,
) -> None:
    error = attrib.pop("error", "red")
    attrib.setdefault("inverted", False)
    attrib.setdefault("offset", False)

    attrib = merge_attr(attrib, klass="coverage")
    g.append(
        hist(
            coverage_hist(rec),
            r0,
            r=r,
            get_angle=get_angle,
            rule=error_rule(error),
            **attrib,
        ),
    )


def gc_histogram(
    rec: SeqRecord,
    g: Element,
    r0: float,
    r: Callable[[float], float],
    get_angle: Callable[[int], float],
    **attrib: Any,
) -> None:
    attrib.setdefault("offset", False)
    attrib.setdefault("inverted", True)
    attrib = merge_attr(attrib, klass="gc")
    d = gc_hist(rec, span=100)
    g.append(hist(d, r0, r=r, get_angle=get_angle, **attrib))
