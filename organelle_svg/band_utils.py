from __future__ import annotations

from collections import defaultdict
from typing import TYPE_CHECKING

from Bio.Seq import UndefinedSequenceError
from sortedcontainers import SortedList

if TYPE_CHECKING:
    from typing import Any
    from typing import Callable
    from typing import Iterator
    from Bio.SeqRecord import SeqRecord
    from Bio.SeqFeature import SeqFeature, SimpleLocation


def has_seq(rec: SeqRecord) -> bool:
    def tryseq() -> bool:
        try:
            if rec.seq is None:
                return False
            _ = rec.seq[0]
            return True
        except UndefinedSequenceError:
            return False

    # if isinstance(rec.seq, UnknownSeq):
    #     return False
    return tryseq()


class intrange:
    __slots__ = ("start", "end")

    def __init__(self, start: int, end: int):
        # all ranges are open [start,end)
        assert end > start, (start, end)
        self.start = start
        self.end = end

    def __len__(self) -> int:
        return self.end - self.start

    def __eq__(self, other: Any) -> bool:
        if not isinstance(other, intrange):
            return False
        return other.start == self.start and other.end == self.end

    def __hash__(self) -> int:
        return hash((self.start, self.end))

    @property
    def center(self) -> int:
        return (self.start + self.end) // 2

    def overlaps(self, ir: intrange) -> bool:
        return (
            (ir.start < self.end and ir.start >= self.start)  # ir.start overlaps self
            or (ir.end > self.start and ir.end <= self.end)  # ir.end overlaps self
            or (ir.start < self.start and ir.end > self.end)  # ir spans self
        )

    def cdistance(self, ir: intrange) -> int:
        """Center to center distance."""
        return self.center - ir.center

    def distance(self, ir: intrange) -> int:
        if self.overlaps(ir) or self.adjacent(ir):
            return 0
        if self.end < ir.start:  # we are before ir
            return self.end - ir.start

        # assert self.start > ir.end
        return self.start - ir.end  # we are after ir

    def merge(self, ir: intrange) -> intrange:
        # assert self.overlaps(ir) or self.adjacent(ir)
        return intrange(min(self.start, ir.start), max(self.end, ir.end))

    def adjacent(self, ir: intrange) -> bool:
        """endpoint adjacent"""
        return self.end == ir.start or self.start == ir.end

    def chop(self, ir: intrange) -> Iterator[intrange]:
        assert self.overlaps(ir)
        # return parts of this range that *don't* overlap with ir
        if ir.end < self.end:
            yield intrange(ir.end, self.end)
        if ir.start > self.start:
            yield intrange(self.start, ir.start)

    def intersect(self, ir: intrange) -> intrange:
        assert self.overlaps(ir)
        if ir.start >= self.start and ir.start < self.end:
            return intrange(ir.start, min(ir.end, self.end))
        if ir.end <= self.end and ir.end > self.start:
            return intrange(max(self.start, ir.start), ir.end)
        if ir.start >= self.start and ir.end <= self.end:
            return ir
        return self

    def __repr__(self) -> str:
        return f"{self.start}..{self.end}"


class intrangeset:
    def __init__(self, args: Iterator[intrange] | None = None):
        # ss will contain the start value for all ranges stored in
        # m (keyed on start values)
        # we keep that state that all adjacent and overlapping ranges are merged
        # so that there are just a sorted list of non-overlapping
        # ranges that don't abut.
        self.ss = SortedList()
        self.m: dict[int, intrange] = {}
        if args is not None:
            for a in args:
                self.add(a)

    def __repr__(self) -> str:
        return " ".join(repr(self.m[s]) for s in self.ss)

    def remove(self, ir: intrange) -> None:
        overlapping = list(self.overlapping(ir))
        for o in overlapping:
            self.ss.remove(o.start)
            del self.m[o.start]
            for n in o.chop(ir):
                self.add(n)

    def add(self, ir: intrange) -> None:
        ss, m = self.ss, self.m
        ip = max(ss.bisect_left(ir.start) - 1, 0)
        merge = []
        for ist in ss.islice(ip):
            ir2 = m[ist]  # fetch actual range
            if ir2.start > ir.end:  # gone too far
                break
            if ir2.adjacent(ir) or ir2.overlaps(ir):
                merge.append(ir2)
        # all ranges either overlap ir or
        # are adjacent to it. so we can merge them
        # into one.
        for i in merge:
            ir = i.merge(ir)
            ss.remove(i.start)
            del m[i.start]  # not needed anymore

        ss.add(ir.start)
        if not isinstance(ir, intrange):
            ir = intrange(ir.start, ir.end)
        m[ir.start] = ir

    def __iter__(self) -> Iterator[intrange]:
        m = self.m
        for s in self.ss:
            yield m[s.start]

    def overlaps(self, ir: intrange) -> bool:
        # ss = self.ss
        # ip = max(ss.bisect_left(ir.start) - 1, 0)
        for _ in self.overlapping(ir):
            return True
        return False

    def overlapping(self, ir: intrange) -> Iterator[intrange]:
        ss, m = self.ss, self.m
        ip = max(ss.bisect_left(ir.start) - 1, 0)
        for ist in ss.islice(ip):
            ir2 = m[ist]
            if ir2.start >= ir.end:
                break
            if ir2.overlaps(ir):
                yield ir2

    def overlaps_mod(self, ir: intrange, mx: int) -> bool:
        if self.overlaps(ir):
            return True
        start = ir.start % mx
        if start == ir.start:
            return False
        i = intrange(start, start + len(ir))
        if self.overlaps(i):
            return True
        return False


def iter_features(rec: SeqRecord) -> Iterator[tuple[int, SeqFeature]]:
    from .config import TYPES

    idx = 0
    for f in rec.features:
        if f.type in TYPES:
            yield idx, f
            idx += 1
        elif f.type == "gene" and hasattr(f, "sub_features"):  # from BCBio.GFF
            for sf in getattr(f, "sub_features"):
                if sf.type in TYPES:
                    yield idx, sf
                idx += 1
        else:
            idx += 1


class overlapset:
    def __init__(self, start: int, end: int, size: int):
        self.start, self.end, self.size = start, end, size
        self.buckets = {i: SortedList() for i in range(start, end, size)}
        self.m: dict[int, set[intrange]] = defaultdict(set)

    def add(self, ir: intrange) -> None:
        start, end, size = self.start, self.end, self.size
        m, buckets = self.m, self.buckets
        s, e = size * (ir.start // size), size * (ir.end // size)
        if e < ir.end:
            e += size
        if e > end:
            e = end
        for i in range(s, e, size):
            buckets[i].add(ir.start)
            m[ir.start].add(ir)
        if ir.end > end:
            s, e = start, start + size * ((ir.end - end) // size)
            if e < start + (ir.end - end):
                e += size
            for i in range(s, e, size):
                buckets[i].add(ir.start)
                m[ir.start].add(ir)

    def overlaps(self, ir: intrange) -> list[intrange]:
        end, size = self.end, self.size
        m, buckets = self.m, self.buckets
        s, e = size * (ir.start // size), size * (ir.end // size)
        if e < ir.end:
            e += size
        if e > end:
            e = end

        def ii() -> Iterator[intrange]:
            for b in range(s, e, size):
                bk = buckets[b]
                iend = bk.bisect_left(ir.end)
                for istart in bk.islice(0, iend + 1):
                    for i in m[istart]:
                        if i.overlaps(ir):
                            yield i

        return list(set(ii()))


class Overlap:
    def __init__(
        self,
        strand: str,
        feat1: SeqFeature,
        part1: SimpleLocation,
        feat2: SeqFeature,
        part2: SimpleLocation,
        overlap: Any,
    ):
        self.strand = strand
        self.feat1 = feat1
        self.feat2 = feat2
        self.part1 = part1
        self.part2 = part2
        self.overlap = overlap


# maybe use https://pypi.org/project/intervaltree/
class Overlapping:
    def __init__(
        self,
        rec: SeqRecord,
        nblocks: int = 100,
        strands: Callable[[Any], Any] = lambda s: s,
    ) -> None:
        self.b, self.fmap, self.collisions = self._mkoverlap(
            rec,
            strands=strands,
            nblocks=nblocks,
        )
        self.strands = strands
        self.rec = rec

    def _mkoverlap(
        self,
        rec: SeqRecord,
        nblocks: int,
        strands: Callable[[Any], Any] = lambda s: s,
    ) -> tuple[
        dict[str, overlapset],
        dict[Any, dict[intrange, tuple[SeqFeature, SimpleLocation]]],
        list[
            tuple[
                tuple[SeqFeature, SimpleLocation],
                tuple[SeqFeature, SimpleLocation],
            ]
        ],
    ]:
        collisions: list[
            tuple[
                tuple[SeqFeature, SimpleLocation],
                tuple[SeqFeature, SimpleLocation],
            ]
        ] = []
        fmap: dict[Any, dict[intrange, tuple[SeqFeature, SimpleLocation]]] = {
            strands(s): {} for s in (1, -1, 0, None)
        }
        size = len(rec) // nblocks
        size = 1 if size == 0 else size
        b = {strands(s): overlapset(0, len(rec), size) for s in (1, -1, 0, None)}
        for _, feat in iter_features(rec):
            if feat.location is None:
                continue
            for part in feat.location.parts:
                start, end = int(part.start), int(part.end)  # type: ignore
                if end <= start:
                    continue
                while start < 0:
                    start += len(rec)
                    end += len(rec)
                s = strands(part.strand)
                i = intrange(start, end)

                b[s].add(i)

                if i in fmap[s]:
                    of, op = fmap[s][i]
                    collisions.append(((feat, part), (of, op)))
                fmap[s][i] = (feat, part)
        return b, fmap, collisions

    def overlaps(self) -> Iterator[Overlap]:
        seen = set()
        rec, strands, b, fmap = self.rec, self.strands, self.b, self.fmap
        for _, feat in iter_features(rec):
            if feat.location is None:
                continue
            for part in feat.location.parts:
                # end must be greater that start
                if part.end <= part.start:  # type: ignore
                    continue
                s = strands(part.strand)
                po = intrange(int(part.start), int(part.end))  # type: ignore
                for o in b[s].overlaps(po):
                    # ofeat, opart = fmap[s][o]
                    if not o == po and o not in seen:
                        ofeat, opart = fmap[s][o]
                        overlap = po.intersect(o)
                        # assert overlap, (po, o)

                        yield Overlap(s, feat, part, ofeat, opart, overlap)
                    seen.add(o)
