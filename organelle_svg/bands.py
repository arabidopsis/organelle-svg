from __future__ import annotations

from collections import Counter
from collections import defaultdict
from typing import TYPE_CHECKING

from Bio.SeqFeature import SimpleLocation

from .band_utils import intrange
from .band_utils import intrangeset
from .band_utils import iter_features

if TYPE_CHECKING:
    from typing import Any
    from typing import Callable
    from typing import Iterator
    from Bio.SeqRecord import SeqRecord
    from Bio.SeqFeature import SeqFeature


def sanitize_option(name: str) -> str:
    # can't use () for options since these are stripped
    # in `Circos::Utils::parse_csv` ....
    # 'LEFT DOUBLE PARENTHESIS' (U+2E28)
    # 'RIGHT DOUBLE PARENTHESIS' (U+2E29)
    name = name.replace(" ", "\u00a0")
    return name.replace("(", "⸨").replace(")", "⸩")


def sanitize_name(name: str) -> str:
    # maybe set tab as separator.... `file_delim = \t`
    # no spaces in names (circos files are space separated)
    # but 'NO-BREAK SPACE' (U+00A0) seems to work
    name = name.replace(" ", "\u00a0")
    # circos svg interprets underscores as subscript for svg
    #  'FULLWIDTH LOW LINE' (U+FF3F)
    name = name.replace("_", "\uff3f")
    # circos seems to barf with quote "'" marks for 3' 5' etc
    # 'RIGHT SINGLE QUOTATION MARK' (U+2019)
    name = name.replace("'", "’")
    return name


class Band:
    STRANDS = {1: "watson", -1: "crick", None: "none", 0: "unknown"}

    def __init__(
        self,
        name: str,
        gene: str,
        start: int,
        end: int,
        strand: int,
        idx: int,
        gene_type: str,
        color: str | None = None,
        properties: dict[str, Any] | None = None,
    ):
        self.name = name
        self.gene = gene
        self.start = start
        self.end = end
        self.strand = strand
        self.idx = idx
        self.type = gene_type
        self.color_ = color
        self.strand_key_: str | None = None
        self.option_name_: str | None = None
        self.option_gene_: str | None = None
        self.properties_ = properties or {}
        self._props_text: str | None = None

    def prop(self, prop: str, default: Any = None) -> Any:
        return self.properties_.get(prop, default)

    @property
    def properties(self) -> str:
        if self._props_text is None:
            self._props_text = ",".join(f"{k}={v}" for k, v in self.properties_.items())
        return self._props_text

    @property
    def option_name(self) -> str:
        # can't use () for options since these are stripped
        # in `Circos::Utils::parse_csv` ....
        if self.option_name_ is None:
            self.option_name_ = sanitize_option(self.name)
        return self.option_name_

    @property
    def option_gene(self) -> str:
        if self.option_gene_ is None:
            self.option_gene_ = sanitize_option(self.gene)
        return self.option_gene_

    @property
    def color(self) -> str | None:
        return self.color_

    def clone(self, start: int, end: int) -> Band:
        ret = self.__class__(
            self.name,
            self.gene,
            start,
            end,
            self.strand,
            self.idx,
            self.type,
            self.color_,
        )
        # ret.color_ = self.color_
        return ret

    def __repr__(self) -> str:
        return f"{self.__class__.__name__}({self.name},{self.gene},[{self.start},{self.end}], {self.strand})"

    @property
    def strand_key(self) -> str:
        if self.strand_key_ is None:
            self.strand_key_ = self.STRANDS.get(self.strand)
        assert isinstance(self.strand_key_, str)
        return self.strand_key_


def get_range_map(rec: SeqRecord) -> dict[int, intrangeset]:
    ranges = {1: intrangeset(), -1: intrangeset()}

    for _, feat in iter_features(rec):
        if feat.location is None:
            continue
        for part in feat.location.parts:
            strand = part.strand if part.strand else 1

            start, end = int(part.start), int(part.end)  # type: ignore
            if end > start:
                ranges[strand].add(intrange(start, end))
    return ranges


def intron_walker(
    parts: list[SimpleLocation],
    ranges: dict[int, intrangeset] | None = None,
) -> Iterator[tuple[str, SimpleLocation]]:
    prev = pstrand = None
    for part in parts:
        assert part.strand in {1, -1, 0, None}, (
            part,
            part.strand,
        )  # can be 0(?) or None too
        strand = part.strand
        start, end = int(part.start), int(part.end)  # type: ignore

        s, e = (prev, start) if strand == 1 else (end, prev)
        if s is None or e is None:
            prev = end if strand == 1 else (start if strand == -1 else None)
            continue

        if ranges and prev is not None and strand == pstrand and s < e:
            # meas = (RangeSet(s, e) & ranges[p.strand]).measure()
            #  meas = len(RangeSet(s, e).intersection(ranges[p.strand]))

            # actual intron that doesn't overlap other genes
            if not ranges[strand or 0].overlaps(intrange(s, e)):
                yield "intron", SimpleLocation(s, e, strand)
        # circos fails if start == end even though this is a 1 bp span
        # see /lib/Circos/Karyotype.pm line 134
        assert end >= start
        # at least for tiles
        # circos seems to use Set:IntSpan https://metacpan.org/pod/Set::IntSpan
        # which seems to be [start, end] inclusive! so we reduce the end by 1
        yield "exon", part
        # we're walking backward for strand -1
        pstrand = strand
        prev = end if strand == 1 else (start if strand == -1 else None)


def get_bands(
    rec: SeqRecord,
    span: Callable[[Any], Any] | bool = False,
    number_gene: bool = True,
    introns: bool = True,
    BandCls: type[Band] = Band,
) -> Iterator[Band]:
    ranges = None
    if introns:
        ranges = get_range_map(rec)

    nc: dict[str, int] = Counter()
    found: list[tuple[int, SeqFeature, int]] = []
    for idx, feat in iter_features(rec):
        q = feat.qualifiers
        if "gene" in q:
            gene = q["gene"][0]
            nc[gene] += 1
            c = nc[gene]
        else:
            nc[rec.name] += 1
            c = nc[rec.name]

        found.append((idx, feat, c))

    for idx, feat, count in found:
        q = feat.qualifiers
        typ = feat.type
        if "gene" in q:
            gene = q["gene"][0]
            name = sanitize_name(gene)  # can't have spaces in circos files
            if number_gene and nc[gene] > 1:
                name = f"{name}/{count}"
            typ = "gene"
        else:
            gene = rec.name
            name = sanitize_name(rec.name)
            if nc[rec.name] > 1:
                name = f"{name}/{count}"
                gene = f"{gene}/{count}"
        if span and callable(span):
            parts = span(feat.location)
        else:
            assert feat.location is not None
            parts = feat.location.parts
        for typ, part in intron_walker(parts, ranges):  # type: ignore[arg-type]
            gene_type = typ if typ == "intron" else feat.type

            yield BandCls(
                # name=f"i-{name}",  # we key on this
                name=name,  # now we key on gene_type
                gene=gene,
                start=int(part.start),  # type: ignore
                end=int(part.end),  # type: ignore
                strand=part.strand,  # type: ignore
                idx=idx,  # used to hover identifier in svg
                gene_type=gene_type,
            )


def get_links(
    rec1: SeqRecord,
    rec2: SeqRecord,
    rec1name: str | None = None,
    rec2name: str | None = None,
    min_dist: int = 0,
    get_color: Callable[[Any], str] | None = None,
) -> Iterator[tuple[str, Band, Band, dict[str, Any] | None]]:
    rec1name_str = rec1name or rec1.id or ""
    rec2name_str = rec2name or rec2.id or ""
    b1: dict[str, list[Band]] = defaultdict(list)
    b2: dict[str, list[Band]] = defaultdict(list)
    for b in get_bands(rec1, span=False, number_gene=True):
        b1[b.name].append(b)
    for b in get_bands(rec2, span=False, number_gene=True):
        b2[b.name].append(b)
    genes = set(b1) & set(b2)
    minhit: float | None
    for gene in genes:
        for exon1 in b1[gene]:
            minhit = exon = None
            for exon2 in b2[gene]:
                rng1 = exon1.start, exon1.end
                rng2 = exon2.start, exon2.end
                if rng1 == rng2 and exon1.strand == exon2.strand:
                    continue
                mismatch: float = abs(exon1.start - exon2.start) + abs(
                    exon1.end - exon2.end,
                )
                if minhit is None or mismatch < minhit:
                    minhit = mismatch
                    exon = exon2

            if exon is not None and minhit is not None and minhit >= min_dist:
                # overwrite
                exon1.name = rec1name_str
                exon.name = rec2name_str
                # holdover from circos
                d = dict(color=f"{get_color(exon1)}_a4") if get_color else None
                yield "matched", exon1, exon, d

    # entries that don't match

    def dounmatched(
        unmatched: set[str],
        d: dict[str, list[Band]],
        name: str,
        key: str,
    ) -> Iterator[tuple[str, Band, Band, dict[str, Any] | None]]:
        for gene in unmatched:
            for exon in d[gene]:
                # c = get_color(exon)
                exon.name = name
                # [sic] reverse start and end...
                exon2 = exon.clone(exon.end, exon.start)
                r = ".9r" if exon.strand > 0 else "1.05r"
                yield key, exon, exon2, dict(color="black_a2", bezier_radius=r)

    yield from dounmatched(set(b2) - genes, b2, rec2name_str, "new")
    yield from dounmatched(set(b1) - genes, b1, rec1name_str, "old")
