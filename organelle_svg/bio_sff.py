from __future__ import annotations

import gzip
import logging
import os
import re
from collections import defaultdict
from datetime import date
from pathlib import Path
from typing import IO
from typing import Iterator
from typing import Literal
from typing import NamedTuple
from typing import TypeAlias

from Bio.Seq import Seq
from Bio.SeqFeature import CompoundLocation
from Bio.SeqFeature import SeqFeature
from Bio.SeqFeature import SimpleLocation
from Bio.SeqRecord import SeqRecord

logger = logging.getLogger("chloe")


def wrap(pos: int, genome_length: int) -> int:
    if pos < 0:
        return genome_length + pos
    if pos >= genome_length:
        return pos - genome_length
    return pos


class SFFSeqRecord(SeqRecord):  # pylint: disable=abstract-method
    is_sffseq = True


NAN = re.compile("^nan$", re.I)

F1 = 0x1
F2 = 0x2
F3 = 0x4
ANY = F1 | F2 | F3


def check_nans(i: int) -> str:
    ret = []
    if i & F1:
        ret.append("relative_length")
    if i & F2:
        ret.append("depth")
    if i & F3:
        ret.append("coverage")
    return ",".join(ret)


# * 5th column: feature phase (number of nucleotides to ignore at the start of
#   the range when translating  the feature; only of relevance to CDS features, should
#   be zero for all non-CDS features)
# * 6th column: relative length (length of feature relative to median length of the
#   same feature in the reference genomes, i.e. range can be from 0-2ish, ideally near 1 )
# * 7th column: depth (average relative depth of aligned reference annotations across
#   the feature; should  be in the range 0-1, ideally near 1)
# * 8th column: coverage (proportion of the feature covered by at least
#   one aligned reference annotation; should be in the range 0-1, ideally near 1)
class SFFLocation(SimpleLocation):
    def __init__(
        self,
        start: int,
        end: int,
        strand: int | None,
        key: str,
        phase: int,
        relative_length: float,
        depth: float,
        coverage: float,
        glm: tuple[float, float],
        isnan: int = 0,
        intron: bool = False,
        comment: str | None = None,
    ):
        super().__init__(start, end, strand)
        # extra stuff
        self.key = key  # original key gene/gene_count/type/exon_count
        self.phase = phase
        self.relative_length = relative_length
        self.depth = depth
        self.coverage = coverage
        self.glm = glm
        self.isnan = isnan  # if any of relative_length,depth,coverage are nans
        self.comment = comment
        self.intron = intron


def maybe_gzip(fname: str | Path, mode: Literal["wt", "rt"] = "rt") -> IO[str]:
    if Path(fname).name.endswith(".gz"):
        ret = gzip.open(fname, mode, encoding="utf-8")
    else:
        ret = open(fname, mode, encoding="utf-8")
    return ret  # t.cast(t.TextIO, ret)


def readsff(
    fname: str | Path,
    skip_negative: bool = False,
    strict: bool = False,
    fail_fast: bool = True,
    include_introns: bool = False,
    expand_features: bool = False,
    wrapped: bool = True,
) -> SFFSeqRecord:
    with maybe_gzip(fname) as fp:
        return readsff_fp(
            fp,
            skip_negative=skip_negative,
            strict=strict,
            fail_fast=fail_fast,
            include_introns=include_introns,
            expand_features=expand_features,
            wrapped=wrapped,
        )


LOC: TypeAlias = tuple[str, int, SFFLocation | None]


def readsff_fp(  # noqa: C901
    fp: IO[str],
    skip_negative: bool = False,
    strict: bool = False,
    fail_fast: bool = True,
    include_introns: bool = False,
    expand_features: bool = False,
    wrapped: bool = True,
) -> SFFSeqRecord:
    it = iter_readsff_fp(
        fp,
        skip_negative=skip_negative,
        fail_fast=fail_fast,
        wrapped=wrapped,
    )
    gdict: dict[tuple[str, int], list[LOC]] = defaultdict(list)
    ncid, length, mean_coverage_str, _, _ = next(it)
    assert length > 0, length

    for gene, gene_count, typ, exon_count, location in it:
        gdict[(gene, gene_count)].append((typ, exon_count, location))

    def mknote(locations: list[SFFLocation | None]) -> list[str] | None:
        dd = defaultdict(list)
        for loc in locations:
            if loc is not None:
                if loc.comment:
                    dd[loc.comment].append(loc)
        if not dd:
            return None

        return [
            "; ".join(
                c + ": " + ", ".join(f"({loc.start + 1}..{loc.end})" for loc in locs)
                for c, locs in dd.items()
            ),
        ]

    def get_date() -> str:
        try:
            mt = os.fstat(fp.fileno()).st_mtime
            return date.fromtimestamp(mt).strftime("%d-%b-%Y").upper()
        except Exception:  # pylint: disable=broad-except
            return date.today().strftime("%d-%b-%Y").upper()

    annotations = {
        "molecule_type": "DNA",
        "topology": "circular",
        # "data_file_division": "PLN",
        "date": get_date(),
        "accessions": [ncid],
        "sequence_version": 1,
        "keywords": [],
        "source": "annotated by Chloë",
        "organism": "Unknown",
        "taxonomy": [],
        "references": [],
        "comment": "chloe annotated .sff file",
        "mean_coverage": mean_coverage_str,  # extra
    }

    def create_feature(
        gene: str,
        gene_count: int,
        llocations: list[LOC],
    ) -> SeqFeature | None:
        if not include_introns:
            llocations = [
                t for t in llocations if t[0] != "intron"
            ]  # remove introns...
            if not llocations:
                if strict:
                    raise ValueError(f"no exons for {ncid}:{gene}")
                logger.warning("No exons for gene %s:%s", ncid, gene)
                return None

        types: list[str]
        orders: list[int]
        locations: list[SFFLocation | None]
        types, orders, locations = zip(*llocations)  # type: ignore

        note = mknote(locations)

        xlocation: SFFLocation | CompoundLocation | None
        if len(locations) > 1:
            xtypes = [t for t in types if t != "intron"]
            if len(set(xtypes)) != 1:
                if strict:
                    raise ValueError(f"No type for gene {ncid}:{gene}")
                logger.warning("No type for gene %s:%s got %s", ncid, gene, types)
                xtypes = ["intron"]
                # return None
            feature_type = xtypes[0]
            xlocation = CompoundLocation(list(locations))
        else:
            feature_type, xlocation = (
                types[0],
                locations[0],
            )
        qualifiers: dict[str, list[str] | list[int]] = dict(
            gene=[gene],
            source=["chloe"],
        )
        if note is not None:
            qualifiers["note"] = note
        phases: set[int] = (
            {loc.phase for loc in xlocation.parts} if xlocation else set()
        )
        if len(phases) == 1:
            phase: int = phases.pop()
            if phase != 0:
                if feature_type == "CDS":
                    qualifiers["codon_start"] = [phase + 1]
                else:
                    qualifiers["phase"] = [phase]
        keys: set[str] = {loc.key for loc in xlocation.parts} if xlocation else set()
        if len(keys) == 1:
            key = keys.pop()
            qualifiers["ID"] = [key]
        return SeqFeature(
            xlocation,
            type=feature_type,
            id=f"{gene}/{gene_count}",
            qualifiers=qualifiers,
        )

    rec = SFFSeqRecord(
        Seq(None, length),
        id=ncid,
        name=ncid,
        annotations=annotations,
        description=f"Chloe Annotation for {ncid}",
    )
    for (gene, gene_count), locations in gdict.items():
        if expand_features:
            for typ, order, loc in locations:
                feat = create_feature(gene, gene_count, [(typ, order, loc)])
                if feat is not None:
                    rec.features.append(feat)
        else:
            feat = create_feature(gene, gene_count, locations)
            if feat is not None:
                rec.features.append(feat)

    return rec


class SFFRow(NamedTuple):
    gene: str
    gene_count: int
    type: str
    exon_count: int
    location: SFFLocation | None = None


def iter_readsff_fp(
    fp: IO[str],
    skip_negative: bool = False,
    strict: bool = False,
    fail_fast: bool = True,
    wrapped: bool = True,  # if we split gene when wrapped
) -> Iterator[SFFRow]:
    # pylint: disable = too-many-locals

    header = [s.strip() for s in fp.readline().strip().split("\t")]

    if len(header) == 3:
        ncid, length_, mean_coverage = header
    else:
        mean_coverage = "0.0"
        ncid, length_, *_rest = header

    length = int(length_)
    assert length > 0, "SFF: sequence length is negative!"
    yield SFFRow(ncid, length, mean_coverage, 0)

    for line in fp:
        line = line.strip()
        l2 = [s.strip() for s in line.strip().split("\t")]
        comment = None
        if len(l2) == 8:  # old version
            l2 = [*l2, "0.0", "0.0"]
        elif len(l2) == 9:  # old version
            l2 = [*l2[:-1], "0.0", "0.0", l2[-1]]

        if len(l2) > 10:
            # raise ValueError(f'ERR[{idx+2}]:{len(l2)} {l}')
            l2, comment = l2[:10], " ".join(l2[10:])

        (
            key,
            strand_,
            start_,
            glen_,
            phase_,
            relative_length_,
            depth_,
            coverage_,
            feature_prob_,
            coding_prob_,
        ) = l2
        gene, gene_count_, typ, exon_count_ = key.split("/")
        gene_count, exon_count = map(int, [gene_count_, exon_count_])
        strand = 1 if strand_ == "+" else -1
        start = int(start_)
        phase = int(phase_)
        glen = int(glen_)

        if glen < 1:
            if not strict and comment and "pseudogene" in comment:
                continue
            if skip_negative:
                continue
            msg = f"{ncid}: negative length! {line}"

            if fail_fast:
                raise ValueError(msg)

            logger.error(msg)
            continue
        # this is from when sffs where being generated by
        # old broken code.....
        isnan = 0
        for i, f in enumerate(
            [relative_length_, depth_, coverage_, feature_prob_, coding_prob_],
        ):
            if NAN.match(f):
                isnan |= 1 << i
        if strict and isnan:
            msg = f"{ncid}: NaNs! {check_nans(isnan)} {line}"
            if fail_fast:
                raise ValueError(msg)
            logger.error(msg)

        relative_length, depth, coverage, feature_prob, coding_prob = map(
            float,
            [relative_length_, depth_, coverage_, feature_prob_, coding_prob_],
        )

        lm1 = glen - 1  # length minus 1
        end = start + lm1
        if strand < 0:
            start = length - end + 1

        while start <= 0:
            start += length

        end = start + lm1
        if end > length and wrapped:  # ok wraps
            location = SFFLocation(
                start - 1,
                length,
                strand,
                key,
                phase,
                relative_length,
                depth,
                coverage,
                glm=(feature_prob, coding_prob),
                isnan=isnan,
                intron=typ == "intron",
                comment=comment,
            )
            yield SFFRow(gene, gene_count, typ, exon_count, location)

            lm1 = end - length - 1
            start = 1
            end = start + lm1
        assert not wrapped or (0 < end <= length and 0 < start <= length), (
            f"{key}: 1 <= [{start}..{end}] <= {length} (strand={strand})"
        )
        location = SFFLocation(
            start - 1,
            end,
            strand,
            key,
            phase,
            relative_length,
            depth,
            coverage,
            glm=(feature_prob, coding_prob),
            isnan=isnan,
            comment=comment,
            intron=typ == "intron",
        )
        yield SFFRow(gene, gene_count, typ, exon_count, location)


def readsff_dir(
    directory: str,
    ext: tuple[str, ...] = (".sff", ".sff.gz"),
) -> Iterator[SFFSeqRecord]:
    for f in os.listdir(directory):
        if f.endswith(ext):
            yield readsff(os.path.join(directory, f))
