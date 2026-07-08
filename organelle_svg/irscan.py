from __future__ import annotations

from typing import NamedTuple
from typing import TYPE_CHECKING


if TYPE_CHECKING:
    from typing import Any
    from Bio.Seq import MutableSeq
    from Bio.Seq import Seq


class IRScanResult(NamedTuple):
    IRA_start: int
    IRA_end: int
    IRB_start: int
    IRB_end: int


WORD_SIZE = 2000
ROUGH = 500


# seems to just as fast as the OGDraw.irscan (possibly when
# external program invocation overhead is taken into account)


class Scanner:
    def __init__(
        self,
        seq: Seq | MutableSeq,
        wordsize: int = WORD_SIZE,
        precision: int = 1,
        rough: int = ROUGH,
    ):
        # don't need unicode use byte comparisons....
        self.seq = str(seq).encode("ascii")
        # memoryview so rev[i:j] just creates a view and does not copy `wordsize` bytes
        self.rev = memoryview(str(seq.reverse_complement()).encode("ascii"))
        self.precision = precision
        self.wordsize = wordsize
        self.rough = rough

    def screen_from_to_pos(
        self,
        pos: int,
        endpos: int,
        step: int,
    ) -> tuple[int, int] | tuple[None, None]:
        seq = self.seq
        rev = self.rev
        wordsize = self.wordsize

        nseq = len(seq)
        start, end = None, None
        while True:
            if pos + wordsize > endpos:
                return None, None

            # test_word = seq[pos : pos + wordsize].reverse_complement()
            test_word = rev[nseq - pos - wordsize : nseq - pos]
            end = seq.find(test_word)
            if end != -1:
                end += wordsize
                start = pos + 1
                break

            pos += step
        return start, end

    def rough_scan(self) -> tuple[int, int, int, int] | None:
        seq = self.seq
        rev = self.rev
        wordsize = self.wordsize
        rough = self.rough

        position = 0
        nseq = len(seq)

        while True:
            if position + wordsize > nseq:
                return None

            # test_word = sequence[position : position + wordsize].reverse_complement()
            test_word = rev[nseq - position - wordsize : nseq - position]
            irb_end = seq.find(test_word)
            if irb_end != -1:
                irb_end += wordsize
                ira_start = position
                break
            position += rough
        else:
            return None

        ssc_center = (ira_start + irb_end) // 2

        irb_start, ira_end = self.screen_from_to_pos(ssc_center, nseq, rough)
        if irb_start is None or ira_end is None:
            return None
        return ira_start, ira_end, irb_start, irb_end

    def irscan(self) -> IRScanResult | None:
        pre_borders = self.rough_scan()
        if pre_borders is None:
            return None
        ira_start, _, irb_start, _ = pre_borders

        IRA_start, IRB_end = self.screen_from_to_pos(
            max(0, ira_start - self.rough),
            len(self.seq),
            self.precision,
        )
        if IRA_start is None or IRB_end is None:
            return None

        IRB_start, IRA_end = self.screen_from_to_pos(
            max(0, irb_start - self.rough),
            len(self.seq),
            self.precision,
        )
        if IRB_start is None or IRA_end is None:
            return None

        return IRScanResult(IRA_start, IRA_end, IRB_start, IRB_end)


class IRScan:
    def doirscan(
        self,
        name: str,
        seq: Seq | MutableSeq | None,
        **_kwargs: Any,
    ) -> IRScanResult | None:
        if seq is None:
            return None
        return Scanner(seq).irscan()
