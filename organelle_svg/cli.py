from __future__ import annotations

import gzip
from pathlib import Path
from typing import TYPE_CHECKING

import click

CHLOE = "Chloë"


if TYPE_CHECKING:
    from Bio.SeqRecord import SeqRecord


def read_rec(fname: str | Path, rec_type: str) -> SeqRecord:
    from Bio import SeqIO

    fname = Path(fname)
    if fname.name.endswith(".gz"):
        with gzip.open(fname, "rt", encoding="utf-8") as fp:
            return next(SeqIO.parse(fp, rec_type))

    with open(fname, encoding="utf-8") as fp:
        return next(SeqIO.parse(fp, rec_type))


def readit(gb_or_sff: str) -> SeqRecord:
    from .bio_sff import readsff
    from BCBio import GFF

    if gb_or_sff.endswith((".sff", ".sff.gz")):
        rec = readsff(gb_or_sff, include_introns=True, expand_features=False)
    elif gb_or_sff.endswith((".gff.gz", ".gff3.gz")):
        with gzip.open(gb_or_sff, "rt", encoding="utf-8") as fp:
            rec = next(GFF.parse(fp))
    elif gb_or_sff.endswith((".gff", ".gff3")):
        with open(gb_or_sff, encoding="utf-8") as fp:
            rec = next(GFF.parse(fp))
    else:
        rec = read_rec(gb_or_sff, "genbank")
    if "source" not in rec.annotations:
        rec.annotations["source"] = f"visualized by {CHLOE} ({rec.id})"
    return rec


def all_options(func):
    func = click.option(
        "-b",
        "--bg",
        "bg_color",
        default="white",
        help='background color (use "none" for transparent)',
        show_default=True,
    )(func)
    func = click.option(
        "-o",
        "--output",
        type=click.Path(dir_okay=False, writable=True),
        help="Output SVG file name (default: stdout)",
    )(func)
    func = click.option(
        "-p",
        "--pretty",
        "pretty_print",
        is_flag=True,
        help="pretty print the svg",
    )(func)
    func = click.option(
        "-r",
        "--rotate",
        "rotate_image",
        is_flag=True,
        help="rotate the image so IR is at start of the circle (if present)",
    )(func)

    return func


def out(output: str | None, svg: str) -> None:
    if output is None:
        click.echo(svg)
    else:
        Path(output).write_text(
            svg,
            encoding="utf-8",
        )


def plot_type(choices: list[str]):
    return click.option(
        "-t",
        "--type",
        "plot_type",
        default=choices[0],
        type=click.Choice(choices),
        help="svg type",
        show_default=True,
    )


@click.group(epilog=click.style("Organelle SVG\n", fg="magenta"))
@click.version_option()
def cli() -> None:
    pass


@cli.command()
@all_options
@plot_type(["ogdraw", "gc", "depth"])
@click.argument(
    "gb_or_gff",
    type=click.Path(exists=True, file_okay=True, dir_okay=False),
)
def single(
    gb_or_gff: str,
    plot_type: str,
    output: str | None,
    bg_color: str,
    rotate_image: bool,
    pretty_print: bool,
) -> None:
    """Generate SVG from a genbank or gff file"""
    from .api import OGDraw, DepthDraw, GCOGDraw

    bg = bg_color if bg_color != "none" else False

    rec = readit(gb_or_gff)
    if plot_type == "ogdraw":
        draw = OGDraw(rec, irscan=True, rotate_image=rotate_image, bg=bg)
    elif plot_type == "gc":
        draw = GCOGDraw(rec, irscan=True, rotate_image=rotate_image, bg=bg)
    else:
        draw = DepthDraw(rec, irscan=True, rotate_image=rotate_image, bg=bg)

    svg = draw.to_string(pretty_print=pretty_print)
    out(output, svg)


@cli.command()
@all_options
@plot_type(["normal", "pairs"])
@click.argument(
    "gb_or_gff_inside",
    type=click.Path(exists=True, file_okay=True, dir_okay=False),
)
@click.argument(
    "gb_or_gff_outside",
    type=click.Path(exists=True, file_okay=True, dir_okay=False),
)
def pairs(
    gb_or_gff_inside: str,
    gb_or_gff_outside: str,
    plot_type: str,
    output: str | None,
    bg_color: str,
    rotate_image: bool,
    pretty_print: bool,
) -> None:
    """Generate SVG from pairs of genbank or gff files"""
    from .api import NormalDraw, PairsDraw, BaseDraw

    bg = bg_color if bg_color != "none" else False
    draw: BaseDraw

    rec_in = readit(gb_or_gff_inside)
    rec_out = readit(gb_or_gff_outside)
    if plot_type == "normal":
        draw = NormalDraw(
            rec_in,
            rec_out,
            irscan=True,
            rotate_image=rotate_image,
            bg=bg,
        )
    else:
        draw = PairsDraw(
            rec_in,
            rec_out,
            irscan=True,
            rotate_image=rotate_image,
            bg=bg,
        )

    svg = draw.to_string(pretty_print=pretty_print)
    out(output, svg)


@cli.command()
@all_options
@click.option(
    "--name",
    help="name for svg file",
)
@click.argument(
    "gb_or_gff",
    type=click.Path(exists=True, file_okay=True, dir_okay=False),
    nargs=-1,
    required=True,
)
def stacked(
    gb_or_gff: list[str],
    name: str | None,
    output: str | None,
    bg_color: str,
    rotate_image: bool,
    pretty_print: bool,
) -> None:
    """Generate SVG from a list of genbank or gff files"""
    from .api import StackedDraw

    bg = bg_color if bg_color != "none" else False
    recs = [readit(f) for f in gb_or_gff]

    draw = StackedDraw(
        name or ", ".join(Path(f).stem for f in gb_or_gff),
        recs,
        irscan=True,
        rotate_image=rotate_image,
        bg=bg,
    )

    svg = draw.to_string(pretty_print=pretty_print)
    out(output, svg)
