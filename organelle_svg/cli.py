from __future__ import annotations

import gzip
from functools import wraps
from pathlib import Path
from typing import TYPE_CHECKING, Any

import click

if TYPE_CHECKING:
    from collections.abc import Callable

CHLOE = "Chloë"


if TYPE_CHECKING:
    from Bio.SeqRecord import SeqRecord

    from .api import DrawStyle


def read_rec(fname: str | Path, rec_type: str) -> SeqRecord:
    from Bio import SeqIO

    fname = Path(fname)
    if fname.name.endswith(".gz"):
        with gzip.open(fname, "rt", encoding="utf-8") as fp:
            return next(SeqIO.parse(fp, rec_type))  # type: ignore

    with open(fname, encoding="utf-8") as fp:
        return next(SeqIO.parse(fp, rec_type))  # type: ignore


def readit(gb_or_sff: str) -> SeqRecord:
    from BCBio import GFF

    from .bio_sff import readsff

    rec: SeqRecord

    if gb_or_sff.endswith((".sff", ".sff.gz")):
        rec = readsff(gb_or_sff, include_introns=True, expand_features=False)
    elif gb_or_sff.endswith((".gff.gz", ".gff3.gz")):
        with gzip.open(gb_or_sff, "rt", encoding="utf-8") as fp:
            rec = next(GFF.parse(fp))  # type: ignore
    elif gb_or_sff.endswith((".gff", ".gff3")):
        with open(gb_or_sff, encoding="utf-8") as fp:
            rec = next(GFF.parse(fp))  # type: ignore
    else:
        rec = read_rec(gb_or_sff, "genbank")
    if "source" not in rec.annotations:
        rec.annotations["source"] = f"visualized by {CHLOE} ({rec.id})"
    return rec


def style_options(func: Callable[..., None]) -> Callable[..., None]:
    from .api import DrawStyle

    @wraps(func)
    @click.option(
        "-b",
        "--bg",
        default="white",
        help='background color (use "none" for no background)',
        show_default=True,
    )
    @click.option(
        "-c",
        "--stroke-circles",
        default="grey",
        help="stroke color for bands and IRA/IRB circles",
        show_default=True,
    )
    @click.option(
        "-s",
        "--stroke",
        default="grey",
        help="stroke color for radial lines",
        show_default=True,
    )
    @click.option(
        "-t",
        "--text-color",
        "text_color",
        default="black",
        help="text color for gene names and labels",
        show_default=True,
    )
    def style_options_inner(
        bg: str,
        text_color: str,
        stroke_circles: str,
        stroke: str,
        **kwargs: Any,
    ) -> None:
        style = DrawStyle(
            bg=bg,
            text_color=text_color,
            stroke_circles=stroke_circles,
            stroke=stroke,
        )
        kwargs["style"] = style

        return func(**kwargs)

    return style_options_inner


def all_options(func: Callable[..., None]) -> Callable[..., None]:
    func = click.option(
        "-o",
        "--output",
        type=click.Path(dir_okay=False, writable=True),
        help="Output SVG file name (default: stdout)",
    )(func)
    func = click.option(
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
    func = click.option(
        "--no-legend",
        "no_legend",
        is_flag=True,
        default=False,
        help="do not add legend",
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


def plot_type(
    choices: list[str],
) -> Callable[[Callable[..., None]], Callable[..., None]]:
    return click.option(
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
@plot_type(["ogdraw", "depth"])
@click.argument(
    "gb_or_gff",
    type=click.Path(exists=True, file_okay=True, dir_okay=False),
)
@style_options
def single(
    gb_or_gff: str,
    plot_type: str,
    output: str | None,
    style: DrawStyle,
    rotate_image: bool,
    pretty_print: bool,
    no_legend: bool,
) -> None:
    """Generate SVG from a genbank or gff file"""
    from .api import BaseDraw, DepthDraw, OGDraw

    if no_legend:
        BaseDraw.add_legend = False

    rec = readit(gb_or_gff)
    if plot_type == "ogdraw":
        draw = OGDraw(rec, irscan=True, rotate_image=rotate_image, style=style)
    else:
        draw = DepthDraw(rec, irscan=True, rotate_image=rotate_image, style=style)

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
@style_options
def pairs(
    gb_or_gff_inside: str,
    gb_or_gff_outside: str,
    plot_type: str,
    output: str | None,
    style: DrawStyle,
    rotate_image: bool,
    pretty_print: bool,
    no_legend: bool,
) -> None:
    """Generate SVG from pairs of genbank or gff files"""
    from .api import BaseDraw, NormalDraw, PairsDraw

    if no_legend:
        BaseDraw.add_legend = False
    draw: BaseDraw

    rec_in = readit(gb_or_gff_inside)
    rec_out = readit(gb_or_gff_outside)
    if plot_type == "normal":
        draw = NormalDraw(
            rec_in,
            rec_out,
            irscan=True,
            rotate_image=rotate_image,
            style=style,
        )
    else:
        draw = PairsDraw(
            rec_in,
            rec_out,
            irscan=True,
            rotate_image=rotate_image,
            style=style,
        )

    svg = draw.to_string(pretty_print=pretty_print)
    out(output, svg)


@cli.command()
@all_options
@click.option(
    "--name",
    help="name for svg file",
)
@style_options
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
    style: DrawStyle,
    rotate_image: bool,
    pretty_print: bool,
    no_legend: bool,
) -> None:
    """Generate SVG from a list of genbank or gff files"""
    from .api import BaseDraw, StackedDraw

    if no_legend:
        BaseDraw.add_legend = False

    recs = [readit(f) for f in gb_or_gff]

    draw = StackedDraw(
        name or ", ".join(Path(f).stem for f in gb_or_gff),
        recs,
        irscan=True,
        rotate_image=rotate_image,
        style=style,
    )

    svg = draw.to_string(pretty_print=pretty_print)
    out(output, svg)
