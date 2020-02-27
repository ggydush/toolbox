#!/usr/bin/env python

import subprocess
import os
import sys
import click


def print_help(msg):
    click.echo(f"ERROR: {msg}\n")
    ctx = click.get_current_context()
    click.echo(ctx.get_help())
    ctx.exit()


@click.command(
    context_settings=dict(
        help_option_names=["-h", "--help"],
        allow_extra_args=True,
        ignore_unknown_options=True,
    )
)
@click.argument("file", default="-")
@click.argument("args", nargs=-1, required=False)
@click.option(
    "-s",
    "--spacing",
    type=int,
    metavar="INT",
    help="Spacing between columns (tab width)",
    default=15,
    show_default=True,
)
def view_bam(file, spacing, args):
    """
    \b
    Command for viewing BAMs as columns that are properly spaced.
    Supports piping (reading from stdin) and all samtools view
    arguments apart from "-s" and "-h"
    """
    args = " ".join(args)
    if file != "-":
        if not os.path.exists(file):
            print_help(f"{file} does not exist")
    else:
        if sys.stdin.isatty():
            print_help("No file or pipe used")

    view_cmd = f"samtools view {args} {file} | expand -t {spacing} | less -S"
    subprocess.call([view_cmd], shell=True)


if __name__ == "__main__":
    view_bam()
