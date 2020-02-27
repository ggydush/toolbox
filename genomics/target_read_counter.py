from collections import defaultdict
import logging

import click
from quicksect import IntervalTree
import pysam

import bbutils


logger = logging.getLogger("TargetCounter")


def read_bed(bedfile):
    """ Creates generator from bed file or interval_list """
    logger.info("Reading region file...")
    interval_list = bedfile.endswith("interval_list")
    with open(bedfile, "r") as bed:
        for line in bed:
            if line.startswith("@"):
                continue
            line = line.strip()
            chrom, start, stop = line.split()[0:3]
            start, stop = int(start), int(stop)
            if interval_list:
                start -= 1
            yield chrom, start, stop


def create_interval_dict_from_bed(bedfile, padding=0):
    """
    Used for marking on/off target fragments by creating interval
    trees for each chromosome in bedfile. Intervals can be padded
    as well to mimic tools like HSMetrics.
    """
    logger.info("Creating interval tree from region file")
    interval_dict = defaultdict(IntervalTree)
    if bedfile:
        for chrom, start, stop in read_bed(bedfile):
            interval_dict[chrom].add(start - padding, stop + padding)
    return interval_dict


def get_unclipped_read_ends(read):
    """
    Accepts read and returns unclipped ends

    Returns
    -------
    0-based start, 1-based end coordinate of read
    """
    # Check if soft clipping (softclipping = 4)
    l_adj = 0 if read.cigartuples[0][0] != 4 else read.cigartuples[0][1]
    r_adj = 0 if read.cigartuples[-1][0] != 4 else read.cigartuples[-1][1]
    unclipped_start = read.reference_start - l_adj
    unclipped_end = read.reference_end + r_adj
    return unclipped_start, unclipped_end


def get_unclipped_fragment_ends(read1, read2):
    """
    Accepts a read pair and returns the fragment unclipped ends
    Useful when tracking down families after Fgbio

    Returns:
    0-based start, 1-based end coordinate of fragment
     """
    unclipped_ends = (
        get_unclipped_read_ends(read1),
        get_unclipped_read_ends(read2),
    )
    unclipped_ends = [x for sublist in unclipped_ends for x in sublist]
    fragment_start = min(unclipped_ends)
    fragment_end = max(unclipped_ends)
    return fragment_start, fragment_end


def read_pair_generator(bamfile):
    """ Generate fragments from BAM file """
    logger.info("Iterating over BAM file")
    bam = pysam.AlignmentFile(bamfile, "rb")
    for read1, read2 in bbutils.iterate_read_pairs(
        bam, noun="read pairs", log_every=1_000_000
    ):
        yield read1, read2


def get_target(read1, read2, interval_dict, padding):
    """
    Checks if a chrom, start, end overlaps an interval.

    Returns:
    -------
    Region string of interval of overlap or 'Off Bait' if no overlap
    """
    chrom = read1.reference_name
    start, end = get_unclipped_fragment_ends(read1, read2)
    target = "off_bait"
    if chrom in interval_dict:
        res = interval_dict[chrom].search(start, end)
        if res:
            tstart, tend = res[0].start, res[0].end
            target = f"{chrom}:{tstart + padding}-{tend - padding}"

    return target


@click.command(context_settings=dict(help_option_names=["-h", "--help"]))
@click.option(
    "-b",
    "--bamfile",
    type=click.Path(exists=True),
    metavar="FILE",
    help="BAM file",
    required=True,
)
@click.option(
    "-r",
    "--regions",
    type=click.Path(exists=True),
    metavar="FILE",
    help="Regions of interest as BED or interval list",
    required=True,
)
@click.option(
    "-o",
    "--output",
    type=click.Path(writable=True),
    metavar="FILE",
    help="Output file",
    required=True,
)
@click.option(
    "-p",
    "--padding",
    type=click.INT,
    metavar="INT",
    default=250,
    help="Padding for each region",
    show_default=True,
)
def read_counter(bamfile, regions, output, padding):
    """ Count number of read pairs overlapping each region in
    bed file or interval list """
    interval_dict = create_interval_dict_from_bed(regions, padding)
    result = defaultdict(lambda: defaultdict(lambda: 0))
    for read1, read2 in read_pair_generator(bamfile):
        target = get_target(read1, read2, interval_dict, padding)
        result[target]["read_pairs"] += 1
        if ~read1.is_duplicate:
            result[target]["unique_read_pairs"] += 1

    headers = ["target", "read_pairs", "unique_read_pairs"]
    with open(output, "w") as outfile:
        outfile.write("\t".join(headers) + "\n")
        for target, read_counts in result.items():
            data = [
                target,
                read_counts["read_pairs"],
                read_counts["unique_read_pairs"],
            ]
            outfile.write("\t".join(map(str, data)) + "\n")


if __name__ == "__main__":
    read_counter()
