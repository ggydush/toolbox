import click
from collections import defaultdict
import pandas as pd
import pysam
from quicksect import IntervalTree

import bbutils


def read_bed(bedfile):
    """ Creates generator from bed file """
    with open(bedfile, "r") as bed:
        for line in bed:
            line = line.strip()
            chrom, start, stop = line.split()[0:3]
            yield chrom, int(start), int(stop)


def create_interval_dict_from_bed(bedfile, padding=0):
    """
    Used for marking on/off target fragments by creating interval
    trees for each chromosome in bedfile. Intervals can be padded
    as well to mimic tools like HSMetrics.
    """
    interval_dict = defaultdict(IntervalTree)
    if bedfile:
        for chrom, start, stop in read_bed(bedfile):
            interval_dict[chrom].add(start - padding, stop + padding)
    return interval_dict


def get_target(read1, read2, interval_dict, full_overlap):
    """
    Checks if a read pair overlaps an interval.

    Returns:
    -------
    Region string of interval of overlap or NA if no overlap
    """
    target = "NA"
    if len(interval_dict) == 0:
        return target
    chrom = read1.reference_name
    start, end = get_unclipped_fragment_ends(read1, read2)

    if chrom in interval_dict:
        res = interval_dict[chrom].search(start, end)
        if res:
            tstart, tend = res[0].start, res[0].end
            target = f"{chrom}:{tstart}-{tend}"
            if full_overlap:
                if start >= tstart or end <= tend:
                    target = "NA"
    return target


def _get_unclipped_read_ends(read):
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

    Returns:
    0-based start, 1-based end coordinate of fragment
     """
    unclipped_ends = (
        _get_unclipped_read_ends(read1),
        _get_unclipped_read_ends(read2),
    )
    unclipped_ends = [x for sublist in unclipped_ends for x in sublist]
    fragment_start = min(unclipped_ends)
    fragment_end = max(unclipped_ends)
    return fragment_start, fragment_end


def get_family_id(read):
    """ Get family ID from read """
    return read.get_tag("MI")


def generate_umi_histogram_dict(bamfile, bedfile=None, full_overlap=False):
    """
    Generate a UMI histogram dict (same as output from GroupReadsByUmi)
    :param bam - pysam.AlignmentFile object
    :returns summary_dict with key=family_size, value=count
    """
    bam = pysam.AlignmentFile(bamfile, "rb")
    interval_dict = create_interval_dict_from_bed(bedfile, padding=0)
    summary_dict = defaultdict(lambda: defaultdict(lambda: 0))
    print("\nBegin iteration through BAM")
    for idx, (read1, read2) in enumerate(
        bbutils.iterate_read_pairs(bam, until_eof=True, log_every=1000000)
    ):
        if idx == 0:
            family_size = 0
            previous_family_id = get_family_id(read1)
            previous_target = get_target(
                read1, read2, interval_dict, full_overlap
            )

        current_family_id = get_family_id(read1)

        if current_family_id == previous_family_id:
            family_size += 1
        else:
            summary_dict[previous_target][family_size] += 1
            family_size = 1
            previous_family_id = current_family_id
            previous_target = get_target(
                read1, read2, interval_dict, full_overlap
            )

    # Account for very last family
    summary_dict[previous_target][family_size] += 1

    # Clean it up
    bam.close()
    return summary_dict


def calculate_statistics(umi_histogram):
    """ Add info to dataframe so that it matches umi_histogram.txt """
    total_count = umi_histogram["count"].sum().astype(float)
    umi_histogram["fraction"] = umi_histogram["count"] / total_count
    fraction_gt_or_eq = 1 - umi_histogram["fraction"].cumsum()[:-1]
    umi_histogram["fraction_gt_or_eq_family_size"] = (
        pd.Series([1]).append(fraction_gt_or_eq, ignore_index=True).values
    )
    target = umi_histogram["target"]
    umi_histogram.drop("target", axis=1, inplace=True)
    umi_histogram["target"] = target
    return umi_histogram


def create_dataframe_from_summary_dict(summary_dict):
    """ Create dataframe similar to umi_histogram.txt from GroupReadsByUmi """
    print("\nCreating summary file")
    umi_histogram = (
        pd.DataFrame.from_dict(summary_dict)
        .reset_index()
        .melt(id_vars=["index"])
    )
    umi_histogram.columns = ["family_size", "target", "count"]
    umi_histogram.dropna(inplace=True)
    res = []
    for gid, gdf in umi_histogram.groupby("target"):
        res.append(calculate_statistics(gdf.copy()))
    return pd.concat(res)


def create_output_file(ctx, param, value):
    """
    Used to validate output argument. If no output is given,
    use given filename
    """
    if value is None:
        bamfile = ctx.params["bamfile"].split("/")[-1]
        value = bamfile.split(".Group")[0] + "_umi_histogram.txt"
    return value


@click.command(context_settings=dict(help_option_names=["-h", "--help"]))
@click.option(
    "-i",
    "--bamfile",
    type=click.Path(exists=True),
    required=True,
    help="BAM file from GroupReadsByUmi",
)
@click.option(
    "-b",
    "--bedfile",
    type=click.Path(exists=True),
    help="BED file specifying targets",
)
@click.option(
    "-o",
    "--output",
    type=click.Path(writable=True),
    callback=create_output_file,
    help="Output file [<sample>_umi_histogram.txt]",
)
@click.option(
    "-f",
    "--full_overlap",
    help="Require full overlap of target interval",
    is_flag=True,
)
def main(bamfile, bedfile, output, full_overlap):
    """
    \b
    Script to create UMI count histogram similar to that in Fgbio
    GroupReadsByUmi. The difference here, is that this can run
    on a GroupReadsByUmi bam (useful when merge/splitting) and
    you can also specify targets via a BED file.

    If target is not in specified targets (or no BED file is given)
    resulting dataframe column target will equal NA
    """
    print(f"Running {__file__} with arguments:")
    args = [f"{k}: {v}" for k, v in click.get_current_context().params.items()]
    print("\n".join(args))
    umi_dict = generate_umi_histogram_dict(bamfile, bedfile, full_overlap)
    umi_hist_df = create_dataframe_from_summary_dict(umi_dict)
    umi_hist_df.to_csv(output, sep="\t", index=None)


if __name__ == "__main__":
    main()
