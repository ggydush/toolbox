import bbutils
import click
import pandas as pd
import pybedtools
import pysam

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


def get_reference_lengths(bam, chromosomes=None):
    """
    By default, all chromosome lengths will be returned.
    If chromosomes are specified (as iterator) only specified
    chromosomes will be returned.
    """
    reference_lengths = zip(bam.references, bam.lengths)
    if "all" in chromosomes:
        return list(reference_lengths)
    elif chromosomes:
        reference_lengths = (
            (x, y) for x, y in reference_lengths if x in chromosomes
        )
    return list(reference_lengths)


def generate_windows(chrom, length, window):
    """ Generates windows based on chromosome length """
    end = 0  # when window > length
    interval_starts = range(0, length - window, window)
    interval_ends = range(window, length, window)
    for start, end in zip(interval_starts, interval_ends):
        yield chrom, start, end
    yield chrom, end, length


def get_read_counts(bam, chromosomes, window, mapq):
    """
    Generate read counts for specified chromosomes and
    window given a minimum mapping quality

    Each window within each chromosome will be given a read count
    """
    reference_lengths = get_reference_lengths(bam, chromosomes)
    print(reference_lengths)
    previous_read_names = set()
    for chromosome, length in reference_lengths:
        for chrom, start, end in generate_windows(chromosome, length, window):
            current_read_names = set()
            read_count = 0
            for read in bbutils.iterate_reads(
                bam, chrom, start, end, log_every=window
            ):
                read_name = f"{read.query_name}_{int(read.is_read2) + 1}"
                if (
                    read.mapping_quality >= mapq
                    and read_name not in previous_read_names
                ):
                    current_read_names.add(read_name)
                    read_count += 1
            previous_read_names = current_read_names
            yield chrom, start, end, read_count


def intersect_dataframes(counts_df, target_df):
    """
    Intersect counts dataframe with target dataframe.
    Outputs bed-like dataframe with columns:
    chrom, start, end, read_counts, num_targets
    """
    read_bed = pybedtools.BedTool.from_dataframe(counts_df)
    target_bed = pybedtools.BedTool.from_dataframe(target_df)
    output_df = (
        read_bed.intersect(target_bed, wa=True, c=True).to_dataframe().iloc[:]
    )
    output_df.rename(
        columns={"name": "read_counts", "score": "num_targets"}, inplace=True
    )
    return output_df


def create_plot(
    counts_df, output_file, on_target_color="red", off_target_color="black"
):
    """ Create genome wide plot coloring on and off target when available """
    counts_df["color"] = (counts_df.num_targets > 0).replace(
        {True: on_target_color, False: off_target_color}
    )
    on_target = counts_df[counts_df.num_targets > 0]
    off_target = counts_df[counts_df.num_targets == 0]

    chrom_min = counts_df.groupby("chrom")["start"].idxmin()
    chrom_max = counts_df.groupby("chrom")["start"].idxmax()
    label_positions = chrom_min + (chrom_max - chrom_min) / 2
    label_positions = sorted(
        zip(label_positions.index, label_positions), key=lambda x: x[1]
    )

    fig, ax = plt.subplots(figsize=(11, 3))
    ax.bar(
        off_target.index,
        off_target.read_counts,
        color=off_target_color,
        edgecolor=off_target_color,
        label="Off Target",
    )
    if not on_target.empty:
        ax.bar(
            on_target.index,
            on_target.read_counts,
            color=on_target_color,
            edgecolor=on_target_color,
            label="On Target",
        )
    ax.set_xticks([x[1] for x in label_positions])
    ax.set_xticklabels([x[0] for x in label_positions])
    ax.set_xlim(0, max(chrom_max) + 50)
    ax.set(
        xlabel="Chromosome",
        ylabel="Read Counts",
        title=output_file.replace(".png", ""),
    )
    ax.legend()

    for xc in chrom_max:
        ax.axvline(x=xc, color="k", linestyle="--", linewidth=0.5, alpha=0.2)
    plt.tight_layout()
    plt.savefig(output_file)


@click.command(
    context_settings=dict(
        help_option_names=["-h", "--help"], max_content_width=100
    )
)
@click.option(
    "-b",
    "--bamfile",
    help="BAM file",
    type=click.Path(exists=True),
    required=True,
)
@click.option(
    "-f",
    "--bedfile",
    help="BED file (used for coloring on-target) [optional]",
    type=click.Path(exists=True),
)
@click.option(
    "-o",
    "--prefix",
    help="Output prefix [inferred from BAM filename]",
    default=None,
)
@click.option(
    "-c",
    "--chromosomes",
    help="Quoted space-separated list of chromsomes [1-22, X, Y]",
    default=None,
)
@click.option(
    "-w",
    "--window",
    help="Window size for counting [1M]",
    type=click.INT,
    metavar="INT",
    default=1000000,
)
@click.option(
    "-m",
    "--mapq",
    help="Mapping quality for read counting [0]",
    type=click.INT,
    metavar="INT",
    default=0,
)
@click.option(
    "-s",
    "--skip",
    help="Skip creating genome-wide plot",
    is_flag=True,
    default=False,
)
def read_counter(bamfile, bedfile, prefix, chromosomes, window, mapq, skip):
    """
    \b
    Script to count reads, create SEG file and genome-wide plot.

    Genome-wide plots require window size to be >= 500kb.
    For reference genomes other than human, list chromosomes or use 'all'
    """
    # Clean up arguments
    if chromosomes:
        chromosomes = chromosomes.split(" ")
    else:
        chromosomes = list(map(str, range(1, 23))) + ["X", "Y"]
    if not prefix:
        prefix = bamfile.split("/")[-1].replace(".bam", "")
    if not skip and window < 500_000:
        print(
            "Plotting only works with large window sizes! Increase window or skip plotting."
        )
        raise click.Abort()

    # Count Reads
    print("Iterating over BAM file")
    bam = pysam.AlignmentFile(bamfile, "rb")
    read_counts = get_read_counts(
        bam, chromosomes=chromosomes, window=window, mapq=mapq
    )
    counts_df = pd.DataFrame(
        [x for x in read_counts],
        columns=["chrom", "start", "end", "read_count"],
    )
    total_counted_reads = counts_df["read_count"].sum()
    print(counts_df)

    # Add on/off target note
    if bedfile:
        interval_list = False
        with open(bedfile, "r") as infile:
            print(next(infile))
            if next(infile).startswith("@"):
                interval_list = True
        if interval_list:
            target_df = pd.read_table(
                bedfile,
                comment="@",
                usecols=[0, 1, 2],
                names=["chrom", "start", "end"],
            )
            target_df["start"] = target_df["start"] - 1
        else:
            target_df = pd.read_table(
                bedfile, usecols=[0, 1, 2], names=["chrom", "start", "end"]
            )
        counts_df = intersect_dataframes(counts_df, target_df)
    else:
        target_df = pd.DataFrame(None)
    counts_df.to_csv(prefix + "_read_counts.seg", sep="\t", index_label="name")
    print(
        f"Counted {total_counted_reads} reads in BAM file for specified chromosomes"
    )

    # Create Genome Wide Plot
    if not skip:
        print("Creating genome-wide plot...")
        counts_df["color"] = (counts_df["num_targets"] > 0).replace(
            {True: "red", False: "black"}
        )
        create_plot(counts_df, prefix + "_genome_wide_plot.png")


if __name__ == "__main__":
    read_counter()
