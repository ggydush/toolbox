""" Used for HiSeqX/HiSeq4000 comparisons """


from collections import defaultdict
import os

from quicksect import IntervalTree
import bbutils
import pysam
import click


class Fragment:
    def __init__(self, read1, read2):
        self.read1 = read1
        self.read2 = read2
        self.family_id = self.get_read_tag("MI")
        self.family_size = self.get_read_tag("cD")
        self.umi = self._get_corrected_umi()
        self.chrom = read1.reference_name
        self.start, self.end = self._get_unclipped_fragment_ends()
        self.insert_size = self.end - self.start
        self.fragment_id = f"{self.chrom}:{self.start}-{self.end}"
        self.mapq1 = self.read1.mapping_quality
        self.mapq2 = self.read2.mapping_quality
        self.target = None  # assigned when check_if_on_target is called

        self.set_mismatches()

    def get_read_tag(self, tag, default="NA"):
        """ Attempts to get read tag from BAM. If it can't find it, the
        tag will be filled in as NA """
        try:
            return self.read1.get_tag(tag)
        except ValueError:
            return default

    @staticmethod
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

    def _get_unclipped_fragment_ends(self):
        """
        Accepts a read pair and returns the fragment unclipped ends
        Useful when tracking down families after Fgbio

        Returns:
        0-based start, 1-based end coordinate of fragment
         """
        read1, read2 = self.read1, self.read2
        unclipped_ends = (
            self._get_unclipped_read_ends(read1),
            self._get_unclipped_read_ends(read2),
        )
        unclipped_ends = [x for sublist in unclipped_ends for x in sublist]
        fragment_start = min(unclipped_ends)
        fragment_end = max(unclipped_ends)
        return fragment_start, fragment_end

    def _get_corrected_umi(self):
        """
        Always get UMI from read1 positive strand orientation

        Note: this only really works for SSC files.
        Still need to figure out a fix for DSC
        """
        read1 = self.read1
        umi = self.get_read_tag("RX")
        if read1.is_reverse and umi != "NA":
            umi = "-".join(umi.split("-")[::-1])
        return umi

    def check_if_on_target(self, interval_dict, padding):
        """
        Checks if a chrom, start, end overlaps an interval.

        Returns:
        -------
        Region string of interval of overlap or NA if no overlap
        """
        chrom, start, end = self.chrom, self.start, self.end
        if chrom in interval_dict:
            res = interval_dict[chrom].search(start, end)
            if res:
                tstart, tend = res[0].start, res[0].end
                self.target = f"{chrom}:{tstart + padding}-{tend - padding}"
                return
        self.target = "NA"

    @staticmethod
    def get_mismatched_counts(read):
        """
        Get all mismatched bases in read.

        Parameters
        ----------
        read : pysam.AlignedSegment
            Pysam Read object

        Returns
        -------
        tuple
            tuple(num_Ns, num_MMs)
        """
        aligned_read_pairs = read.get_aligned_pairs(with_seq=True)
        arr = [
            x for x in aligned_read_pairs if x[-1] and x[0]
        ]  # Remove indels/soft
        mismatches = [x for x in arr if not x[-1].isupper()]  # remove matches
        count_ns = 0
        count_mm = 0
        read_sequence = read.query_sequence
        for idx, start, refbase in mismatches:
            alt_base = read_sequence[idx]
            if alt_base == "N":
                count_ns += 1
            else:
                count_mm += 1
        return count_ns, count_mm

    def set_mismatches(self):
        """
        Add mismatch attributes for number Ns and number of mismatches
        in read
        """
        r1_ns, r1_mm = self.get_mismatched_counts(self.read1)
        r2_ns, r2_mm = self.get_mismatched_counts(self.read2)
        self.r1_ns = r1_ns
        self.r1_mm = r1_mm
        self.r2_ns = r2_ns
        self.r2_mm = r2_mm


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


def fragment_generator(bamfile):
    """ Generate fragments from BAM file """
    bam = pysam.AlignmentFile(bamfile, "rb")
    for read1, read2 in bbutils.iterate_read_pairs(
        bam, noun="consensus families", log_every=10000
    ):
        yield Fragment(read1, read2)


def write_fragments(bamfile, bedfile, output_file, padding):
    """ Summarize fragments found in BAM and write to file """
    interval_dict = create_interval_dict_from_bed(bedfile, padding=padding)
    headers = [
        "family_id",
        "fragment_id",
        "insert_size",
        "family_size",
        "umi",
        "mapq1",
        "mapq2",
        "target",
    ]
    with open(output_file, "w") as outfile:
        outfile.write("\t".join(headers) + "\n")
        fragments = fragment_generator(bamfile)
        for fragment in fragments:
            fragment.check_if_on_target(interval_dict, padding)
            fragment_info = [fragment.__dict__.get(x, None) for x in headers]
            outfile.write("\t".join(map(str, fragment_info)) + "\n")


@click.command(
    help=(
        "Create dataframe of all DSC or SSC families for plotting."
        "Accepts BAM file or BAM file list and a BED file "
        "for annotating on/off target families"
    )
)
@click.option("-b", "--bam", metavar="FILE", help="BAM file")
@click.option(
    "-f",
    "--filelist",
    metavar="FILE",
    help="File list containing BAMs (one per line)",
)
@click.option(
    "-i",
    "--bedfile",
    metavar="FILE",
    help="BED file to annotate on/off targets [optional]",
)
@click.option(
    "-p",
    "--padding",
    metavar="INT",
    help="Padding for determining on target",
    default=250,
    show_default=True,
)
@click.option(
    "-d",
    "--output_dir",
    metavar="DIR",
    help="Directory for output files",
    default="./",
    show_default=True,
)
@click.option(
    "-o",
    "--output_file",
    metavar="FILE",
    help="Output file (if running on single BAM)",
)
def create_summary(bam, filelist, bedfile, padding, output_dir, output_file):
    suffix = "_consensus_summary.txt"
    if bam:
        filelist = [bam]
        output_files = [
            os.path.join(output_dir, output_file)
            if output_file
            else output_dir + os.path.basename(bam).replace(".bam", suffix)
        ]
    else:
        filelist = [x.strip() for x in open(filelist).readlines()]
        output_files = [
            output_dir + os.path.basename(x).replace(".bam", suffix)
            for x in filelist
        ]
    for bamfile, output_file in zip(filelist, output_files):
        print(f"Working on {os.path.basename(bamfile)}")
        write_fragments(bamfile, bedfile, output_file, padding)


if __name__ == "__main__":
    create_summary()
