from collections import namedtuple
import pandas as pd

import bbutils
import click
import pysam


class TargetMolecule:
    def __init__(
        self, pileupread, target_pos, target_ref_base, target_alt_base
    ):
        self.pileupread = pileupread
        self.target_pos = target_pos
        self.target_ref_base = target_ref_base
        self.target_alt_base = target_alt_base
        self.read = pileupread.alignment
        self.chrom = self.read.reference_name
        self.read_start, self.read_end = self._get_unclipped_ends(
            self.read.reference_start, self.read.cigarstring
        )
        self.read_sequence = self.read.query_sequence
        self.base = self.read.query_sequence[pileupread.query_position]
        self.start, self.end = self._get_unclipped_fragment_ends()
        self.mapping_quality = self.read.mapping_quality
        self.base_quality = self.read.query_qualities[pileupread.query_position]
        self.family_id = self.get_read_tag("MI")
        self.family_size = self.get_read_tag("cD")
        self.umi = self.get_read_tag(
            "RX"
        )  # note this should be fixed when Fgbio folks respond
        self.overlap = False
        self.discordant = False
        self.mismatch_count = self.get_read_tag(
            "NM"
        ) - pileupread.alignment.query_alignment_sequence.count(
            "N"
        )  # added by Chris
        self.N_count = pileupread.alignment.query_alignment_sequence.count(
            "N"
        )  # added by Chris

    def __str__(self):
        args = [
            f'{k}="{v}"'
            for k, v in self.__dict__.items()
            if k not in ["pileupread", "read"]
        ]
        args = ", ".join(args)
        return f"{self.__class__.__name__}({args})"

    def get_read_tag(self, tag, default="NA"):
        try:
            return self.read.get_tag(tag)
        except KeyError:
            return default

    def get_mapping_qualities(self):
        """ Returns tuple: (read1.mapping_quality, read2.mapping_quality) """
        res = (self.read.mapping_quality, self.get_read_tag("MQ"))
        if not self.read.is_read1:
            res = res[::-1]
        return res

    @staticmethod
    def _get_cigar_tuples(cigar):
        """ Convert cigar string to cigar tuples (num bases, match type) """
        total = ""
        res = []
        for x in cigar:
            if x.isdigit():
                total += x
            else:
                res.append((int(total), x))
                total = ""
        return res

    def _get_unclipped_ends(self, start, cigar):
        """ Pysam reports reference positions of reads in terms of
        alignment positions whereas Fgbio reports read ends regardless
        of alignment (aka softclippings included). This function takes
        in a start position and cigar and will return the unclipped
        start and unclipped end.
        """
        end = start
        cigar_tuples = self._get_cigar_tuples(cigar)
        if cigar_tuples[0][1] == "S":
            start -= cigar_tuples[0][0]
            cigar_tuples.pop(0)
        for l, kind in cigar_tuples:
            if kind in ["M", "D", "S"]:
                end += l
        return start, end

    def _get_unclipped_fragment_ends(self):
        """ Same as _get_unclipped_ends but returns fragment level start
        and end rather than read level """
        a_start = self.read.reference_start
        a_cigar = self.read.cigarstring
        b_start = self.read.next_reference_start
        b_cigar = self.get_read_tag("MC")
        template_length = abs(self.read.template_length)

        if "S" not in a_cigar and "S" not in b_cigar:
            start = min([a_start, b_start])
            start, end = start, start + template_length
        else:
            unclipped_ends = (
                self._get_unclipped_ends(a_start, a_cigar),
                self._get_unclipped_ends(b_start, b_cigar),
            )
            unclipped_ends = [x for sublist in unclipped_ends for x in sublist]
            start, end = min(unclipped_ends), max(unclipped_ends)
        return start, end

    @property
    def aD_bD(self):
        """ aD and bD tags (sorted so max is first) """
        ad, bd = self.get_read_tag("aD"), self.get_read_tag("bD")
        ad, bd = sorted([int(ad), int(bd)], reverse=True)
        return f"{ad}/{bd}"

    @property
    def insert_size(self):
        """ Return insert size with unclipped ends taken into account """
        return self.end - self.start

    @property
    def read_position(self):
        """ Returns 5' position of base in read """
        if self.read.is_reverse:
            pos = self.read.infer_read_length() - self.pileupread.query_position
        else:
            pos = self.pileupread.query_position + 1
        return pos

    @property
    def distance_from_5_prime(self):
        """ Returns distance from 5' end of fragment. In cases
        where reads overlap, 5' distance from read may not equal 5'
        distance from end of fragment """
        return min(
            [self.read_position, abs(self.insert_size) - self.read_position + 1]
        )

    @property
    def target_site(self):
        """ Returns string representing target site """
        return f"{self.chrom}:{self.target_pos}"

    @property
    def fragment_id(self):
        """ Returns string representing fragment_id """
        return f"{self.chrom}:{self.start + 1}-{self.end}"

    @property
    def molecule_class(self):
        """ Return string representing molecule class """
        if self.base == self.target_ref_base:
            return "REF"
        elif self.base == self.target_alt_base:
            return "ALT"
        else:
            return "OTHER"


def parse_maf(maf_file):
    """ Reads in MAF file and yields a Variant namedtuple """
    infile = pd.read_table(maf_file)
    headers = infile.columns
    Variant = namedtuple("VariantSite", headers)
    for idx, row in infile.iterrows():
        yield Variant(*row)


def pass_filter(pileupread):
    """ Make sure read can be used """
    if (
        not pileupread.is_del
        and not pileupread.is_refskip
        and not pileupread.alignment.is_secondary
        and not pileupread.alignment.is_supplementary
        and pileupread.alignment.is_paired
        and not pileupread.alignment.is_unmapped
        and not pileupread.alignment.mate_is_unmapped
    ):
        return True
    return False


def check_discordance(molecule, other_molecule):
    """ Check discordance between two molecules and update """
    molecule.overlap = True
    other_molecule.overlap = True
    if molecule.base != other_molecule.base:
        molecule.discordant = True
        molecule.base = "N"
        other_molecule.discordant = True
        other_molecule.base = "N"


def target_molecule_generator(
    variant_site, bam, mapping_quality, base_quality, output_both=False
):
    """ For a given variant site, a pileup will be created from the BAM.
    Each read overlapping the variant site will be yielded if it
    passes filters
    Output both: output both reads if both overlap base """
    chrom = str(variant_site.Chromosome)
    start = int(variant_site.Start_position)
    end = start
    ref_base = variant_site.Reference_Allele
    alt_base = variant_site.Tumor_Seq_Allele2
    molecule_dict = {}
    for pileupread in bbutils.iterate_pileup_reads(
        bam,
        chrom,
        start - 1,
        end,
        stepper="nofilter",
        min_mapping_quality=mapping_quality,
        min_base_qual=base_quality,
        truncate=True,
        max_depth=1000000,
        noun="pileup reads",
        verb="processed",
        log_every=100000,
    ):
        read_name = pileupread.alignment.query_name
        if pass_filter(pileupread):
            molecule = TargetMolecule(pileupread, end, ref_base, alt_base)
            if read_name not in molecule_dict:
                molecule_dict[read_name] = [molecule]
            else:
                other_molecule = molecule_dict[read_name][0]
                check_discordance(molecule, other_molecule)
                molecule_dict[read_name].append(molecule)

    for _, molecules in molecule_dict.items():
        if output_both:
            for molecule in molecules:
                yield molecule
        else:
            yield molecules[0]


@click.command(context_settings=dict(help_option_names=["-h", "--help"]))
@click.option(
    "-m",
    "--maf_file",
    type=click.Path(exists=True),
    help="MAF describing variants",
    required=True,
)
@click.option(
    "-b",
    "--bam_file",
    type=click.Path(exists=True),
    help="Input DSC/SSC consensus BAM file",
)
@click.option(
    "-f",
    "--file_list",
    type=click.Path(exists=True),
    help="Run on list of BAM files",
)
@click.option(
    "-o",
    "--output",
    type=click.Path(writable=True),
    help="Output file  [default: BAM prefix]",
    default=None,
)
@click.option(
    "--mapq",
    type=int,
    help="Minimum mapping quality",
    default=60,
    show_default=True,
)
@click.option(
    "--baseq",
    type=int,
    help="Minimum base quality",
    default=89,
    show_default=True,
)
def create_summary_file(bam_file, maf_file, output, mapq, baseq, file_list):
    """ Summarize consensus molecules at variant sites. """
    if file_list:
        bam_files = [x.strip() for x in open(file_list, "r").readlines()]
        outputs = [
            bam.split("/")[-1].replace(".bam", "_pileup_summary.txt")
            for bam in bam_files
        ]
    else:
        bam_files = [bam_file]
        outputs = [
            output
            if output
            else bam_file.split("/")[-1].replace(".bam", "_pileup_summary.txt")
        ]

    for bam_file, output in zip(bam_files, outputs):
        print(f"Working on: {bam_file}")
        outfile = open(output, "w")
        headers = [
            "family_id",
            "family_size",
            "umi",
            "fragment_id",
            "target_site",
            "molecule_class",
            "read_position",
            "distance_from_5_prime",
            "discordant",
            "mismatch_count",
            "N_count",
        ]
        outfile.write("\t".join(headers) + "\n")

        bam = pysam.AlignmentFile(bam_file)
        for variant_site in parse_maf(maf_file):
            for target_molecule in target_molecule_generator(
                variant_site, bam, mapq, baseq
            ):
                outfile.write(
                    "\t".join(
                        str(target_molecule.__getattribute__(x))
                        for x in headers
                    )
                    + "\n"
                )
        total_reads = bbutils.iterate_pileup_reads.counts
        print(f"Processed {total_reads} consensus reads")
        bam.close()


if __name__ == "__main__":
    create_summary_file()
