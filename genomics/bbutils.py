"""
These are common utility functions used for Blood Biopsy Team
"""
import gzip
import os

import pandas as pd


try:
    from bam_iterator import (
        iterate_reads,
        iterate_pileup_reads,
        iterate_read_pairs,
        progress_logger,
    )
except SyntaxError:
    pass


class FastaReader:
    def __init__(self, fasta, fasta_idx=None):
        """Reference fasta reader

        Parameters
        ----------
        fasta : path
            Path to reference fasta
        fasta_idx : path, optional
            Path to reference fasta index

        Raises
        ------
        FileNotFoundError
            Fasta index is required, raises error if not found
        """
        if not fasta_idx:
            fasta_idx = fasta + ".fai"

        if not os.path.exists(fasta_idx):
            raise FileNotFoundError("Fasta file must be indexed (.fasta.fai)")

        self.__fasta_idx_df = pd.read_table(
            fasta_idx,
            header=None,
            dtype={"contig": str},
            names=[
                "contig",
                "length",
                "start",
                "bases_per_line",
                "bytes_per_line",
            ],
        )
        self.__fasta_idx_df.set_index("contig", inplace=True)
        self.__fasta_fh = open(fasta, "r")
        self.fasta_path = fasta

    def get(self, contig, seq_start, seq_end):
        """Get sequence from fasta

        Parameters
        ----------
        contig : str
            Chromosome
        seq_start : int
            Start position [0-indexed]
        seq_end : int
            End position [1-indexed]

        Returns
        -------
        str
            Returns sequence string
        """
        contig = str(contig)
        if seq_start < 0:
            raise ValueError("Sequence start must be greater than 0")
        bases_per_line = self.__fasta_idx_df.loc[contig, "bases_per_line"]
        bytes_per_line = self.__fasta_idx_df.loc[contig, "bytes_per_line"]
        seq_length = seq_end - seq_start
        contig_start = self.__fasta_idx_df.loc[contig, "start"]
        num_newlines = seq_start // bases_per_line

        adj_seq_start = contig_start + seq_start + num_newlines
        max_seq_newlines = (seq_length // bytes_per_line) + 2
        max_seq_length = seq_length + max_seq_newlines

        self.__fasta_fh.seek(adj_seq_start)
        raw_seq = self.__fasta_fh.read(max_seq_length)
        raw_seq = raw_seq.replace("\n", "")
        seq = raw_seq[:seq_length]

        return seq

    def close(self):
        self.__fasta_fh.close()


class VCFReader(object):
    """Reader for VCF files."""

    def __init__(self, filename):
        gzipped = filename.endswith(".gz")
        if gzipped:
            self.__reader = gzip.open(filename, mode="rt")
        else:
            self.__reader = open(filename, "rt")

        line = next(self.__reader)
        while line.startswith("##"):
            line = next(self.__reader)

    def __iter__(self):
        return self

    def __next__(self):
        raw_line = next(self.__reader)
        line = raw_line.strip().split("\t")

        chrom = line[0]
        pos = int(line[1])

        if line[2] == ".":
            ident = None
        else:
            ident = line[2]

        ref = line[3]
        alt = line[4].split(",")

        try:
            qual = int(line[5])
        except ValueError:
            try:
                qual = float(line[5])
            except ValueError:
                qual = None

        if line[6] == ".":
            filt = []
        else:
            filt = line[6].split(";")

        record = VCFRecord(chrom, pos, ident, ref, alt, qual, filt, raw_line)

        return record

    def close(self):
        self.__reader.close()

    next = __next__  # Python 2/3 compatibility


class VCFRecord(object):
    """Record from a VCF file"""

    def __init__(self, chrom, pos, ident, ref, alt, qual, filt, raw_line):
        self.chrom = chrom
        self.pos = pos
        self.id = ident
        self.ref = ref
        self.alt = alt
        self.qual = qual
        self.filter = filt
        self.raw = raw_line

    def __str__(self):
        self.id = "." if not self.id else self.id
        self.qual = "." if not self.qual else self.qual
        self.filter = ["."] if not self.filter else self.filter

        return "{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}".format(
            self.chrom,
            self.pos,
            self.id,
            self.ref,
            ",".join(self.alt),
            self.qual,
            ";".join(self.filter),
        )


class FastqRecord:
    """Fastq Record"""

    def __init__(self, seq_id, seq, quals, raw_lines):
        self.seq_id = seq_id
        self.seq = seq
        self.quals = quals
        self.raw_lines = raw_lines

    def __str__(self):
        return "{0}\t{1}\t{2}".format(self.seq_id, self.seq, self.quals)


def reverse_complement(sequence):
    """ Reverse complement a given sequence. Unknown bases are kept in output """
    compdict = {"A": "T", "T": "A", "G": "C", "C": "G"}
    revseq = reversed(sequence)
    return "".join([compdict.get(base.upper(), base) for base in revseq])


def vcf_to_df(vcf_file):
    """Convert a VCF file into a data frame."""
    vcf_reader = VCFReader(vcf_file)

    mut_matrix = []
    for mut in vcf_reader:
        for alt in mut.ALT:
            mut_matrix.append([mut.CHROM, mut.POS, mut.ID, mut.REF, alt])

    df = pd.DataFrame(
        mut_matrix,
        columns=["contig", "position", "id", "ref_allele", "alt_allele"],
    )

    return df


def interval_list_to_df(interval_list_file):
    """Convert interval list to data frame."""
    df = pd.read_table(
        interval_list_file,
        header=None,
        names=["contig", "start", "end", "strand", "id"],
        comment="@",
        dtype={"contig": str},
    )

    return df


def maf_to_df(maf_file):
    """Convert a MAF file into a data frame."""
    try:
        maf_df = pd.read_table(maf_file, comment="#")
    except UnicodeDecodeError:
        maf_df = pd.read_table(maf_file, comment="#", encoding="latin-1")

    return maf_df


def write_fasta(seq_dict, output_file, line_length=80):
    """ Write fasta file containing sequences in seq_dict (name: sequence) """
    if isinstance(output_file, str):
        outfile = open(output_file, "w")
    else:
        outfile = output_file

    for seq_id in seq_dict:
        outfile.write(">{0}\n".format(seq_id))
        seq = seq_dict[seq_id]

        # split seq into lines
        seq_lst = [
            seq[i : i + line_length] for i in range(0, len(seq), line_length)
        ]

        outfile.write("\n".join(seq_lst) + "\n")

    if isinstance(output_file, str):
        outfile.close()
