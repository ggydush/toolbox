#!/usr/bin/env python
import argparse
import sys
import os

import bbutils

fname = os.path.basename(__file__)


def parse_cl_arguments():
    parser = argparse.ArgumentParser(usage=f"{fname} bedfile > output.fasta")
    parser.add_argument(
        "bedfile", help="Input bedfile", nargs="?", default=os.getcwd()
    )
    parser.add_argument(
        "-c",
        "--coordinate",
        metavar="",
        help="Get sequence of coordinate chr:start-end [start, end)",
    )
    parser.add_argument(
        "-r",
        "--reference",
        default="/seq/references/Homo_sapiens_assembly19/v1/Homo_sapiens_assembly19.fasta",
        metavar="",
        help="Reference fasta file for converting BED [Broad HG19]",
    )
    return parser


def reverse_complement(sequence):
    complements = {"A": "T", "T": "A", "G": "C", "C": "G"}
    return "".join(complements.get(base, "N") for base in reversed(sequence))


def main():
    parser = parse_cl_arguments()
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit()

    args = parser.parse_args()
    bedfile = args.bedfile
    reference = args.reference
    coordinate = args.coordinate
    fasta = bbutils.FastaReader(reference)

    if coordinate:
        chrom, start, stop = coordinate.strip().replace(":", "-").split("-")
        start, stop = int(start), int(stop)
        print(fasta.get(chrom, start, stop))
        sys.exit()

    cords = [x.strip() for x in open(bedfile, "r").readlines()]
    seq_dict = {}
    for cord in cords:
        info = cord.split()
        chrom, start, end = info[0:3]
        start, end = int(start), int(end)
        seq = fasta.get(chrom, start, end)
        if len(info) >= 6 and info[5] == "-":
            seq = reverse_complement(seq)
        name = f"{chrom}:{start + 1}-{end}"
        seq_dict[name] = seq

    bbutils.write_fasta(seq_dict, sys.stdout)


if __name__ == "__main__":
    main()
