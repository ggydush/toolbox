#!/usr/bin/env python
import argparse
import sys


def parse_cl_arguments():
    parser = argparse.ArgumentParser(description="Reverse complement sequence")
    parser.add_argument(
        "sequence", help="Sequence to be reverse complemented (use - for stdin)"
    )
    parser.add_argument(
        "-c",
        "--comp_only",
        action="store_true",
        help="Return complement only (no reverse)",
    )
    return parser


def get_reverse_complement(sequence, comp_only=False):
    complements = {"A": "T", "T": "A", "G": "C", "C": "G"}
    if not comp_only:
        sequence = reversed(sequence)
    return "".join(complements.get(base, "N") for base in sequence)


def main():
    parser = parse_cl_arguments()
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit()
    args = parser.parse_args()

    sequences = [args.sequence]

    if args.sequence == "-":
        sequences = sys.stdin.read().split("\n")

    for sequence in sequences:
        sys.stdout.write(
            get_reverse_complement(sequence, args.comp_only) + "\n"
        )


if __name__ == "__main__":
    main()
