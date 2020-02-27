#!/usr/bin/env python

from collections import Counter
from datetime import datetime
import argparse
import glob
import gzip
import os
import sys


file_path = os.path.abspath(__file__)


def parse_cl_arguments():
    parser = argparse.ArgumentParser(
        usage="""

Process sample:
Countbarcodes -f1 <fastq1> -f2 <fastq2> -o <outputprefix> [-c]

Process batch:
Countbarcodes -b <batchfile> [-c]

Process directory:
Countbarcodes -d <directory> [-c]
"""
    )

    parser.add_argument(
        "-f1",
        "--fastq1",
        metavar="",
        help="Fastq file 1/2 (gz or uncompressed)",
    )
    parser.add_argument(
        "-f2",
        "--fastq2",
        metavar="",
        help="Fastq file 2/2 (gz or uncompressed)",
    )
    parser.add_argument(
        "-o",
        "--outputfile",
        metavar="",
        help="Output file prefix (<output>_counts.txt)",
    )
    parser.add_argument(
        "-c",
        "--counts",
        metavar="",
        type=int,
        default=500,
        help=("Top number of hits to return. " "0 = all hits [500]"),
    )
    parser.add_argument(
        "-s",
        "--shell",
        action="store_true",
        help="Create shell script for launching with UGER",
    )
    batchparser = parser.add_argument_group("Process batch")
    batchparser.add_argument(
        "-b",
        "--batchfile",
        metavar="",
        help=("Batch file: include headers for " "fastq1, fastq2, sampleid"),
    )
    directoryparser = parser.add_argument_group("Process directory")
    directoryparser.add_argument(
        "-d",
        "--directory",
        metavar="",
        help=("Process all unmatched barcode files in" " a given directory"),
    )
    return parser.parse_args()


def collect_barcodes(fastq1, fastq2):
    print(f"Working on {fastq1} and {fastq2}")
    if fastq1.endswith("gz"):
        fq1 = gzip.open(fastq1, "rt")
    else:
        fq1 = open(fastq1, "r")
    if fastq2.endswith("gz"):
        fq2 = gzip.open(fastq2, "rt")
    else:
        fq2 = open(fastq2, "r")
    next(fq1)
    next(fq2)
    count = 0
    count_every = 1000000
    time_start = datetime.now()
    for i, (code1, code2) in enumerate(zip(fq1, fq2)):
        if i == 0 or i % 4 == 0:
            count += 1
            if count % count_every == 0:
                print(
                    f"Processed {count} reads so far. Time for last {count_every}: {datetime.now() - time_start}"
                )
                time_start = datetime.now()
            barcode = code1.strip() + "_" + code2.strip()
            yield barcode
    fq1.close()
    fq2.close()


def return_top_hits(barcode_generator, hits):
    if hits == 0:
        hits = None
    barcode_counter = Counter(barcode_generator)
    top_hits = barcode_counter.most_common(hits)
    return top_hits


def write_top_hits(top_hits, outputfile):
    with open(outputfile, "w") as outfile:
        outfile.write("barcode\tcount\n")
        for (barcode, count) in top_hits:
            outfile.write("{}\t{}\n".format(barcode, count))
    return outputfile


def process_batch_file(batchfile):
    listfiles = []
    with open(batchfile, "r") as infile:
        next(infile)  # skip header
        for line in infile:
            f1, f2, s = line.strip().split("\t")
            listfiles.append((f1, f2, s))
    return listfiles


def process_directory(directory):
    dirpath = os.path.abspath(directory)
    filepath = os.path.join(dirpath, "*unmatched.barcode*")
    filelist = sorted(glob.glob(filepath))
    listfiles = []
    if len(filelist) % 2 != 0:
        sys.exit("Not all fastq files have pairs in directory... exiting")
    for f1, f2 in zip(filelist[::2], filelist[1::2]):
        lanenumber = os.path.basename(f1)[0]
        s = "Lane" + lanenumber
        listfiles.append((f1, f2, s))
    return listfiles


def main():
    args = parse_cl_arguments()
    batchfile = args.batchfile
    fastq1 = args.fastq1
    fastq2 = args.fastq2
    counts = args.counts
    directory = args.directory

    if not fastq1 and not batchfile and not directory:
        sys.exit("Specifiy your inputs or use -h, --help")

    if args.batchfile:
        listfiles = process_batch_file(batchfile)
    elif args.directory:
        listfiles = process_directory(directory)
    else:
        listfiles = [(args.fastq1, args.fastq2, args.outputfile)]

    if args.shell:
        with open("count_barcodes.sh", "w") as outfile:
            for (fastq1, fastq2, samplename) in listfiles:
                outputfile = samplename + "_counts.txt"
                script = f"python {file_path} -f1 {fastq1} -f2 {fastq2} -o {samplename} -c {args.counts}"
                outfile.write(script + ";\n")
            sys.exit()
    for (fastq1, fastq2, samplename) in listfiles:
        outputfile = samplename + "_counts.txt"
        barcode_generator = collect_barcodes(fastq1, fastq2)
        top_hits = return_top_hits(barcode_generator, counts)
        write_top_hits(top_hits, outputfile)
        print("Counts for {} created successfully".format(samplename))


if __name__ == "__main__":
    main()
