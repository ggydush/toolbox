""" Provides functions for iterating over BAMs while logging progress

Examples:

from bam_iterator import iterate_reads, iterate_pileup_reads, iterate_read_pairs

bam = pysam.AlignmentFile(bamfile)

for read in iterate_reads(bam):
    # do stuff

for pileupread in iterate_pileup_reads(bam, '1', 100, 101):
    # do stuff

for read1, read2 in iterate_read_pairs(bam):
    # do stuff

Total counts can be accessed by function attributes:
total_counts = iterate_read_pairs.counts

"""
from collections import defaultdict
from datetime import datetime
import functools
import logging

import pysam


logging.basicConfig(
    format="%(levelname)s\t%(asctime)s\t%(name)s\t%(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
    level=logging.INFO,
)
logger = logging.getLogger("bam_iterator")


def _get_position(obj):
    """ Return short region string (chrom:start) of object
    if it's a read, otherwise return '*' """
    try:
        if isinstance(obj, pysam.libcalignedsegment.PileupRead):
            obj = obj.alignment
        elif isinstance(obj, tuple):
            obj = obj[0]
        last_position = f"{obj.reference_name}:{obj.reference_start}"
    except AttributeError:
        last_position = "*"
    return last_position


def _get_time_difference(start_time, elapsed=False):
    """ Get time between start_time and now. If elapsed, return h:mm:ss,
    otherwise return seconds """
    difference = datetime.now() - start_time
    if elapsed:
        return str(difference).split(".")[0]
    else:
        return f"{difference.total_seconds():0.3f}"


def progress_logger(func):
    """ Decorate generators with this function to keep counts/timers in log file

    Simple Example:
    @progress_logger
    def test_function(iterable, *args, **kwargs):
        for x in iterable:
            yield x

    for x in test_function(range(1, 100), log_every=10):
        pass
    """

    @functools.wraps(func)
    def wrapped(
        *args, verb="processed", noun="reads", log_every=1000000, **kwargs
    ):
        wrapped.timer = getattr(wrapped, "timer", datetime.now())
        gen = func(*args, **kwargs)
        for entry in gen:
            wrapped.counts += 1
            if wrapped.counts % log_every == 0:
                # Time everything
                time_for_last_str = _get_time_difference(
                    wrapped.timer, elapsed=False
                )
                elapsed_time_str = _get_time_difference(
                    wrapped.start, elapsed=True
                )
                wrapped.timer = datetime.now()

                # Output to log
                last_position = _get_position(entry)
                output_string = (
                    f"{verb}\t"
                    f"{wrapped.counts:,} {noun}\t"
                    f"Elapsed time: {elapsed_time_str}s\t"
                    f"Time for last {log_every:,}: {time_for_last_str}s\t"
                    f"Last position: {last_position}"
                )
                logger.info(output_string)
            yield entry

    wrapped.counts = 0
    wrapped.start = datetime.now()
    return wrapped


@progress_logger
def iterate_reads(bam, *args, **kwargs):
    """ Fetch reads from BAM file. Accepts all arguments of bam.fetch().

    param: bam (pysam.AlignmentFile) - bam object to iterate over

    Additional kwargs:
    param: verb (str) - verb to print in log
    param: noun (str) - noun to print in log
    param: log_every (int) - print a log statement every log_every records

    yields read
    """
    for read in bam.fetch(*args, **kwargs):
        yield read


@progress_logger
def iterate_pileup_reads(bam, *args, **kwargs):
    """ Get pileupreads from BAM file. Accepts all arguments of bam.pileup().

    param: bam (pysam.AlignmentFile) - bam object to iterate over

    Additional kwargs:
    param: verb (str) - verb to print in log
    param: noun (str) - noun to print in log
    param: log_every (int) - print a log statement every log_every records

    yields pileupread
    """
    for pileupcolumn in bam.pileup(*args, **kwargs):
        for pileupread in pileupcolumn.pileups:
            yield pileupread


def iterate_read_pairs(
    bam, *args, verb="processed", noun="read pairs", log_every=100000, **kwargs
):
    """ Get read pairs from BAM file. Accepts all arguments of bam.fetch().

    param: bam (pysam.AlignmentFile) - bam object to iterate over

    Additional kwargs:
    param: verb (str) - verb to print in log
    param: noun (str) - noun to print in log
    param: log_every (int) - print a log statement every log_every records

    yields tuple(read1, read2)
    """
    for read1, read2 in _read_pair_generator(
        bam, *args, verb=verb, noun=noun, log_every=log_every, **kwargs
    ):
        yield read1, read2
    iterate_read_pairs.counts = _read_pair_generator.counts


def _keep_read(read):
    """ Only keep reads if they are mapped and paired """
    if not read.is_unmapped:
        if not read.is_secondary and not read.is_supplementary:
            if read.is_paired:
                return read
    return None


@progress_logger
def _read_pair_generator(bam, *args, **kwargs):
    """ Generate read pairs in a BAM file. Accepts all arguments of bam.fetch()

    param: bam (pysam.AlignmentFile) - bam object to iterate over

    yields tuple(read1, read2)
    """
    read_dict = defaultdict(lambda: [None, None])
    for read in bam.fetch(*args, **kwargs):
        if not _keep_read(read):
            continue
        qname = read.query_name
        if qname not in read_dict:
            if read.is_read1:
                read_dict[qname][0] = read
            else:
                read_dict[qname][1] = read
        else:
            if read.is_read1:
                yield read, read_dict[qname][1]
            else:
                yield read_dict[qname][0], read
            del read_dict[qname]


if __name__ == "__main__":
    pass
