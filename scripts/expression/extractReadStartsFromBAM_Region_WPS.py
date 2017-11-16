#!/usr/bin/env python
"""Original author: Martin Kircher / mkircher@uw.edu / *03.06.2014"""

from sys import stderr
from os.path import exists
import gzip
from pysam import Samfile
from random import random
from optparse import OptionParser
from collections import defaultdict
from bx.intervals.intersection import Intersecter, Interval

VALID_CHROMS = set(map(str,list(range(1,23))+["X","Y"]))

def is_soft_clipped(cigar):
    for op, count in cigar:
        if op in [4, 5, 6]:
            return True
    return False

def aln_length(cigarlist):
    tlength = 0
    for operation, length in cigarlist:
        if operation == 0 or operation == 2 or operation == 3 or operation >= 6:
            tlength += length
    return tlength

def log_action(message, verbosity=True):
    if verbosity:
        print(message, file=stderr)

def get_arguments():
    parser = OptionParser()
    parser.add_option("-i", "--input", help="Regions transcript file")
    parser.add_option("-l", "--length_sr", help="Full read length (default=76)", default=76, type="int")
    parser.add_option("-m", "--merged", help="Assume reads are merged",default=False,action="store_true")
    parser.add_option("-t", "--trimmed", help="Assume reads are trimmed",default=False,action="store_true")
    parser.add_option("-w", "--protection", help="Base pair protection (default=120)",default=120,type="int")
    parser.add_option("-o", "--outfile", help="Output file(s) prefix")
    parser.add_option("-e", "--empty", help="Keep files for empty blocks",default=False,action="store_true")
    parser.add_option("--min_ins_size", help="Minimum read length threshold to consider (default=None)",default=-1,type="int")
    parser.add_option("--max_ins_size", help="Minimum read length threshold to consider (default=None)",default=-1,type="int")
    parser.add_option("--max_length", help="Maximum insert size (default=1000)",default=1000,type="int")
    parser.add_option("--downsample", help="Ratio to down sample reads",default=None,type="float")
    parser.add_option("-v","--verbose", help="Be verbose",default=False,action="store_true")
    options, args = parser.parse_args()
    options.outfile = options.outfile.strip("""\'""")
    bamfiles = [bf.strip("""\'""") for bf in args]
    for bamfile in bamfiles:
        if not exists(bamfile):
            raise IOError("File not found: {}".format(bamfile))
        if not exists(bamfile + ".bai"):
            if not exists(bamfile.replace(".bam", ".bai")):
                raise IOError("Index for file {} not found".format(bamfile))
    return options, bamfiles

def get_ins_sizes(options):
    if 0 < options.min_ins_size < options.max_ins_size:
        log_mask = "Using min/max length cutoffs: {}/{}"
        log_action(log_mask.format(options.min_ins_size, options.max_ins_size))
        return options.min_ins_size, options.max_ins_size
    else:
        return None, None

def valid_regions(line_iterator):
    for line in line_iterator:
        cid, chrom, start, end, strand = line.split()
        if (chrom in VALID_CHROMS) and (int(start) >= 1):
            yield cid, chrom, int(start), int(end), strand

def get_prefix(references):
    for tchrom in references:
        if tchrom.startswith("chr"):
            return "chr"
    else:
        return ""

def filter_reads(sam, chrom, region_start, region_end, protection, options):
    span_start = region_start - protection - 1
    span_end = region_end + protection + 1
    for read in sam.fetch(chrom, span_start, span_end):
        if (not options.downsample) or (random() < options.downsample):
            if not (read.is_duplicate or read.is_qcfail or read.is_unmapped):
                if not is_soft_clipped(read.cigar):
                    yield read

def is_valid_paired(r, region_start, protection, min_ins_size, max_ins_size):
    if r.is_paired and (not r.mate_is_unmapped) and (r.rnext == r.tid):
        span_start = region_start - protection - 1
        if r.is_read1 or (r.is_read2 and r.pnext+r.qlen < span_start):
            if r.isize != 0:
                if (min_ins_size == None):
                    return True
                elif min_ins_size < abs(r.isize) < max_ins_size:
                    return True
    return False

def is_valid_single(r):
    if (min_ins_size == None):
        return True
    elif min_ins_size < aln_length(r.cigar) < max_ins_size:
        return True
    else:
        return False

def as_merged(r, options):
    return options.merged or r.qname.startswith("M_")

def as_trimmed(r, options):
    if (options.trimmed or r.qname.startswith("T_")):
        return r.qlen <= options.length_sr - 10
    return False

def update_pr_region(pos_range, rstart, rend, region_start, region_end):
    for i in range(rstart, rend + 1):
        if region_start <= i <= region_end:
            pos_range[i][0] += 1

def update_pr_at(pos_range, at, region_start, region_end):
    if region_start <= at <= region_end:
        pos_range[at][1] += 1

def get_reads_and_ranges(bamfiles, cid, chrom, region_start, region_end, strand, options):
    pos_range = defaultdict(lambda: [0,0])
    filtered_reads = Intersecter()
    for bamfile in bamfiles:
        log_action("Reading {}".format(bamfile), options.verbose)
        with Samfile(bamfile, "rb") as sam:
            prefix = get_prefix(sam.references)
            for read in filter_reads(sam, prefix+chrom, region_start, region_end, protection, options):
                if is_valid_paired(read, region_start, protection, min_ins_size, max_ins_size):
                    rstart = min(read.pos, read.pnext) + 1
                    rend = rstart + abs(read.isize) - 1
                    filtered_reads.add_interval(Interval(rstart, rend))
                    update_pr_region(pos_range, rstart, rend, region_start, region_end)
                    update_pr_at(pos_range, rstart, region_start, region_end)
                    update_pr_at(pos_range, rend, region_start, region_end)
                elif is_valid_single(read):
                    rstart = read.pos + 1
                    rend = rstart + aln_length(read.cigar) - 1
                    filtered_reads.add_interval(Interval(rstart, rend))
                    update_pr_region(pos_range, rstart, rend, region_start, region_end)
                    if as_merged(read, options) or as_trimmed(read, options):
                        update_pr_at(pos_range, rstart, region_start, region_end)
                        update_pr_at(pos_range, rend, region_start, region_end)
                    elif read.is_reverse:
                        update_pr_at(pos_range, rend, region_start, region_end)
                    else:
                        update_pr_at(pos_range, rstart, region_start, region_end)
    return filtered_reads, pos_range

def calculate_wps(filtered_reads, pos_range, cid, chrom, region_start, region_end, strand, options):
    cov_sites = 0
    wps_list = []
    for pos in range(region_start, region_end + 1):
        rstart, rend = pos - protection, pos + protection
        gcount, bcount = 0, 0
        for read in filtered_reads.find(rstart, rend):
            if (read.start > rstart) or (read.end < rend):
                bcount += 1
            else:
                gcount += 1
        cov_count, start_count = pos_range[pos]
        cov_sites += cov_count
        wps_list.append([chrom, pos, cov_count, start_count, gcount - bcount])
    return wps_list, cov_sites

if __name__ == "__main__":
    options, bamfiles = get_arguments()
    min_ins_size, max_ins_size = get_ins_sizes(options)
    protection = options.protection//2
    log_action("Reading {}".format(options.input), options.verbose)
    with open(options.input) as region_file:
        for region in valid_regions(region_file):
            filtered_reads, pos_range = get_reads_and_ranges(bamfiles, *region, options)
            wps_list, cov_sites = calculate_wps(filtered_reads, pos_range, *region, options)
            if cov_sites or options.empty:
                if region[4] == "-":
                    wps_list = reversed(wps_list)
                with gzip.open(options.outfile%region[0], "wt") as wps_handle:
                    for line in wps_list:
                        print(*line, sep="\t", file=wps_handle)
