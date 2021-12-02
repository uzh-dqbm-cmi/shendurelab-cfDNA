#!/usr/bin/env python
"""Original author: Martin Kircher / mkircher@uw.edu / *03.06.2014"""
#from IPython import embed
from sys import stderr
from os.path import exists
from gzip import open as gzopen
from pysam import Samfile
from random import random
from optparse import OptionParser
from collections import namedtuple, defaultdict
from bx.intervals.intersection import Intersecter, Interval
from joblib import Parallel, delayed

ARG_RULES = {
    ("-i", "--input"): {"help": "Regions transcript file"},
    ("-o", "--outfile"): {"help": "Output file(s) prefix"},
    ("--chrom_prefix",): {"default": ""},
    ("-l", "--length_sr"): {
        "help": "Full read length (default: 76)",
        "default": 76,
        "type": "int"
    },
    ("-m", "--merged"): {
        "help": "Assume reads are merged",
        "default": False,
        "action": "store_true"
    },
    ("-t", "--trimmed"): {
        "help": "Assume reads are trimmed",
        "default": False,
        "action": "store_true"
    },
    ("-w", "--protection"): {
        "help": "Base pair protection (default: 120)",
        "default": 120,
        "type": "int"
    },
    ("-j", "--jobs"): {
        "help": "Number of pooled processes",
        "default": 1,
        "type": "int"
    },
    ("-e", "--empty"): {
        "help": "Keep files for empty blocks",
        "default": False,
        "action": "store_true"
    },
    ("--min_ins_size",): {
        "help": "Minimum read length threshold to consider (default: None)",
        "default": -1,
        "type": "int"
    },
    ("--max_ins_size",): {
        "help": "Minimum read length threshold to consider (default: None)",
        "default": -1,
        "type": "int"
    },
    ("--max_length",): {
        "help": "Maximum insert size (default: 1000)",
        "default": 1000,
        "type": "int"
    },
    ("--downsample",): {
        "help": "Ratio to down sample reads",
        "default": None,
        "type": "float"
    },
    ("-v", "--verbose"): {
        "default": False,
        "action": "store_true",
    }
}

# VALID_CHROMS = set(map(str,list(range(1,23))+["X","Y"]))
VALID_CHROMS = {"chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12" ,
                "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX", "chrY"}

Region = namedtuple("Region", [
    "cid", "chrom", "region_start", "region_end", "strand"
])
Read = namedtuple("Read", [
    "cigar", "end", "is_duplicate", "isize", "is_paired", "is_qcfail",
    "is_read1", "is_read2", "is_reverse", "is_unmapped", "mate_is_unmapped",
    "pnext", "pos", "qlen", "qname", "rnext", "start", "tid"
])


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


def get_arguments(arg_rules):
    parser = OptionParser()
    for rule_unnamed, rule_named in arg_rules.items():
        parser.add_option(*rule_unnamed, **rule_named)
    options, args = parser.parse_args()
    options.prot_radius = int(options.protection / 2)
    options.outfile = options.outfile.strip("""\'""")
    bamfile = args[0].strip("""\'""")
    if not exists(bamfile):
        raise IOError("File not found: {}".format(bamfile))
    if not exists(bamfile + ".bai"):
        if not exists(bamfile.replace(".bam", ".bai")):
            raise IOError("Index for file {} not found".format(bamfile))
    if 0 < options.min_ins_size < options.max_ins_size:
        log_mask = "Using min/max length cutoffs: {}/{}"
        log_action(log_mask.format(options.min_ins_size, options.max_ins_size))
    else:
        options.min_ins_size, options.max_ins_size = None, None
    return options, bamfile


def pickleable_region(bam_fetch):
    return [
        Read(*[getattr(read, attr, None) for attr in Read._fields])
        for read in bam_fetch
    ]


def valid_regions(anno, bam, options):
    for line in anno:
        cid, chrom, start, end, strand = line.split()
        region_start, region_end = int(start), int(end)
        span_start = region_start - options.prot_radius - 1
        span_end = region_end + options.prot_radius + 1
        if (chrom in VALID_CHROMS) and (region_start >= 1):
            region = Region(cid, chrom, int(start), int(end), strand)
            bam_fetch = bam.fetch(chrom, span_start, span_end)
            bam_region = pickleable_region(bam_fetch)
            yield region, bam_region


def filter_reads(bam_region, chrom, region_start, region_end, options):
    for read in bam_region:
        if (not options.downsample) or (random() < options.downsample):
            if not (read.is_duplicate or read.is_qcfail or read.is_unmapped):
                if not is_soft_clipped(read.cigar):
                    yield read


def is_valid_paired(r, region_start, options):
    if r.is_paired and (not r.mate_is_unmapped) and (r.rnext == r.tid):
        span_start = region_start - options.prot_radius - 1
        if r.is_read1 or (r.is_read2 and r.pnext+r.qlen < span_start):
            if r.isize != 0:
                if (options.min_ins_size == None):
                    return True
                elif options.min_ins_size < abs(r.isize) < options.max_ins_size:
                    return True
    return False


def is_valid_single(r, options):
    if (options.min_ins_size == None):
        return True
    elif options.min_ins_size < aln_length(r.cigar) < options.max_ins_size:
        return True
    else:
        return False


def as_merged(r, options):
    return options.merged or r.qname.startswith("M_")


def as_trimmed(r, options):
    if (options.trimmed or r.qname.startswith("T_")):
        return r.qlen <= options.length_sr - 10
    return False


def inc_pr(pos_range, rstart, rend, region_start, region_end):
    for i in range(rstart, rend + 1):
        if region_start <= i <= region_end:
            pos_range[i][0] += 1


def inc_pr_at(pos_range, at, region_start, region_end):
    if region_start <= at <= region_end:
        pos_range[at][1] += 1


def get_reads_and_ranges(bam_region, cid, chrom, region_start, region_end, strand, options):
    pos_range = defaultdict(lambda: [0,0])
    filtered_reads = Intersecter()
    read_iterator = filter_reads(
        bam_region, options.chrom_prefix + chrom,
        region_start, region_end, options
    )
    for read in read_iterator:
        if is_valid_paired(read, region_start, options):
            rstart = min(read.pos, read.pnext) + 1
            rend = rstart + abs(read.isize) - 1
            filtered_reads.add_interval(Interval(rstart, rend))
            inc_pr(pos_range, rstart, rend, region_start, region_end)
            inc_pr_at(pos_range, rstart, region_start, region_end)
            inc_pr_at(pos_range, rend, region_start, region_end)
        elif is_valid_single(read, options):
            rstart = read.pos + 1
            rend = rstart + aln_length(read.cigar) - 1
            filtered_reads.add_interval(Interval(rstart, rend))
            inc_pr(pos_range, rstart, rend, region_start, region_end)
            if as_merged(read, options) or as_trimmed(read, options):
                inc_pr_at(pos_range, rstart, region_start, region_end)
                inc_pr_at(pos_range, rend, region_start, region_end)
            elif read.is_reverse:
                inc_pr_at(pos_range, rend, region_start, region_end)
            else:
                inc_pr_at(pos_range, rstart, region_start, region_end)
    return filtered_reads, pos_range


def get_wps(filtered_reads, pos_range, cid, chrom, region_start, region_end, strand, options):
    cov_sites = 0
    wps_list = []
    for pos in range(region_start, region_end + 1):
        rstart, rend = pos - options.prot_radius, pos + options.prot_radius
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


def generate_region_file(bam_region, region, options):
    reads, pos_range = get_reads_and_ranges(bam_region, *region, options)
    wps_list, cov_sites = get_wps(reads, pos_range, *region, options)
    if cov_sites or options.empty:
        if region.strand == "-":
            wps_list = reversed(wps_list)
        with gzopen(options.outfile.replace('*', region.cid), "wt") as wps_handle:
            for line in wps_list:
                print(*line, sep="\t", file=wps_handle)


if __name__ == "__main__":
    options, bamfile = get_arguments(ARG_RULES)
    with Samfile(bamfile,check_sq = False) as bam, open(options.input) as anno:
        Parallel(n_jobs=options.jobs)(
            delayed(generate_region_file)(bam_region, region, options)
            for region, bam_region in valid_regions(anno, bam, options)
        )
