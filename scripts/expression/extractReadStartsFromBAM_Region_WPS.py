#!/usr/bin/env python

"""

:Author: Martin Kircher
:Contact: mkircher@uw.edu
:Date: *03.06.2014
"""

import sys, os
from optparse import OptionParser
import gzip
import pysam
import random

from collections import defaultdict
from bx.intervals.intersection import Intersecter, Interval

def is_soft_clipped(cigar):
    for op, count in cigar:
        if op in [4,5,6]:
            return True
    return False

def aln_length(cigarlist):
    tlength = 0
    for operation, length in cigarlist:
        if operation == 0 or operation == 2 or operation == 3 or operation >= 6:
            tlength += length
    return tlength

parser = OptionParser()
parser.add_option("-i","--input", help="Use regions transcript file (def transcriptAnno.tsv)",default="transcriptAnno.tsv")
parser.add_option("-l","--length_sr", help="Length of full reads (default 76)",default=76,type="int")
parser.add_option("-m","--merged", help="Assume reads are merged (default Off)",default=False,action="store_true")
parser.add_option("-t","--trimmed", help="Assume reads are trimmed (default Off)",default=False,action="store_true")
parser.add_option("-w","--protection", help="Base pair protection assumed for elements (default 120)",default=120,type="int")
parser.add_option("-o","--outfile", help="Outfile prefix (def 'block_%s.tsv.gz')",default='block_%s.tsv.gz')
parser.add_option("-e","--empty", help="Keep files of empty blocks (def Off)",default=False,action="store_true")
parser.add_option("--min_ins_size", help="Minimum read length threshold to consider (def None)",default=-1,type="int")
parser.add_option("--max_ins_size", help="Minimum read length threshold to consider (def None)",default=-1,type="int")
parser.add_option("--max_length", help="Assumed maximum insert size (default 1000)",default=1000,type="int")
parser.add_option("--downsample", help="Ratio to down sample reads (default OFF)",default=None,type="float")
parser.add_option("-v","--verbose", help="Turn debug output on",default=False,action="store_true")
(options, args) = parser.parse_args()

options.outfile = options.outfile.strip("""\'""")
protection = options.protection//2
VALID_CHROMS = set(map(str,list(range(1,23))+["X","Y"]))
min_ins_size, max_ins_size = None, None

if options.min_ins_size > 0 and options.max_ins_size > 0 and options.min_ins_size < options.max_ins_size:
    min_ins_size = options.min_ins_size
    max_ins_size = options.max_ins_size
    sys.stderr.write("Using min/max length cutoffs: %d/%d\n"%(min_ins_size,max_ins_size))

if os.path.exists(options.input):
    if options.verbose: sys.stderr.write("Reading %s\n"%options.input)
    infile = open(options.input)
    for line in infile:
        if options.verbose: sys.stderr.write(line)
        cid, chrom, start, end, strand = line.split() # positions should be 0-based and end non-inclusive
        if chrom not in VALID_CHROMS:
            continue
        region_start, region_end = int(start), int(end)
        if region_start < 1:
            continue
     
        pos_range = defaultdict(lambda:[0,0])
        filtered_reads = Intersecter()
     
        for bamfile in args:
            if options.verbose: sys.stderr.write("Reading %s\n"%bamfile)
            bamfile = bamfile.strip("""\'""")
            if not os.path.exists(bamfile):
                raise IOError("File not found: {}".format(bamfile))
            if not (os.path.exists(bamfile.replace(".bam",".bai")) or os.path.exists(bamfile+".bai")):
                raise IOError("Index for file {} not found".format(bamfile))
            input_file = pysam.Samfile( bamfile, "rb" )
            prefix = ""
            for tchrom in input_file.references:
                if tchrom.startswith("chr"): 
                    prefix = "chr"
                    break
            for read in input_file.fetch(prefix+chrom,region_start-protection-1,region_end+protection+1):
                if read.is_duplicate or read.is_qcfail or read.is_unmapped: continue
                if is_soft_clipped(read.cigar):
                    continue
                if read.is_paired:
                    if read.mate_is_unmapped:
                        continue
                    if read.rnext != read.tid:
                        continue
                    if read.is_read1 or (read.is_read2 and read.pnext+read.qlen < region_start-protection-1):
                        if read.isize == 0:
                            continue
                        if options.downsample != None and random.random() >= options.downsample:
                            continue
                        rstart = min(read.pos,read.pnext)+1 # 1-based
                        lseq = abs(read.isize)
                        rend = rstart+lseq-1 # end included
                        if min_ins_size != None and ((lseq < min_ins_size) or (lseq > max_ins_size)):
                            continue
                        filtered_reads.add_interval(Interval(rstart,rend))
                        for i in range(rstart,rend+1):
                            if i >= region_start and i <= region_end:
                                pos_range[i][0]+=1
                        if rstart >= region_start and rstart <= region_end:
                            pos_range[rstart][1]+=1
                        if rend >= region_start and rend <= region_end:
                            pos_range[rend][1]+=1
                else:
                    if options.downsample != None and random.random() >= options.downsample:
                        continue
                    rstart = read.pos+1 # 1-based
                    lseq = aln_length(read.cigar)
                    rend = rstart+lseq-1 # end included
                    if min_ins_size != None and ((lseq < min_ins_size) or (lseq > max_ins_size)):
                        continue
                    filtered_reads.add_interval(Interval(rstart,rend))
                    for i in range(rstart,rend+1):
                        if i >= region_start and i <= region_end:
                            pos_range[i][0]+=1
                    if ((options.merged or read.qname.startswith('M_')) or ((options.trimmed or read.qname.startswith('T_')) and read.qlen <= options.length_sr-10)):
                        if (rstart >= region_start and rstart <= region_end):
                            pos_range[rstart][1]+=1
                        if rend >= region_start and rend <= region_end:
                            pos_range[rend][1]+=1
                    elif read.is_reverse:
                        if rend >= region_start and rend <= region_end:
                            pos_range[rend][1]+=1
                    else:
                        if (rstart >= region_start and rstart <= region_end):
                            pos_range[rstart][1]+=1
         
        filename = options.outfile%cid
        outfile = gzip.open(filename,'wt')
        cov_sites = 0
        out_lines = []
        for pos in range(region_start, region_end+1):
            rstart, rend = pos-protection, pos+protection
            gcount, bcount = 0, 0
            for read in filtered_reads.find(rstart,rend):
                if (read.start > rstart) or (read.end < rend):
                    bcount += 1
                else:
                    gcount += 1
            cov_count, start_count = pos_range[pos]
            cov_sites += cov_count
            out_lines.append("%s\t%d\t%d\t%d\t%d\n"%(chrom,pos,cov_count,start_count,gcount-bcount))

        if strand == "-": out_lines = out_lines[::-1]
        for line in out_lines: outfile.write(line)
        outfile.close()

        if cov_sites == 0 and not options.empty:
            os.remove(filename)
