#!/usr/bin/env python

"""

:Author: Martin Kircher
:Contact: mkircher@uw.edu
:Date: *09.01.2015
"""

import sys
from os.path import exists
from optparse import OptionParser
from pysam import Fastafile, Samfile
from collections import defaultdict
from tqdm import tqdm

table = str.maketrans("TGCA", "ACGT")


def is_soft_clipped(cigar):
    for op, count in cigar:
        if op in [4,5,6]:
            return True
    return False


def aln_len(cigarlist):
    tlength = 0
    for operation, length in cigarlist:
        if operation == 0 or operation == 2 or operation == 3 or operation >= 6:
            tlength += length
    return tlength


def has_index(bamfile):
    return exists(bamfile.replace(".bam", ".bai")) or exists(bamfile + ".bai")


def get_args():
    parser = OptionParser()
    parser.add_option(
        "-r", "--region",
        dest="region",
        help="Region to be looked up (default = 12:34,443,233-34,453,733)",
        default="12:34,443,233-34,453,733"
    )
    parser.add_option(
        "-f", "--fasta",
        dest="reference",
        help="Fasta index reference genome"
    )
    parser.add_option(
        "-m", "--max",
        dest="maxreads",
        help="Maximum number of reads to consider (default 0 = Off)",
        default=0,
        type="int"
    )
    parser.add_option(
        "-d", "--dinucleotides",
        dest="dinucleotides",
        help="Additionaly extract dinucleotide patterns (def Off)",
        default=False,
        action="store_true"
    )
    parser.add_option(
        "-v", "--verbose",
        dest="verbose",
        help="Turn debug output on",
        default=False,
        action="store_true"
    )
    parser.add_option(
        "-o", "--outfile",
        dest="outfile",
        help="Write output to file (def STDOUT)",
        default="/dev/stdout"
    )
    options, args = parser.parse_args()
    if not exists(options.reference) or not exists(options.reference + ".fai"):
        raise IOError("Fasta indexed reference genome is not available")
    options.region = options.region.strip("""\'""").strip()
    if options.region.upper() == "ALL":
        options.region = ""
    return options, [bamfile for bamfile in args if exists(bamfile)]


options, args = get_args()
reference = Fastafile(options.reference)

chrom, start, end, outchrom = None, None, None, None

if len(options.region):
    try:
        chrom = options.region.split(":")[0]
        start, end = map(int, options.region.split(":")[1].replace(",","").split("-"))
        outchrom = chrom
        outchrom = outchrom.replace("chr", "")
        if outchrom.startswith("gi|"):
            NCfield = outchrom.split("|")[-2]
            if NCfield.startswith("NC"):
                outchrom = "%d"%(int(NCfield.split(".")[0].split("_")[-1]))
        if options.verbose:
            sys.stderr.write("Coordinates parsed: Chrom %s Start %d End %d\n"%(chrom,start,end))
    except:
        sys.stderr.write("Error: Region not defined in a correct format!\n")
        sys.exit()

data = {
    True: dict((nuc, defaultdict(int)) for nuc in list("ACGT")),
    False: dict((nuc, defaultdict(int)) for nuc in list("ACGT"))
}

dinucleotides = {
    True: {},
    False: {}
}

if options.dinucleotides:
    for base1 in list("ACGT"):
        for base2 in list("ACGT"):
            dinucleotides[True][base1+base2]=defaultdict(int)
            dinucleotides[False][base1+base2]=defaultdict(int)


def is_meaningful(alignment):
    return not (
        alignment.is_unmapped or alignment.is_duplicate or
        alignment.is_qcfail or (alignment.cigar == None) or
        is_soft_clipped(alignment.cigar)
    )

readCount = 0
for bamfile in args:
    if (chrom != None) and not has_index(bamfile):
        continue
    with Samfile(bamfile, "rb") as input_file:
        BAMreferences = dict(enumerate(input_file.references))
        if chrom == None:
            BAMiter = input_file
        else:
            BAMiter = input_file.fetch(chrom, start, end)
        for alignment in tqdm(filter(is_meaningful, BAMiter), desc=bamfile):
            if chrom == None:
                try:
                    outchrom = BAMreferences[alignment.tid]
                    outchrom = outchrom.replace("chr", "")
                except:
                    continue
            posList = []
            if alignment.is_paired:
                if alignment.mate_is_unmapped:
                    continue
                if alignment.rnext != alignment.tid:
                    continue
                if alignment.is_reverse:
                    posList.append((
                        alignment.pos + aln_len(alignment.cigar),
                        True,
                        alignment.is_read1
                    ))
                else:
                    posList.append((alignment.pos, False, alignment.is_read1))
            else:
                posList.append((alignment.pos, False, not(alignment.is_reverse)))
                posList.append((
                    alignment.pos + aln_len(alignment.cigar),
                    True,
                    alignment.is_reverse
                ))
            readCount += 1
            for pos, orient, threePrime in posList:
                try:
                    seq = reference.fetch(outchrom, pos-20, pos+21)
                except:
                    continue
                if orient:
                    rseq = seq.translate(table)[::-1]
                    for ind, base in enumerate(rseq[1:]):
                        try:
                            data[threePrime][base][ind] += 1
                        except:
                            pass
                        if options.dinucleotides:
                            try:
                                dinucleotides[threePrime][base+seq[ind+1]][ind] += 1
                            except:
                                pass
                else:
                    for ind,base in enumerate(seq[:-1]):
                        try:
                            data[threePrime][base][ind] += 1
                        except:
                            pass
                        if options.dinucleotides:
                            try:
                                dinucleotides[threePrime][base+seq[ind+1]][ind] += 1
                            except:
                                pass
            if (options.maxreads > 0) and (readCount > options.maxreads):
                sys.stderr.write("Maximum number of reads reached.\n")
                break


def get_outstr(threePrime, helper, base, subDict):
    primity = "ThreePrime" if threePrime else "FivePrime"
    tail = [
        "NA" if helper[ind] == 0 else "%.4f"%(x/float(helper[ind]))
        for ind, x in sorted(subDict.items())
    ]
    return [primity, base] + tail


with open(options.outfile, "w") as outfile:
    ranges = list(range(-20, 0)) + list(range(1, 21))
    print("Type", "Base", *ranges, file=outfile)
    for threePrime in [True, False]:
        helper = list(map(lambda k,v:v, sorted(data[threePrime]["A"].items())))
        for base in ["C", "G", "T"]:
            for ind, val in data[threePrime][base].items():
                helper[ind] += val
        for base, subDict in sorted(data[threePrime].items()):
            outstr = get_outstr(threePrime, helper, base, subDict)
            print(*outstr, sep="\t", file=outfile)
    if options.dinucleotides:
        for threePrime in [True, False]:
            helper = list(map(lambda k, v: v, sorted(dinucleotides[threePrime]["AA"].items())))
            for base in dinucleotides[threePrime].keys():
                if base == "AA":
                    continue
                for ind, val in dinucleotides[threePrime][base].items():
                    helper[ind] += val
            for base, subDict in sorted(dinucleotides[threePrime].items()):
                outstr = get_outstr(threePrime, helper, base, subDict)
                print(*outstr, sep="\t", file=outfile)
