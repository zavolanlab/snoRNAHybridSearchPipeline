#!/usr/bin/env python
"""
Compare output results with original data
"""

__date_ = "2014-07-01"
__author__ = "Rafal Gumienny"
__email__ = "r.gumienny@unibas.ch"
__license__ = "GPL"

# imports
import sys
import csv
import time
from argparse import ArgumentParser

parser = ArgumentParser(description=__doc__)
parser.add_argument("-v",
                    "--verbose",
                    dest="verbose",
                    action="store_true",
                    default=False,
                    help="Be loud!")
parser.add_argument("--input",
                    dest="input",
                    required=True,
                    help="Bed file with special fields")
parser.add_argument("--only-chrom",
                    dest="only_chrom",
                    action="store_true",
                    default=False,
                    help="If there is a bed file with only chromosome information use this flag")

try:
    options = parser.parse_args()
except Exception, e:
    parser.print_help()

# redefine a functions for writing to stdout and stderr to save some writting
syserr = sys.stderr.write
sysout = sys.stdout.write


def main():
    """Main logic of the script"""
    match_count = 0
    all_count = 0
    with open(options.input) as bed:
        if not options.only_chrom:
            for contig, start, end, name, score, strand, miRNA in csv.reader(bed, delimiter='\t'):
                my_name, my_chr, my_strand, my_start, my_end, my_1, my_2 = contig.split("|")
                or_chr, or_start, or_mir, or_name, or_mut = name.split(":")
                if my_chr == or_chr and or_start >= my_start and or_start <= my_end:
                    match_count += 1
                all_count += 1
        else:
            for contig, start, end, name, score, strand, miRNA in csv.reader(bed, delimiter='\t'):
                or_chr, or_start, or_mir, or_name, mymir, or_mut = name.split(":")
                if or_chr == contig:
                    match_count += 1
                else:
                    print contig, start, end, name, score, strand
                all_count += 1
    if options.verbose:
        syserr("Found %i sites of wich %i match the original (%.2f%%)\n" % (all_count,
                                                                            match_count,
                                                                            100*(match_count/float(all_count))))

    # Cluster_246110_57_4|chr10|-|101156587|101156743|0|0 50  70  chr10:101156673:hsa-miR-24:human07603:14    40  +
    # Cluster_246177_75_38|chr10|-|101370722|101370896|0|0    89  109 chr10:101370787:hsa-miR-10a:human01320:8    40  +
    # Cluster_246259_119_18|chr10|-|101636627|101636845|0|0   105 124 chr10:101636721:hsa-miR-28-5p:human03016:9  38  +

if __name__ == '__main__':
    try:
        if options.verbose:
            start_time = time.time()
            start_date = time.strftime("%d-%m-%Y at %H:%M:%S")
            syserr("############## Started script on %s ##############\n" % start_date)
        main()
        if options.verbose:
            syserr("### Successfully finished in %i seconds, on %s ###\n" % (time.time() - start_time, time.strftime("%d-%m-%Y at %H:%M:%S")))
    except KeyboardInterrupt:
        syserr("Interrupted by user after %i seconds!\n" % (time.time() - start_time))
        sys.exit(-1)
