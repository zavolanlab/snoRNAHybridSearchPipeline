#!/usr/bin/env python
"""
"""

__date_ = "2014-06-05"
__author__ = "Rafal Gumienny"
__email__ = "r.gumienny@unibas.ch"
__license__ = "GPL"

# imports
import sys
import csv
import time
import swalign
import pandas as pd
from ushuffle import shuffle
from Bio import SeqIO
from collections import OrderedDict, defaultdict
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
                    help="Input table")
parser.add_argument("--match",
                    dest="match",
                    type=int,
                    default=2,
                    help="Match score, defaults to 2")
parser.add_argument("--mismatch",
                    dest="mismatch",
                    type=int,
                    default=-5,
                    help="Mismatch penalty, defaults to -5")
parser.add_argument("--gap-open",
                    dest="gap_open",
                    type=int,
                    default=-6,
                    help="Open gap penalty, defaults to -6")
parser.add_argument("--gap-extend",
                    dest="gap_extend",
                    type=int,
                    default=-4,
                    help="Gap extension penalty, defaults to -4")
parser.add_argument("--output",
                    dest="output",
                    help="Output table")


try:
    options = parser.parse_args()
except Exception, e:
    parser.print_help()

# redefine a functions for writing to stdout and stderr to save some writting
syserr = sys.stderr.write
sysout = sys.stdout.write


def main():
    """Main logic of the script"""
    #
    # prepare swalign parameters
    #
    scoring = swalign.NucleotideScoringMatrix(match=options.match,
                                              mismatch=options.mismatch)
    sw = swalign.LocalAlignment(scoring,
                                gap_penalty=options.gap_open,
                                gap_extension_penalty=options.gap_extend)

    #
    # Iterate reads to find matches
    #
    # counter = 0
    # names = ["coords",
             # "snorid",
             # "seq",
             # "hg19seq",
             # "hg18seq"]
    df = pd.read_table(options.input, header=None)
    hg19_score = []
    hg18_score = []
    assembly = []
    for idx, row in df.iterrows():
        # import ipdb; ipdb.set_trace()
        try:
            h19 = sw.align(row[6], row[20]).score
        except AttributeError:
            h19 = 0.0
        try:
            h18 = sw.align(row[6], row[21]).score
        except AttributeError:
            h18 = 0.0
        hg19_score.append(h19)
        hg18_score.append(h18)

        if h19 > h18:
            assembly.append("hg19")
        elif h19 < h18:
            assembly.append("hg18")
        else:
            assembly.append("ambiguous")
    df["hg19score"] = hg19_score
    df["hg18score"] = hg18_score
    df["assembly"] = assembly
    df.to_csv(options.output, sep="\t", index=False)


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
