#!/usr/bin/env python
"""
Translate coordinates of one rRNA to other rRNA from different species
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
from Bio import SeqIO, AlignIO
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
                    help="Alignment file from Needle (EMBOSS format)")
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
    with open(options.input) as infile:
        alignment = AlignIO.read(infile, 'emboss')

    human_seq = str(alignment[0].seq)
    mouse_seq = str(alignment[1].seq)
    if human_seq[-1] == '-':
        human_seq = human_seq[:-1]
        mouse_seq = mouse_seq[:-1]
    assert len(human_seq) == len(mouse_seq), "Alignments should be equal"
    human_counter = 0
    mouse_counter = 0
    translation = []
    for hpos, mpos in zip(human_seq, mouse_seq):
        if hpos != '-' and mpos != '-':
            human_counter += 1
            mouse_counter += 1
            translation.append((human_counter, mouse_counter))
        elif hpos != '-' and mpos == '-':
            human_counter += 1
            translation.append((human_counter, "NA"))
        elif hpos == '-' and mpos != '-':
            mouse_counter += 1
            translation.append((human_counter, mouse_counter))
    print translation


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
