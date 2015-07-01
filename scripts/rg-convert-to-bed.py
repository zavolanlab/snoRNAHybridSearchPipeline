#!/usr/bin/env python
"""
Convert Probability results into bed for annotations
"""

__date__ = "2015-05-31"
__author__ = "Rafal Gumienny"
__email__ = "r.gumienny@unibas.ch"
__license__ = "GPL"

# imports
import sys
import time
import pandas as pd
from contextlib import contextmanager
from argparse import ArgumentParser, RawTextHelpFormatter

parser = ArgumentParser(description=__doc__, formatter_class=RawTextHelpFormatter)
parser.add_argument("-v",
                    "--verbose",
                    dest="verbose",
                    action="store_true",
                    default=False,
                    help="Be loud!")
parser.add_argument("--input",
                    dest="input",
                    required=True,
                    help="Input file")
parser.add_argument("--output",
                    dest="output",
                    required=True,
                    help="Output file")

try:
    options = parser.parse_args()
except Exception, e:
    parser.print_help()

# redefine a functions for writing to stdout and stderr to save some writting
syserr = sys.stderr.write
sysout = sys.stdout.write


def main():
    """Main logic of the script"""
    df = pd.read_table(options.input)
    df["beg"] = df["Modification"] - 1
    df["end"] = df["Modification"]
    columns = ["chrom",
               "beg",
               "end",
               "snoRNAs",
               "count",
               "strand",
               "score",
               "logsitespec",
               "modified_nucleotide",
               "alignment",
               "Probability",
               "Modification",
               ]
    df[columns].to_csv(options.output, header=None, index=None, sep="\t")


if __name__ == '__main__':
    try:
        if options.verbose:
            start_time = time.time()
            start_date = time.strftime("%d-%m-%Y at %H:%M:%S")
            syserr("############## Started script on %s ##############\n" %
                   start_date)
        main()
        if options.verbose:
            syserr("### Successfully finished in %i seconds, on %s ###\n" %
                   (time.time() - start_time,
                    time.strftime("%d-%m-%Y at %H:%M:%S")))
    except KeyboardInterrupt:
        syserr("Interrupted by user after %i seconds!\n" %
               (time.time() - start_time))
        sys.exit(-1)
