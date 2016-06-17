#!/usr/bin/env python
"""
convert unmapped sequences to fasta
"""

__date_ = "2014-07-02"
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
                    help="Coma separated list of files")
parser.add_argument("--output",
                    dest="output",
                    required=True,
                    help="Output name")

try:
    options = parser.parse_args()
except Exception, e:
    parser.print_help()

# redefine a functions for writing to stdout and stderr to save some writting
syserr = sys.stderr.write
sysout = sys.stdout.write


def main():
    """Main logic of the script"""
    with open(options.output, 'w') as out:
        for unmapped_file in options.input.split(","):
            with open(unmapped_file) as f:
                if options.verbose:
                    syserr("Processing %s\n" % unmapped_file)
                for row in csv.reader(f, delimiter="\t"):
                    if row[0] != 'id' and len(row[1]) >= 25:
                        out.write(">%s_%s\n%s\n" % (row[0], row[2], row[1]))

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
