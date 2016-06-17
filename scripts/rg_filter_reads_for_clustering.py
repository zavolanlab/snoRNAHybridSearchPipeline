#!/usr/bin/env python
"""
Filter reads based on annotation in the last column
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
                    help="Input table")
parser.add_argument("--output",
                    dest="output",
                    required=True,
                    help="Output table")
parser.add_argument("--annotations",
                    dest="annotations",
                    default="None",
                    help="Coma separated list of annotations to consider, defaults to None")


# redefine a functions for writing to stdout and stderr to save some writting
syserr = sys.stderr.write
sysout = sys.stdout.write


def main(options):
    """Main logic of the script"""
    if options.annotations != "None":
        annots = options.annotations.split(",")
        with open(options.output, 'w') as out, open(options.input) as inputfile:
            for row in csv.reader(inputfile, delimiter='\t'):
                if row[-1] in annots:
                    out.write("\t".join(row) + "\n")
    else:
        with open(options.output, 'w') as out, open(options.input) as inputfile:
            for row in csv.reader(inputfile, delimiter='\t'):
                out.write("\t".join(row) + "\n")


if __name__ == '__main__':
    try:
        try:
            options = parser.parse_args()
        except Exception, e:
            parser.print_help()
            sys.exit()
        if options.verbose:
            start_time = time.time()
            start_date = time.strftime("%d-%m-%Y at %H:%M:%S")
            syserr("############## Started script on %s ##############\n" % start_date)
        main(options)
        if options.verbose:
            syserr("### Successfully finished in %i seconds, on %s ###\n" % (time.time() - start_time, time.strftime("%d-%m-%Y at %H:%M:%S")))
    except KeyboardInterrupt:
        syserr("Interrupted by user after %i seconds!\n" % (time.time() - start_time))
        sys.exit(-1)
