#!/usr/bin/env python
"""
Convert result to asmbed and in the same time extend sequences to be equal
desired length
"""

__date_ = "2014-08-11"
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
                    default="output.asmbed",
                    help="Output asmbed file , defaults to output.asmbed")
parser.add_argument("--length",
                    dest="length",
                    type=int,
                    default=50,
                    help="Desired read length, defaults to 50")


# redefine a functions for writing to stdout and stderr to save some writting
syserr = sys.stderr.write
sysout = sys.stdout.write


def main(options):
    """Main logic of the script"""
    with open(options.input) as infile, open(options.output, 'w') as outfile:
        for row in csv.reader(infile, delimiter='\t'):
            assert int(row[2]) - int(row[1]) <= 50
            lind, upind = get_indices(options.length, int(row[1]), int(row[2]))
            outfile.write("%s\t%s\t%s\t%i\t%i\n" % (row[3],
                                                    row[0],
                                                    row[5],
                                                    lind + 1,
                                                    upind))

def get_indices(desired_length, beg, end):
    """Calculate upper and lower index based on what size we want and what are the coordinates"""
    lower_index = beg - ((desired_length - (end - beg))//2 + (desired_length - (end - beg))%2)
    upper_index = end + ((desired_length - (end - beg))//2)
    return lower_index, upper_index


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
