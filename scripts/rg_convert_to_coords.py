#!/usr/bin/env python
"""
Convert result to coordinate file
"""

__date_ = "2014-08-11"
__author__ = "Rafal Gumienny"
__email__ = "r.gumienny@unibas.ch"
__license__ = "GPL"

# imports
import sys
import csv
import time
from Bio import SeqIO
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
                    help="Input result file")
parser.add_argument("--sequences",
                    dest="sequences",
                    required=True,
                    help="File with sequences")
parser.add_argument("--output",
                    dest="output",
                    default="coords.tab",
                    help="Output coordinate file , defaults to coords.tab")

try:
    options = parser.parse_args()
except Exception, e:
    parser.print_help()

# redefine a functions for writing to stdout and stderr to save some writting
syserr = sys.stderr.write
sysout = sys.stdout.write


def main():
    """Main logic of the script"""
    with open(options.sequences) as s:
        seqs = {str(rec.id).split("|")[0]: str(rec.seq).upper() for rec in SeqIO.parse(s, 'fasta')}

    with open(options.input) as infile, open(options.output, 'w') as outfile:
        for row in csv.reader(infile, delimiter='\t'):
            outfile.write("%s\t%s\t%s\t%s\t%s\n" % (row[0],
                                              row[0].split(":")[-2],
                                              row[3],
                                              row[4],
                                              seqs[row[0]]))

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
