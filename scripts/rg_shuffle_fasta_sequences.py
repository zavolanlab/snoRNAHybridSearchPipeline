#!/usr/bin/env python
"""
Shuffle fasta sequences in the file
"""

__date_ = "2014-08-14"
__author__ = "Rafal Gumienny"
__email__ = "r.gumienny@unibas.ch"
__license__ = "GPL"

# imports
import sys
import time
import ushuffle
from Bio import SeqIO
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
                    help="Input fasta file")
parser.add_argument("--output",
                    dest="output",
                    default="output_shuffled.fa",
                    help="Output fasta file , defaults to output_shuffled.fa")
parser.add_argument("--let-size",
                    dest="let_size",
                    type=int,
                    default=2,
                    help="Let size to preserve, defaults to 2")

try:
    options = parser.parse_args()
except Exception, e:
    parser.print_help()

# redefine a functions for writing to stdout and stderr to save some writting
syserr = sys.stderr.write
sysout = sys.stdout.write


def main():
    """Main logic of the script"""
    with open(options.input) as infile, open(options.output, 'w') as outfile:
        for rec in SeqIO.parse(infile, 'fasta'):
            shuffled_seq = ushuffle.shuffle(str(rec.seq), len(rec.seq), options.let_size)
            outfile.write(">%s\n%s\n" % (rec.id, shuffled_seq))

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
