#!/usr/bin/env python
"""
From the bed from FilterBed step get the positions of the found target sites
in terms of real chromosomes not clusters.
"""

__date__ = "2015-05-27"
__author__ = "Rafal Gumienny"
__email__ = "r.gumienny@unibas.ch"
__license__ = "GPL"

# imports
import sys
import time
import pandas as pd
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
                    default=sys.stdin,
                    help="Input file in special bed format. Defaults to sys.stdin.")
parser.add_argument("--output",
                    dest="output",
                    default=sys.stdout,
                    help="Output file in special bed format. Defaults to sys.stdout.")

try:
    options = parser.parse_args()
except Exception, e:
    parser.print_help()

# redefine a functions for writing to stdout and stderr to save some writting
syserr = sys.stderr.write
sysout = sys.stdout.write


def main():
    """Main logic of the script"""
    mapped = pd.read_table(options.input, names=['chr', 'beg', 'end', 'name', 'score', 'strand'])
    mapped.index = [":".join(i.split(":")[:-1]) for i in mapped.name]
    mapped['true_beg'] = [get_beg(int(i.split('|')[3]),int(i.split('|')[4]),j,k,i.split('|')[2]) for i,j,k in zip(mapped.chr, mapped.beg, mapped.end)]
    mapped['true_end'] = [get_end(int(i.split('|')[3]),int(i.split('|')[4]),j,k,i.split('|')[2]) for i,j,k in zip(mapped.chr, mapped.beg, mapped.end)]
    mapped['true_strand'] = [i.split('|')[2] for i in mapped.chr]
    mapped['true_chr'] = [i.split('|')[1] for i in mapped.chr]

    mapped = mapped[['true_chr', 'true_beg', 'true_end', 'name', 'score', 'true_strand']]
    mapped.columns = ['chr', 'beg', 'end', 'name', 'score', 'strand']
    df = mapped
    df.to_csv(options.output, sep="\t", index=False, header=False)



def get_end(cluster_beg, cluster_end, beg, end, strand):
    if strand == "+":
        return cluster_beg + end - 1
    else:
        return cluster_end - beg


def get_beg(cluster_beg, cluster_end, beg, end, strand):
    if strand == "+":
        return cluster_beg - 1 + beg
    else:
        return cluster_end - end


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
