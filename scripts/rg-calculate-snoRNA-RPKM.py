#!/usr/bin/env python
"""
Based on annotations calculate RPKM values for each snoRNA
and filter all that falls below given quantile.

RPKM = (10^9 * C)/(N * L)

where:
    C = Number of reads mapped to a gene
    N = Total mapped reads in the experiment (library size)
    L = Length of the feature (in this case snoRNA length)
"""

__date_ = "2014-10-01"
__author__ = "Rafal Gumienny"
__email__ = "r.gumienny@unibas.ch"
__license__ = "GPL"

# imports
import sys
import time
import csv
import pandas as pd
from Bio import SeqIO
from collections import Counter
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
                    help="Part of the library that is annotated as snoRNA")
parser.add_argument("--output",
                    dest="output",
                    help="Output file in tab format.")
parser.add_argument("--library",
                    dest="library",
                    required=True,
                    help="Library from which the annotations were generated (in bed format)")
parser.add_argument("--snoRNAs",
                    dest="snoRNAs",
                    required=True,
                    help="BED file with snoRNAs")
parser.add_argument("--quantile",
                    dest="quantile",
                    type=float,
                    default=0.25,
                    help="Quantile for the expression cut-off, defaults to 0.25")
parser.add_argument("--type",
                    dest="type",
                    default="CD",
                    choices=("CD", "HACA"),
                    help="Type of snoRNA, defaults to CD")

try:
    options = parser.parse_args()
except Exception, e:
    parser.print_help()

# redefine a functions for writing to stdout and stderr to save some writting
syserr = sys.stderr.write
sysout = sys.stdout.write


def main():
    """Main logic of the script"""
    library_count = 0
    with open(options.library) as f:
        for line in f:
            library_count += 1

    df = pd.read_table(options.input, header=None)
    df[8] = [i.split(";")[0] for i in df[8]]
    ndf = pd.DataFrame({'count': Counter(df[8].tolist())})
    with open(options.snoRNAs) as f:
        snor_lengths = {rec[3].split(":")[-1]: float(int(rec[2]) - int(rec[1])) for rec in csv.reader(f, delimiter="\t")}
    ndf['snor_length'] = [snor_lengths[i] for i in ndf.index]
    ndf['library_size'] = float(library_count)
    ndf['rpkm'] = (ndf['count'].astype(float)*1000000000)/(ndf['snor_length']*ndf['library_size'])
    ndf = ndf[ndf['rpkm']>=ndf['rpkm'].quantile(options.quantile)]
    ndf = ndf.sort('rpkm', ascending=False)
    rpkms = ndf['rpkm']
    with open(options.output, "w") as out:
        if options.type == "HACA":
            for idx, row in rpkms.iteritems():
                out.write("%s_stem1\t%f\n" % (idx, row))
                out.write("%s_stem2\t%f\n" % (idx, row))
        elif options.type == "CD":
            for idx, row in rpkms.iteritems():
                out.write("%s\t%f\n" % (idx, row))
        else:
            raise Exception("Not such a snoRNA type: %s" % options.type)
    # rpkms.to_csv(options.output, sep="\t", header=None)



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
