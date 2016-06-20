#!/usr/bin/env python
"""
Append rpkm values to the plexy predictions
"""

__date_ = "2014-10-01"
__author__ = "Rafal Gumienny"
__email__ = "r.gumienny@unibas.ch"
__license__ = "GPL"

# imports
import sys
import csv
import time
import pandas as pd
import numpy as np
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
                    help="Input file in tab format.")
parser.add_argument("--output",
                    dest="output",
                    help="Output file in tab format.")
parser.add_argument("--rpkm",
                    dest="rpkm",
                    required=True,
                    help="File with rpkms of snoRNAs")
parser.add_argument("--annotated-reads",
                    dest="annotated_reads",
                    required=True,
                    help="Mapped reads annotated as snoRNAs")
parser.add_argument("--type",
                    dest="type",
                    default="CD",
                    choices=("CD", "HACA"),
                    help="Type of snoRNAs , defaults to CD")


# redefine a functions for writing to stdout and stderr to save some writting
syserr = sys.stderr.write
sysout = sys.stdout.write


def main(options):
    """Main logic of the script"""
    raw_count = pd.read_table(options.annotated_reads, header=None)
    raw_count[8] = [i.split(";")[0] for i in raw_count[8]]
    raw_count = raw_count.groupby(8).size()
    if options.type == "HACA":
        new_raw_count = {}
        for idx, count in raw_count.iteritems():
            new_raw_count[idx + "_stem1"] = count
            new_raw_count[idx + "_stem2"] = count
        raw_count = pd.Series(new_raw_count)
    ndf = pd.read_table(options.rpkm,
                        names=['rpkm'],
                        index_col=0)
    rpkms = ndf['rpkm']
    df = pd.read_table(options.input, header=None)
    snor_count = []
    for read_id in df[3]:
        c = float(read_id.split(":")[0].split("-")[-1])
        snor_count.append(1.0/c)
    df['count'] = snor_count
    specs = {}
    for name, group in df.groupby(6):
        specs[name] = group.groupby([0, 12]).apply(lambda x: len(x)/group["count"].sum())

    specificity = []
    for chrom, snor, pos in zip(df[0], df[6], df[12]):
        if not pd.isnull(pos):
            specificity.append(specs[snor][chrom][pos])
        else:
            specificity.append(np.nan)

    snor_spec = {}
    for name, group in df.groupby(6):
        snor_spec[name] = len(group)/(group["count"].sum() + raw_count[name])

    df['rpkm'] = [rpkms[i] if i in rpkms else np.nan for i in df[6]]
    df['specificity'] = specificity
    df['snor_spec'] = [snor_spec[snor] for snor in df[6]]
    columns_to_save = [col for col in df.columns if col != 'count']
    df[columns_to_save].to_csv(options.output, header=None, index=None, sep="\t")

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
