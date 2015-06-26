#!/usr/bin/env python
"""
Divide plexy output into positives and negatives set
"""

__date_ = "2014-09-26"
__author__ = "Rafal Gumienny"
__email__ = "r.gumienny@unibas.ch"
__license__ = "GPL"

# imports
import re
import sys
import time
import numpy as np
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
                    required=True,
                    help="Input file in Tab format.")
parser.add_argument("--output",
                    dest="output",
                    help="Output file in Tab format.")
parser.add_argument("--threshold",
                    dest="threshold",
                    type=float,
                    default=-1.0,
                    help="Threshold for the site, defaults to -1.0")

try:
    options = parser.parse_args()
except Exception, e:
    parser.print_help()

# redefine a functions for writing to stdout and stderr to save some writting
syserr = sys.stderr.write
sysout = sys.stdout.write


def main():
    """Main logic of the script"""
    df = pd.read_table(options.input, header=None)
    df = df.dropna()
    df = df[df[8] <= options.threshold]
    snor_count = []
    for read_id in df[3]:
        c = float(read_id.split(":")[0].split("-")[-1])
        snor_count.append(1.0/c)
    df['count'] = snor_count
    df['modification'] = [re.search("[ACTGactg]m", i).group()[0] if not pd.isnull(i) else np.nan for i in df[10]]
    df = df.groupby([0, 6, 12]).agg({14: max, 'count': sum, 8: min, 5: max, 15: np.mean, 16: np.mean,
                                     "modification": lambda x: ":".join(str(i) for i in set(x)),
                                     9: lambda x: [str(i) for i in set(x)][0]}).reset_index()
    df['ratio'] = np.log(df['count']/df[14])
    # ndf = df.groupby([0, 12]).agg({'ratio': np.max, 8: np.min, 6: lambda x: ":".join(str(i) for i in x), 5: max}).reset_index()
    df['beg'] = (df[12] - 3).astype(int)
    df['end'] = (df[12] + 4).astype(int)
    # columns: chrom, snoRNA, beg, end, strand, Plexy, log Ratio, count,
    # modified site, alignment
    df[[0, 6, 'beg', 'end', 5, 8, 'ratio', 'count', 'modification', 9, 15, 16]].sort([0, 8]).to_csv(options.output, sep="\t", header=None, index=None)



if __name__ == '__main__':
    try:
        if options.verbose:
            start_time = time.time()
            start_date = time.strftime("%d-%m-%Y at %H:%M:%S")
            syserr("############## Started script on %s ##############\n" % start_date)
        df = main()
        if options.verbose:
            syserr("### Successfully finished in %i seconds, on %s ###\n" % (time.time() - start_time, time.strftime("%d-%m-%Y at %H:%M:%S")))
    except KeyboardInterrupt:
        syserr("Interrupted by user after %i seconds!\n" % (time.time() - start_time))
        sys.exit(-1)
