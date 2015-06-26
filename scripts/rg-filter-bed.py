#!/usr/bin/env python
"""
Filter bed file based on the alignment
score/number of reads in cluster/number of mutations
"""

__date_ = "2014-06-30"
__author__ = "Rafal Gumienny"
__email__ = "r.gumienny@unibas.ch"
__license__ = "GPL"

# imports
import sys
import time
import pandas as pd
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
                    help="Input bed file with special fields")
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
    df = pd.read_table(options.input, names=['chr', 'beg', 'end', 'name', 'score', 'strand'])
    df['mutpos'] = [ i.split(":")[-1] for i in df.name ]
    df['id'] = [":".join(i.split(":")[:-1])for i in df.name]
    with open(options.output, 'w') as o:
        for name, group in df.groupby("id"):
            tmpdf_score = group[group.score == group.score.max()]
            for gname, gdf in tmpdf_score.groupby(["chr", "beg"]):
                o.write("\t".join(str(i) for i in gdf.iloc[0].tolist()[:-2]) + "\n")

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
