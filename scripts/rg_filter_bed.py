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
parser.add_argument("--filter-multimappers",
                    dest="filter_multimappers",
                    action="store_true",
                    default=False,
                    help="Filter chimeras that can be mapped to multiple places in the genome (with exception of mapping to cannonical targets)")


# redefine a functions for writing to stdout and stderr to save some writting
syserr = sys.stderr.write
sysout = sys.stdout.write


def main(options):
    """Main logic of the script"""
    df = pd.read_table(options.input, names=['chr', 'beg', 'end', 'name', 'score', 'strand'])
    df['mutpos'] = [ i.split(":")[-1] for i in df.name ]
    df['id'] = [":".join(i.split(":")[:-1])for i in df.name]
    df['true_chrom'] = [i.split("|")[1] for i in df.chr]
    canonical_targets = ['RNA18S', 'RNA28S', 'RNA5.8S', "U1", "U1.2", "U1.3", "U1.4", "U1.5", "U1.7",
                         "U1.9", "U1.10", "U2", "U2.1", "U2.2", "U2.3", "U4", "U4.2", "U4atac.1",
                         "U4atac.2", "U4atac.3", "U5", "U5.2", "U5.3", "U5.4", "hsa_U5.5", "U6",
                         "U6.2", "U6.3", "U6.4", "U6.5", "U6.6", "U6.7", "U6atac.1", "U7", "U7.2",
                         "U11", "U12"]
    with open(options.output, 'w') as o:
        for name, group in df.groupby("id"):
            tmpdf_score = group[group.score == group.score.max()]
            is_canonical_target = False
            if any(tmpdf_score.true_chrom.isin(canonical_targets)):
                tmpdf_score = tmpdf_score[tmpdf_score.true_chrom.isin(canonical_targets)]
                is_canonical_target = True


            # number_of_targets = len(tmpdf_score.groupby(["chr", "beg"]))
            if options.filter_multimappers:
                if not is_canonical_target and len(tmpdf_score.groupby(["chr", "beg"]).size()) > 1:
                    continue
            for gname, gdf in tmpdf_score.groupby(["chr", "beg"]):
                # o.write("\t".join(str(i) for i in gdf.iloc[0].tolist()[:-3]) + "\t" + str(number_of_targets) + "\n")
                o.write("\t".join(str(i) for i in gdf.iloc[0].tolist()[:-3]) + "\n")

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
