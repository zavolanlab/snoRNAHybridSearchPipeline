#!/usr/bin/env python
"""
Calculate probability of snoRNA methylation being functional
"""

__date__ = "2015-05-31"
__author__ = "Rafal Gumienny"
__email__ = "r.gumienny@unibas.ch"
__license__ = "GPL"

# imports
import sys
import time
import pandas as pd
import numpy as np
from contextlib import contextmanager
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
                    help="Input file in tab format. Defaults to sys.stdin.")
parser.add_argument("--output",
                    dest="output",
                    required=True,
                    help="Output file in tab format. Defaults to sys.stdout.")
parser.add_argument("--accessibility",
                    dest="accessibility",
                    required=True,
                    help="File with calculated accessibility")
parser.add_argument("--flanks",
                    dest="flanks",
                    required=True,
                    help="File with calculated flanks composition")
parser.add_argument("--model",
                    dest="model",
                    required=True,
                    help="Statsmodel binary file with the model for snoRNA")

try:
    options = parser.parse_args()
except Exception, e:
    parser.print_help()

# redefine a functions for writing to stdout and stderr to save some writting
syserr = sys.stderr.write
sysout = sys.stdout.write


def main():
    """Main logic of the script"""
    accessibility = pd.read_table(options.accessibility, index_col=0)
    flanks = pd.read_table(options.flanks, index_col=0)
    df = pd.read_table(options.input, header=None)
    df.index = [",".join(str(a) for a in [i, j, k, l]) for i, j, k, l in zip(df[0], df[1], df[2], df[3])]
    df.columns = ['chrom',
                  'snoRNAs',
                  'beg',
                  'end',
                  'strand',
                  'score',
                  'logratio',
                  'count',
                  'modified_nucleotide',
                  'alignment',
                  'sitespec',
                  'snorspec']
    df['logsitespec'] = np.log(df['sitespec'])

    ndf = df.join(accessibility)
    ndf = ndf.join(flanks)
    ndf['const'] = 1.0
    ndf = ndf.dropna()
    model = pd.read_pickle(options.model)
    features_for_model = ['const', 'score', 'logsitespec', 'flanksA', 'ContraScore']
    ndf['Probability'] = model.predict(ndf[features_for_model].astype(np.float64), transform=False)
    ndf['Modification'] = ndf['beg'] + 3
    names_to_keep = ['chrom',
                     'snoRNAs',
                     'beg',
                     'end',
                     'strand',
                     'count',
                     'score',
                     'logsitespec',
                     'flanksA',
                     'ContraScore',
                     'Probability',
                     'Modification',
                     'modified_nucleotide',
                     'alignment']
    ndf = ndf.sort("Probability", ascending=False)[names_to_keep]
    ndf.to_csv(options.output, sep='\t', index=False)




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
