#!/usr/bin/env python
"""
Make some plots of the results
"""

__date_ = "2014-08-12"
__author__ = "Rafal Gumienny"
__email__ = "r.gumienny@unibas.ch"
__license__ = "GPL"

# imports
import sys
import csv
import time
import pylab as pl
import numpy as np
import pandas as pd
import seaborn as sns
from collections import Counter
from scipy import stats
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
parser.add_argument("--clustered",
                    dest="clustered",
                    action="store_true",
                    default=False,
                    help="Is the result clustered?")
parser.add_argument("--expressions",
                    dest="expressions",
                    required=True,
                    help="File with miRNA expression")
parser.add_argument("--level",
                    dest="level",
                    type=int,
                    default=0,
                    help="Expression level (in log scale), defaults to 0")
parser.add_argument("--top",
                    dest="top",
                    type=int,
                    default=20,
                    help="Show top mirnas and number of hybrids found, defaults to 20")

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
    if options.clustered:
        ndf = pd.DataFrame({'hybrids': Counter(df[3].tolist())})
        top_mirnas = Counter(df[3].tolist())
    else:
        ndf = pd.DataFrame({'hybrids': Counter(df[6].tolist())})
        top_mirnas = Counter(df[6].tolist())
    ndf.hybrids = np.log(ndf.hybrids)
    exp = {}
    with open(options.expressions) as e:
        for row in csv.reader(e, delimiter='\t'):
            if row[0] != 'Mirnas/Samples':
                exp[row[0].replace("(star)", "*")] = np.log(np.mean([float(value) for value in row[1:]]))

    ndf['expr'] = [exp[i] for i in ndf.index]
    ndf = ndf[ndf['expr'] >= options.level]
    # ndf = ndf[ndf['hybrids'] <= 7.0]
    if options.verbose:
        syserr("Number of miRNAs: %i\n" % len(ndf.index))
        for val in top_mirnas.most_common()[:options.top]:
            syserr("%s\t%i\n" % (val[0], val[1]))
    sns.set(style='white', font='serif')
    grid = sns.JointGrid('hybrids', 'expr', data=ndf, dropna=True, size=8, ratio=10,
                                 space=0.2)
    grid.plot_marginals(sns.distplot, color="firebrick", rug=False,
                                hist=True, kde=False)
    grid.plot_joint(sns.regplot, scatter_kws={"color": "slategray"}, line_kws={"linewidth": 1, "color": "firebrick"})
    grid.annotate(stats.pearsonr, stat="pearson r", template="{stat} = {val:.2g}\np-value = {p:.2g}")
    grid.set_axis_labels('log(Number of hybrids)', 'log(Expression value)')
    pl.savefig("hybrids_vs_expression.png")
    pl.clf()
    sns.set(style='ticks', font='serif')
    sns.distplot(ndf.hybrids, kde=False)
    sns.despine(trim=True)
    sns.offset_spines()
    pl.xlabel("log(Number of hybrids)")
    pl.ylabel("Fraction")
    pl.title("Number of hybrids found per miRNA")
    pl.tight_layout()
    pl.savefig("hist_of_hybrid_number.png")

    pl.clf()
    sns.set(style='ticks', font='serif')
    sns.distplot(ndf.expr, kde=False)
    sns.despine(trim=True)
    sns.offset_spines()
    pl.xlabel("log(Expression value)")
    pl.ylabel("Fraction")
    pl.title("Expression of miRNAs")
    pl.tight_layout()
    pl.savefig("hist_of_expressions.png")

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
