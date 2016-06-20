#!/usr/bin/env python
"""
Make statistic, prepare plots and evaluate thresholds
"""

__date__ = "2015-02-23"
__author__ = "Rafal Gumienny"
__email__ = "r.gumienny@unibas.ch"
__license__ = "GPL"

# imports
import matplotlib
matplotlib.use('Agg')
import os
import sys
import time
import pandas as pd
import numpy as np
import pylab as pl
import seaborn as sns
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
parser.add_argument("--dir",
                    dest="dir",
                    default="Plots",
                    help="Directory to store the plots , defaults to Plots")
parser.add_argument("--length",
                    dest="length",
                    type=int,
                    default=15,
                    help="Threshold for length of the target site, defaults to 15")
parser.add_argument("--fpr",
                    dest="fpr",
                    type=float,
                    default=0.05,
                    help="False positive rate threshold, defaults to 0.05")


# redefine a functions for writing to stdout and stderr to save some writting
syserr = sys.stderr.write
sysout = sys.stdout.write


def main(options):
    """Main logic of the script"""
    pl.rcParams['figure.figsize'] = (14, 10)
    pl.rcParams['ytick.labelsize'] = 20
    pl.rcParams['xtick.labelsize'] = 20
    pl.rcParams['axes.labelsize'] = 23
    pl.rcParams['legend.fontsize'] = 20
    sns.set_style('ticks')
    c1, c2, c3 = sns.color_palette("Set1", 3)
    names_string = "read_id\tname\tscore\tbeg\tend\twhich_end\ttarget\tqbeg\tqend\tisambig\t"
    names_string += "matches:mismatches:identity\tshuf_snor\tshuf_score\tshuf_target_len\tshuf_qbeg\tshuf_qend"
    names = names_string.split("\t")
    df = pd.read_table(options.input, header=None, names=names)
    df["target_length"] = df["target"].str.len()
    ndf = df[df["target_length"] >= options.length]
    #
    # Statistics
    #
    library_size = len(df.index)
    found_sites = len(df["score"].dropna().index)
    found_sites_shuf = len(df["shuf_score"].dropna().index)
    found_sites_thres = len(df["score"].dropna().index)

    df["length"] = df["end"] - df["beg"]
    df["identity"] = [float(i.split(":")[-1]) if not pd.isnull(i) else i for i in df["matches:mismatches:identity"]]
    df["mismatch"] = [float(i.split(":")[-2]) if not pd.isnull(i) else i for i in df["matches:mismatches:identity"]]

    fig, ax  = plot_stats(df,
                          "identity",
                          "Identity",
                          "Identity between snoRNA and read vs. score of alignment")
    fig.savefig(os.path.join(options.dir, "search_anchors_identity_vs_score.pdf"))
    fig.savefig(os.path.join(options.dir, "search_anchors_identity_vs_score.png"), dpi=300)

    fig, ax = plot_stats(df,
                         "length",
                         "Length",
                         "Length of alignment vs. score of alignment")
    fig.savefig(os.path.join(options.dir, "search_anchors_length_vs_score.pdf"))
    fig.savefig(os.path.join(options.dir, "search_anchors_length_vs_score.png"), dpi=300)

    fig, ax = plot_stats(df,
                         "mismatch",
                         "Number of mismatches",
                         "Number of mismatches vs. score of alignment")
    fig.savefig(os.path.join(options.dir, "search_anchors_length_vs_score.pdf"))
    fig.savefig(os.path.join(options.dir, "search_anchors_length_vs_score.png"), dpi=300)

    fig, ax = pl.subplots()
    df["score"].dropna().hist(bins=range(20,70), lw=5, histtype="step", ax=ax, label="Reads")
    df["shuf_score"].dropna().hist(bins=range(20,70), lw=5, histtype="step", ax=ax, label="Shuffled")
    sns.despine(offset=10)
    pl.xlabel("Score")
    pl.ylabel("Count")
    pl.legend()
    pl.title("Score histogram for reads and shuffled control")
    fig.savefig(os.path.join(options.dir, "search_anchors_score_histogram.pdf"))
    fig.savefig(os.path.join(options.dir, "search_anchors_score_histogram.png"), dpi=300)

    props = df.groupby("shuf_score").size()/df.groupby("score").size()
    props = props.fillna(0.0)
    fig, ax = pl.subplots()
    props.plot(lw=5, ax=ax)
    ax.hlines(options.fpr, min(props.index), max(props.index), color='r', linestyle="--")
    sns.despine()
    ax.set_xlabel("Score")
    ax.set_ylabel("False Positive Rate")
    # ax.set_xlim(20, 70)
    fig.savefig(os.path.join(options.dir, "search_anchors_porportions.pdf"))
    fig.savefig(os.path.join(options.dir, "search_anchors_porportions.png"), dpi=300)

    score_threshold = None
    for ind, value in props.iteritems():
        if value <= options.fpr:
            score_threshold = ind
            break

    ambigous = df.groupby("isambig").size()
    ends = df.groupby("which_end").size()
    with open(options.output, "w") as out:
        out.write("Library size\t%i\n" % (library_size))
        out.write("Number of found sites\t%i\n" % (found_sites))
        out.write("Number of found sites with target site over %i nuc\t%i\n" % (options.length, found_sites_thres))
        out.write("Number of unique sites (over threshold)\t%i\n" % (ambigous["Unique"]))
        try:
            out.write("Number of ambiguous sites (over threshold)\t%i\n" % (ambigous["Ambiguous"]))
        except KeyError:
            out.write("Number of ambiguous sites (over threshold)\t0\n")
        out.write("Number of 5' end sites (over threshold)\t%i\n" % (ends['5p']))
        out.write("Number of 3' end sites (over threshold)\t%i\n" % (ends['3p']))
        out.write("Averaged length of alignment:\t%.2f +- %.2f\n" % (df.length.mean(), df.length.std()))
        out.write("Averaged number of mismatches in alignment:\t%.2f +- %.2f\n" % (df.mismatch.mean(), df.mismatch.std()))
        out.write("Averaged identity of alignment:\t%.2f +- %.2f\n" % (df.identity.mean(), df.identity.std()))
        if score_threshold:
            out.write("Score threshold\t%f\n" % score_threshold)
        else:
            out.write("Score threshold\t0.0\n")


def plot_stats(df, col, ylabel, title):
    fig, ax = pl.subplots()
    ndf = df.groupby("score").agg({col: [np.mean, np.std]})
    ax.errorbar(ndf[col]["mean"].index, ndf[col]["mean"].values,
                            ndf[col]["std"].values, linestyle="--", marker='o',
                            markersize=10)
    pl.grid(True)
    sns.despine()
    ax.set_xlabel("Score")
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    return fig, ax

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
            syserr("############## Started script on %s ##############\n" %
                   start_date)
        main(options)
        if options.verbose:
            syserr("### Successfully finished in %i seconds, on %s ###\n" %
                   (time.time() - start_time,
                    time.strftime("%d-%m-%Y at %H:%M:%S")))
    except KeyboardInterrupt:
        syserr("Interrupted by user after %i seconds!\n" %
               (time.time() - start_time))
        sys.exit(-1)
