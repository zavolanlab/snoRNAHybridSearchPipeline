#!/usr/bin/env python
"""

"""

__date__ = "2015-03-03"
__author__ = "Rafal Gumienny"
__email__ = "r.gumienny@unibas.ch"
__license__ = "GPL"

# imports
import os
import sys
import csv
import time
import HTSeq
import pandas as pd
import numpy as np
import seaborn as sns
import pylab as pl
from Bio import SeqIO
from collections import OrderedDict, defaultdict
from argparse import ArgumentParser, RawTextHelpFormatter
#
# add path to sys and import snoRNA module
#
file_dir = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.join(file_dir, "/import/bc2/home/zavolan/gumiennr/Pipelines/Pipelines/pipeline_snoRNASearch_snOPY/modules"))
from snoRNA import CD_snoRNA
#
# Params for figures
#
pl.rcParams['figure.figsize'] = (14, 10)
pl.rcParams['ytick.labelsize'] = 20
pl.rcParams['xtick.labelsize'] = 20
pl.rcParams['axes.labelsize'] = 23
pl.rcParams['legend.fontsize'] = 20
sns.set_style('ticks')
c1, c2, c3 = sns.color_palette("Set1", 3)

parser = ArgumentParser(description=__doc__, formatter_class=RawTextHelpFormatter)
parser.add_argument("-v",
                    "--verbose",
                    dest="verbose",
                    action="store_true",
                    default=False,
                    help="Be loud!")
parser.add_argument("--snoRNAs",
                    dest="snoRNAs",
                    required=True,
                    help="Table with snoRNAs")
parser.add_argument("--input",
                    dest="input",
                    required=True,
                    help="Input file in tab format.")
parser.add_argument("--output",
                    dest="output",
                    help="Output file in tab format.")
parser.add_argument("--type",
                    dest="type",
                    required=True,
                    choices=("CD", "HACA"),
                    help="Type of snoRNA")
parser.add_argument("--window",
                    dest="window",
                    type=int,
                    default=100,
                    help="Window, defaults to 100")
parser.add_argument("--smooth-window",
                    dest="smooth_window",
                    type=int,
                    default=1,
                    help="Smoothing window length, defaults to 1")
parser.add_argument("--dir",
                    dest="dir",
                    default="Plots",
                    help="Direcory for plots, defaults to Plots")

try:
    options = parser.parse_args()
except Exception, e:
    parser.print_help()

# redefine a functions for writing to stdout and stderr to save some writting
syserr = sys.stderr.write
sysout = sys.stdout.write

class Signal:
    def __init__(self, chrom, pos, strand):
        self.chrom = chrom
        self.pos = pos
        self.strand = strand


def main():
    """Main logic of the script"""
    if options.verbose:
        syserr("Reading snoRNAs\n")
    snoRNAs = read_snoRNAs_to_dict(options.snoRNAs, options.type)
    signals = {snorid: Signal(snorid, snor.d_box[0], "+") for snorid, snor in snoRNAs.iteritems()}
    names_string = "read_id\tsnor_id\tscore\tbeg\tend\twhich_end\ttarget\tqbeg\tqend\tisambig\t"
    names_string += "matches:mismatches:identity\tshuf_snor\tshuf_score\tshuf_target_len\tshuf_qbeg\tshuf_qend"
    names = names_string.split("\t")
    my_table = pd.read_table(options.input, header=None, names=names)
    my_table["id"] = [i.split(":")[0] for i in my_table["read_id"]]
    scores = my_table.dropna(subset=["score"])
    shuf_scores = my_table.dropna(subset=["shuf_score"])
    #
    # Make coverage for each snoRNA
    #
    if options.verbose:
        syserr("Making coverage\n")
    coverage = HTSeq.GenomicArray("auto", stranded=True, typecode="d")
    shuf_coverage = HTSeq.GenomicArray("auto", stranded=True, typecode="d")
    #
    # coverage for actual read
    #
    i = 0
    scoreslen = len(scores.index)
    for ind, row in scores.iterrows():
        syserr(" - normal: %i out of %i\r" % (i, scoreslen))
        iv = HTSeq.GenomicInterval(chrom=row.snor_id,
                                   start=int(row.qbeg),
                                   end=int(row.qend),
                                   strand="+")
        coverage[iv] += 1
        i += 1
    syserr("\n")
    #
    # coverage for shuffled read
    #
    i = 0
    shufscoreslen = len(shuf_scores.index)
    for ind, row in shuf_scores.iterrows():
        syserr(" - shuffled: %i out of %i\r" % (i, shufscoreslen))
        iv = HTSeq.GenomicInterval(chrom=row.shuf_snor,
                                   start=int(row.shuf_qbeg),
                                   end=int(row.shuf_qend),
                                   strand="+")
        shuf_coverage[iv] += 1
        i += 1
    syserr("\n")
    #
    # Make profile
    #
    if options.verbose:
        syserr("Making profile\n")
    profile = make_profile(signals, coverage, options.window)
    profile_count = get_aggregated_profile(profile, "sum")
    profile = get_normalized_profile(profile, len(my_table))
    shuf_profile = make_profile(signals, shuf_coverage, options.window)
    shuf_profile_count = get_aggregated_profile(shuf_profile, "sum")
    shuf_profile = get_normalized_profile(shuf_profile, len(my_table))
    #
    # Plot normalized
    #
    df = pd.DataFrame([np.arange(-options.window, options.window + 1),
                       profile,
                       shuf_profile,
                       profile_count,
                       shuf_profile_count]).transpose()
    df.index = df[0]
    df.to_csv(options.output, sep='\t', header=False, index=False)
    ndf = pd.rolling_mean(df[[1, 2, 3, 4]], options.smooth_window, center=True)
    ndf = ndf.dropna()
    ndf.index.name = "position"
    ndf.columns = ["normed_count", "normed_count_shuffled", "count", "count_shuffled"]

    if options.verbose:
        syserr("Plotting normalized profile\n")
    fig, ax = pl.subplots()
    sns.set_style("ticks")
    ax.plot(ndf.index, ndf["normed_count"].values, c=c1, lw=4, label="Reads")
    ax.plot(ndf.index, ndf["normed_count_shuffled"].values, c=c2, lw=4, label="Shuffled Reads")
    sns.despine()
    ax.set_xlim((-options.window, options.window))
    ax.set_ylabel("Normalized count", fontsize=18)
    ax.set_xlabel("Position around signal", fontsize=18)
    pl.legend()
    fig.tight_layout()
    fig.savefig(os.path.join(options.dir, "profile_around_d_box_normed.pdf"))
    fig.savefig(os.path.join(options.dir, "profile_around_d_box_normed.png"), dpi=300)
    #
    # plot count
    #
    if options.verbose:
        syserr("Plotting count profile\n")
    fig, ax = pl.subplots()
    sns.set_style("ticks")
    ax.plot(ndf.index, ndf["count"].values, c=c1, lw=4, label="Reads")
    ax.plot(ndf.index, ndf["count_shuffled"].values, c=c2, lw=4, label="Shuffled Reads")
    sns.despine()
    ax.set_xlim((-options.window, options.window))
    ax.set_ylabel("Count", fontsize=18)
    ax.set_xlabel("Position around signal", fontsize=18)
    pl.legend()
    fig.tight_layout()
    fig.savefig(os.path.join(options.dir, "profile_around_d_box.pdf"))
    fig.savefig(os.path.join(options.dir, "profile_around_d_box.png"), dpi=300)

def make_profile(signals, coverage, window_length, stranded=True):
    # profile = np.zeros(2*window_length, dtype='i')
    profile = {}
    for signalname, signal in signals.iteritems():
        try:
            sig_start = signal.pos - window_length
            if sig_start < 0:
                if stranded:
                    window = HTSeq.GenomicInterval(signal.chrom,
                                                   0,
                                                   signal.pos + window_length + 1,
                                                   signal.strand)
                else:
                    window = HTSeq.GenomicInterval(signal.chrom,
                                                   0,
                                                   signal.pos + window_length + 1,
                                                   ".")
                mywin = list(np.fromiter(coverage[window], dtype="i"))
                wincvg = np.fromiter([0]*abs(sig_start) + mywin,
                                     dtype='i',
                                     count=2*window_length + 1)
            else:
                if stranded:
                    window = HTSeq.GenomicInterval(signal.chrom,
                                                   signal.pos - window_length,
                                                   signal.pos + window_length + 1,
                                                   signal.strand)
                else:
                    window = HTSeq.GenomicInterval(signal.chrom,
                                                   signal.pos - window_length,
                                                   signal.pos + window_length + 1,
                                                   ".")
                wincvg = np.fromiter(coverage[window], dtype='i', count=2*window_length + 1)
            if signal.strand == "+":
                # profile += wincvg
                # profile.append(wincvg)
                profile[signalname] = wincvg
            else:
                # profile += wincvg[::-1]
                # profile.append(wincvg[::-1])
                profile[signalname] = wincvg[::-1]
        except Exception, e:
            print e
            import ipdb; ipdb.set_trace()
    return profile


def get_aggregated_profile(profile, type_of_aggregation):
    aggregated_profile = []
    for name, prof in profile.iteritems():
        aggregated_profile.append(prof)
    if type_of_aggregation == 'mean' or type_of_aggregation == 'average':
        return np.asarray(aggregated_profile).mean(axis=0)
    elif type_of_aggregation == 'sum':
        return np.asarray(aggregated_profile).sum(axis=0)


def get_normalized_profile(profile, library_size):
    profile_normed = []
    for name, prof in profile.iteritems():
        if prof.sum() == 0.0:
            profile_normed.append(prof)
        else:
            profile_normed.append(prof/float(prof.sum()))
            # profile_normed.append((((prof/float(prof.sum())))/library_size)*1000000.0)
    profile_normed = np.asarray(profile_normed).mean(axis=0)
    return profile_normed

def read_fasta_to_dict(path_to_file):
    """Read fasta file into dictionary

    Args:
        path_to_file (str): path to file

    Returns: dictionary

    """
    if options.verbose:
        syserr("Reading sequences from %s \n" % (path_to_file))
    try:
        seq_obj = open(path_to_file, 'Ur')
        seqs = OrderedDict()
        for seq in SeqIO.parse(seq_obj, 'fasta'):
            read = str(seq.seq).upper().replace("U", "T")
            seqs[str(seq.id)] = read
    except IOError:
        raise IOError('Cannot read from %s' % (path_to_file))

    return seqs


def read_snoRNAs_to_dict(path, type_of_snor):
    names = ["chrom",
             "start",
             "end",
             "snor_id",
             "mod_type",
             "strand",
             "sequence",
             "box_d",
             "box_c",
             "box_h",
             "box_aca",
             "alias",
             "gene_name",
             "accession",
             "mod_site",
             "host_gene",
             "host_id",
             "organization",
             "organism",
             "note"]
    if type_of_snor == "CD":
        snoRNAs = pd.read_table(path, names=names)
        if len(set(snoRNAs.mod_type)) != 1:
            Exception("More than one type of snoRNAs detected: %s" % str(set(snoRNAs.mod_type)))
        counter = 0
        snor_dict = {}
        for ind, snor in snoRNAs.iterrows():
            try:
                s = CD_snoRNA(snor_id=snor.snor_id,
                              organism=snor.organism,
                              chrom=snor.chrom,
                              start=snor.start,
                              end=snor.end,
                              strand=snor.strand,
                              sequence=snor.sequence,
                              snor_type=snor.mod_type,
                              d_boxes=snor.box_d,
                              c_boxes=snor.box_c,
                              switch_boxes=True,
                              alias=snor.alias,
                              gene_name=snor.gene_name,
                              accession=snor.accession,
                              modified_sites=snor.mod_site,
                              host_id=snor.host_id,
                              organization=snor.organization,
                              note=snor.note)
                counter += 1
                snor_dict[s.snor_id] = s
            except Exception, e:
                sys.stderr.write(str(e) + "\n")
        return snor_dict

    elif options.type == "HACA":
        raise NotImplementedError("Pipeline was not implemented for HACA box snoRNAs yet")

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
