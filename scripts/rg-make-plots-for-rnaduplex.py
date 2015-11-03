#!/usr/bin/env python
"""
Make some useful plots for RNAduplex results
"""

__date__ = "2015-10-28"
__author__ = "Rafal Gumienny"
__email__ = "r.gumienny@unibas.ch"
__license__ = "GPL"

# imports
import matplotlib
matplotlib.use('Agg')
import sys
import time
from contextlib import contextmanager
from argparse import ArgumentParser, RawTextHelpFormatter


import os
import HTSeq
import pylab as pl
import pandas as pd
import seaborn as sns
import numpy as np
import HTSeq
from collections import Counter
from Bio import SeqIO

file_dir = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.join(file_dir, "../modules"))
sys.path.append("/scicore/home/zavolan/gumiennr/PythonModules/MetaProfile/")
from snoRNA import read_snoRNAs_from_table
import MetaProfile

pl.rcParams['figure.figsize'] = (14, 10)
pl.rcParams['ytick.labelsize'] = 20
pl.rcParams['xtick.labelsize'] = 20
pl.rcParams['axes.labelsize'] = 23
pl.rcParams['legend.fontsize'] = 20
sns.set_style('ticks')
colors = ["windows blue", "amber", "green", "tomato", "dusty purple"]
c1, c2, c3, c4, c5 = sns.xkcd_palette(colors)


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
                    help="Input file in TAB")
parser.add_argument("--snoRNAs",
                    dest="snoRNAs",
                    required=True,
                    help="Table with snoRNAs")
parser.add_argument("--type",
                    dest="type",
                    required=True,
                    choices=("CD", "HACA"),
                    help="Type of snoRNA")
parser.add_argument("--dir",
                    dest="dir",
                    default="Plots",
                    help="Directory to store plots, defaults to Plots")
parser.add_argument("--threshold",
                    dest="threshold",
                    type=float,
                    default=-25.0,
                    help="Threshold for RNAduplex energy, defaults to -25.0")


try:
    options = parser.parse_args()
except Exception, e:
    parser.print_help()

# redefine a functions for writing to stdout and stderr to save some writting
syserr = sys.stderr.write
sysout = sys.stdout.write

class Signal(MetaProfile.Signal):
    def __init__(self, name, coverage):
        self.name = name
        self.coverage = coverage
        self.stranded = True


def main():
    """Main logic of the script"""
    snoRNAs = read_snoRNAs_from_table(options.snoRNAs, options.type, True)
    mod_sites = get_target_sites(snoRNAs, get_chomosomes_mapping(snoRNAs))
    real_mod_sites = get_real_target_sites(snoRNAs, get_chomosomes_mapping(snoRNAs))

    df = pd.read_table(options.input, header=None)
    df = df[~df[21].map(str).str.contains("repeat")]
    df['en_normed'] = df[13] / df[14]
    df["is_known"] = [is_known(chrom, snor_id, pos, mod_sites)
                            for chrom, snor_id, pos in zip(df[0], df[3], df[10])]
    df["is_orphan"] = [snoRNAs[snor_id].is_orphan() for snor_id in df[3]]
    cds_sites = df[(df[14] < options.threshold) & (df[25].map(str) == 'CDS')]
    UTR3_sites = df[(df[14] < options.threshold) & (df[25].map(str) == 'three_prime_UTR')]
    intron_sites = df[(df[14] < options.threshold) & (df[24].map(str) == 'intron')]
    non_canonical_sites = df[(df[14] < options.threshold) & (df[6] > -7)]
    canonical_sites = df[(df[6] < -15.0)]


    #
    # Duplex energy vs random
    #
    fig, ax = pl.subplots()
    df[14].hist(bins=np.arange(-100, 0, 0.1), normed=True, histtype='step', lw=2, cumulative=True, ax=ax, label='RNAduplex energy - normal')
    df[15].hist(bins=np.arange(-100, 0, 0.1), normed=True, histtype='step', lw=2, cumulative=True, ax=ax, label='RNAduplex energy - random snoRNA')
    df[16].hist(bins=np.arange(-100, 0, 0.1), normed=True, histtype='step', lw=2, cumulative=True, ax=ax, label='RNAduplex energy - shuffled target')
    pl.xlim(-50, -5)
    pl.ylabel("Empirical CDF")
    pl.xlabel('RNAduplex Energy')
    sns.despine()
    pl.legend(loc='upper left', frameon=True, fancybox=True)
    fig.savefig(os.path.join(options.dir, "RNAduplex_energy_cdf.pdf"))
    fig.savefig(os.path.join(options.dir, "RNAduplex_energy_cdf.png"), dpi=300)
    pl.clf()

    #
    # Duplex energy vs GC content
    #
    g = sns.jointplot(x=13, y=14, kind='reg', data=df, size=10)
    pl.xlabel("GC content")
    pl.ylabel("Duplex energy")
    pl.savefig(os.path.join(options.dir, "RNAduplex_energy_vs_gc_content.pdf"))
    pl.savefig(os.path.join(options.dir, "RNAduplex_energy_vs_gc_content.png"), dpi=300)
    pl.clf()


    #
    # Duplex structures for each snoRNA
    #
    if options.verbose:
        syserr("Plotting duplex structures:\n")
    plot_duplex_structures(df, snoRNAs, 'all')
    pl.clf()
    plot_duplex_structures(cds_sites, snoRNAs, 'CDS')
    pl.clf()
    plot_duplex_structures(intron_sites, snoRNAs, 'intron')
    pl.clf()
    plot_duplex_structures(UTR3_sites, snoRNAs, '3UTR')
    pl.clf()
    plot_duplex_structures(non_canonical_sites, snoRNAs, 'noncanonical_sites')
    pl.clf()
    plot_duplex_structures(canonical_sites, snoRNAs, 'canonical_sites')
    pl.clf()

    #
    # Coverage of paired nucleotides
    #
    coverage = HTSeq.GenomicArray(chroms='auto', typecode='d')
    coverage_random = HTSeq.GenomicArray(chroms='auto', typecode='d')
    for idx, row in df.iterrows():
        for iv in get_positions_from_structure(row[3], row[17]):
            coverage[iv] += 1
        for iv in get_positions_from_structure(row[3], row[20]):
            coverage_random[iv] += 1
    signals = [Signal('structure', coverage), Signal('structure_random', coverage_random)]
    with open(os.path.join(options.dir, 'snornas_dbox_pos.bed'), 'w') as out:
        for snor_name, snor in snoRNAs.iteritems():
            out.write("%s\t%i\t%i\t%s\t%s\t+\n" % (snor_name, snor.d_box[0] - 1, snor.d_box[0], snor_name, '1'))
    pc = MetaProfile.MetaProfiler(signals=signals, windows=[MetaProfile.Window('d_box',
                                                                               os.path.join(options.dir, 'snornas_dbox_pos.bed'),
                                                                               pseudocount=0,
                                                                               window_length=40)])
    pc.create_profiles()
    fig, ax = pl.subplots()
    for prof_name, col in zip(['structure on d_box', 'structure_random on d_box'], [c1, c2]):
        pc.profiles[prof_name].plot_line('normalized_to_gene', ax=ax, color=col)
    pl.legend(ax.lines, ['Structure on snoRNA', 'Structure on random snoRNA'],
                        frameon=True, fancybox=True)
    sns.despine()
    pl.xlabel("Position around D-box")
    pl.ylabel("Normalized count")
    fig.savefig(os.path.join(options.dir, "Dbox_structure_coverage.pdf"))
    fig.savefig(os.path.join(options.dir, "Dbox_structure_coverage.png"), dpi=300)


def plot_duplex_structures(df, snoRNAs, suffix=''):
    for snor in set(df[3]):
        if options.verbose:
            syserr(" - %s\n" % snor)
        hybrids = np.asarray([list(i) for i in df[df[3]==snor][17]])
        plot_hybrid_statistics_with_boxes(snor, snoRNAs, hybrids)
        pl.savefig(os.path.join(options.dir, "Structures/structure_%s_%s_%s.pdf" % (snor, str(snoRNAs[snor].snor_name), suffix)))
        pl.savefig(os.path.join(options.dir, "Structures/structure_%s_%s_%s.png" % (snor, str(snoRNAs[snor].snor_name), suffix)), dpi=300)
        fig = pl.gcf()
        pl.close(fig)


def is_known(chrom, snor_id, pos, modsites):
    if pd.isnull(pos):
        return False
    else:
        return "%s:%s:%s" % (chrom, snor_id, int(pos)) in modsites


def plot_hybrid_statistics_with_boxes(snor, snoRNAs, hybrids):
    ax = plot_hybrid_statistics(hybrids, title=snor + ", " + str(snoRNAs[snor].snor_name) + " (" + str(len(hybrids)) + ")")
    ax.get_xaxis().set_visible(False)
    line_dbox = pl.Line2D([snoRNAs[snor].d_box[0] - 0.3 - 1,
                           snoRNAs[snor].d_box[1] + 0.7 - 1],
                          [-0.02, -0.02],
                        lw=12., color=c3, solid_capstyle='butt')
    line_dbox.set_clip_on(False)

    line_modification_dbox = pl.Line2D([snoRNAs[snor].d_box[0] - 0.2 - 6,
                                        snoRNAs[snor].d_box[0] + 0.6 - 6],
                                       [-0.02, -0.02],
                        lw=12., color='b', solid_capstyle='butt')
    line_modification_dbox.set_clip_on(False)

    try:
        line_cbox = pl.Line2D([snoRNAs[snor].c_box[0] - 0.3 - 1,
                               snoRNAs[snor].c_box[1] + 0.7 - 1],
                              [-0.02, -0.02],
                        lw=12., color=c4, solid_capstyle='butt')
        line_cbox.set_clip_on(False)
    except TypeError:
        line_cbox = None

    try:
        line_dpbox = pl.Line2D([snoRNAs[snor].dprime_box[0] - 0.3 - 1,
                                snoRNAs[snor].dprime_box[1] + 0.7 - 1],
                               [-0.02, -0.02],
                        lw=12., color=c3, alpha=0.5, solid_capstyle='butt')
        line_dpbox.set_clip_on(False)

        line_modification_dpbox = pl.Line2D([snoRNAs[snor].dprime_box[0] - 0.2 - 6,
                                             snoRNAs[snor].dprime_box[0] + 0.6 - 6],
                                            [-0.02, -0.02],
                        lw=12., color='b', alpha=0.5, solid_capstyle='butt')
        line_modification_dpbox.set_clip_on(False)
    except TypeError:
        line_dpbox = None
    try:
        line_cpbox = pl.Line2D([snoRNAs[snor].cprime_box[0] - 0.3 - 1, snoRNAs[snor].cprime_box[1] + 0.7 - 1], [-0.02, -0.02],
                        lw=12., color=c4, alpha=0.5, solid_capstyle='butt')
        line_cpbox.set_clip_on(False)
    except TypeError:
        line_cpbox = None

    ax.add_line(line_dbox)
    ax.add_line(line_modification_dbox)
    if line_cbox is not None:
        ax.add_line(line_cbox)
    if line_dpbox is not None:
        ax.add_line(line_dpbox)
        ax.add_line(line_modification_dpbox)
    if line_cpbox is not None:
        ax.add_line(line_cpbox)
    return ax


def plot_hybrid_statistics(hybrids, title=""):
    with sns.axes_style("white"):
        cl_dict = {}
        for i in range(len(hybrids[0])):
            cl_dict[i+1] = Counter(hybrids[:,i])
        my_colors = ["DarkOrange",
                     "DodgerBlue",
                     "Red",
                     "green"]
        cl_data = pd.DataFrame(cl_dict).transpose()
        cl_data = cl_data.fillna(0)
        cl_sum = cl_data.sum(axis=1)[1]
        cl_data = cl_data/cl_sum
        cl_data = cl_data.ix[:,["(", ")", ".", "#"]]
        ax = cl_data.plot(kind="bar",
                     stacked=True,
                     grid=False,
                     rot=0,
                     color=my_colors)
        legend = pl.legend(frameon=True, shadow=True, loc="upper left")
        texts = ["Left Pair",
                 "Right Pair",
                 "Unpaired",
                 "Dangling End"]
        for lab, leg_text in zip(texts, legend.get_texts()):
            leg_text.set_text(lab)
        for container in ax.containers:
            pl.setp(container, width=0.92)
        for tick in ax.xaxis.get_ticklabels():
            tick.set_fontsize(18)
        ax.set_ylabel("Frequency", fontsize=25)
        ax.set_xlabel("Position", fontsize=25)
        sns.despine(bottom=False, offset=0)
        ax.set_title(title, fontsize=22)
        return ax


def get_positions_from_structure(chrom, struc, strand='+'):
    for i, pos in enumerate(struc):
        if pos == '(':
            yield HTSeq.GenomicInterval(chrom, i, i + 1, strand=strand)


def get_target_sites(snors, chrom_map):
    my_sites = []
    for snor_id, snor in snors.iteritems():
        if snor.modified_sites:
            # print snor.modified_sites, chrom
            for chrom, sites in snor.modified_sites.iteritems():
                for site in sites:
                    my_sites.append("%s:%s:%i" % (chrom_map[chrom], snor_id, site[0]))
    return my_sites


def get_real_target_sites(snors, chrom_map):
    my_sites = []
    for snor_id, snor in snors.iteritems():
        if snor.modified_sites:
            # print snor.modified_sites, chrom
            for chrom, sites in snor.modified_sites.iteritems():
                for site in sites:
                    my_sites.append("%s:%i" % (chrom_map[chrom], site[0]))
    return my_sites


def get_chomosomes_mapping(snor_dict):
    mod_chroms = {}
    for sn, s in snor_dict.iteritems():
        if s.modified_sites:
            for chrom in s.modified_sites.keys():
                if chrom in ['18S', '28S']:
                    mod_chroms[chrom] = "RNA" + chrom
                elif chrom == '5.8S':
                    mod_chroms[chrom] = "RNA5S"
                else:
                    mod_chroms[chrom] = chrom
    return mod_chroms


# this function is also defined in utils but I put it here to avoid
# unnecessary module import that might not be available everywhere as
# it is my own module
@contextmanager
def smart_open(filepath, mode='r'):
    """Open file intelligently depending on the source

    :param filepath: can be both path to file or sys.stdin or sys.stdout
    :param mode: mode can be read "r" or write "w". Defaults to "r"
    :yield: context manager for file handle

    """
    if mode == 'r':
        if filepath is not sys.stdin:
            fh = open(filepath, 'r')
        else:
            fh = filepath
        try:
            yield fh
        except IOError as e:
            if fh is not sys.stdin:
                fh.close()
            elif e.errno == errno.EPIPE:
                pass
        finally:
            if fh is not sys.stdin:
                fh.close()
    elif mode == 'w':
        if filepath is not sys.stdout:
            fh = open(filepath, 'w')
        else:
            fh = filepath
        try:
            yield fh
        finally:
            if fh is not sys.stdout:
                fh.close()
    else:
        raise NoSuchModeException("No mode %s for file" % mode)


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
