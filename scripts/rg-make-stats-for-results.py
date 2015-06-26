#!/usr/bin/env python
"""
Make some useful plots for results
"""

__date__ = "2015-02-25"
__author__ = "Rafal Gumienny"
__email__ = "r.gumienny@unibas.ch"
__license__ = "GPL"

# imports
import sys
import re
import os
import csv
import glob
import time
import HTSeq
import pandas as pd
import numpy as np
import pylab as pl
import seaborn as sns
from random import *
from collections import OrderedDict
from Bio import SeqIO
from argparse import ArgumentParser, RawTextHelpFormatter
# add path to sys and import snoRNA module
file_dir = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.join(file_dir, "../modules"))
from snoRNA import CD_snoRNA

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
parser.add_argument("--results-probability-complex",
                    dest="results_probability_complex",
                    required=True,
                    help="Main part of the results")
# parser.add_argument("--results-probability-simple",
#                     dest="results_probability_simple",
#                     required=True,
#                     help="Results obtained from simple model of interaction")
parser.add_argument("--results-raw",
                    dest="results_raw",
                    required=True,
                    help="Row results")
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
parser.add_argument("--genome-dir",
                    dest="genome_dir",
                    required=True,
                    help="Path to genome directory where the chromosomes are stored")


try:
    options = parser.parse_args()
except Exception, e:
    parser.print_help()

# redefine a functions for writing to stdout and stderr to save some writting
syserr = sys.stderr.write
sysout = sys.stdout.write

class Point:
    def __init__(self, x, y):
        self.x = x
        self.y = y


def main():
    """Main logic of the script"""
    if options.verbose:
        syserr("Reading snoRNAs\n")
    snoRNAs = read_snoRNAs_to_dict(options.snoRNAs, options.type)
    chromosome_mapping = get_chomosomes_mapping(snoRNAs)

    if options.verbose:
        syserr("Reading models\n")
    # simple_model = pd.read_table(options.results_probability_simple)
    complex_model = pd.read_table(options.results_probability_complex)

    # simple_model_stats = get_predictor_paramters(simple_model, chromosome_mapping, snoRNAs)
    # sfig, sax = plot_predictor_parameters(simple_model_stats, "Simple model")
    # sfig.savefig(os.path.join(options.dir, "predictor_simple_model.pdf"))
    # sfig.savefig(os.path.join(options.dir, "predictor_simple_model.png"),
                 # dpi=300)
    complex_model_stats = get_predictor_paramters(complex_model, chromosome_mapping, snoRNAs)
    cfig, cax = plot_predictor_parameters(complex_model_stats, "Complex model")
    cfig.savefig(os.path.join(options.dir, "predictor_complex_model.pdf"))
    cfig.savefig(os.path.join(options.dir, "predictor_complex_model.png"),
                 dpi=300)

    if options.verbose:
        syserr("Reading chromosomes\n")
    chromosome_lengths = get_chromosome_lengths(options.genome_dir, "*.fa")
    chromosome_lengths = {k.split("|")[0]: v for k, v in chromosome_lengths.iteritems()}

    if options.verbose:
        syserr("Calculating coverage\n")
    coverage, chromosomes = get_coverage(options.results_raw, chromosomes=chromosome_lengths)

    if options.verbose:
        syserr("Plotting probabilities\n")
    for mod_chrom_alt, mod_chrom in chromosome_mapping.iteritems():
        # if mod_chrom not in ["RNA18S"]:
            # continue
        # if options.verbose:
        #     syserr(" - %s (simple model)\n" % mod_chrom)
        # try:
        #     sim_fig, sim_ax = plot_chomosome(mod_chrom,
        #                                      chromosome_lengths,
        #                                      simple_model,
        #                                      get_target_sites(snoRNAs, mod_chrom_alt),
        #                                      0.0)
        #     sns.despine()
        #     sim_fig.savefig(os.path.join(options.dir,
        #                                  "probabilities_simple_model_%s.pdf" % mod_chrom))
        #     sim_fig.savefig(os.path.join(options.dir,
        #                                  "probabilities_simple_model_%s.png" % mod_chrom),
        #                     dpi=300)
        #     sim_fig, sim_ax = plot_chomosome(mod_chrom,
        #                                      chromosome_lengths,
        #                                      simple_model,
        #                                      get_target_sites(snoRNAs, mod_chrom_alt),
        #                                      0.5)
        #     sns.despine()
        #     sim_fig.savefig(os.path.join(options.dir,
        #                                  "probabilities_simple_model_%s_05.pdf" % mod_chrom))
        #     sim_fig.savefig(os.path.join(options.dir,
        #                                  "probabilities_simple_model_%s_05.png" % mod_chrom),
        #                     dpi=300)
        # except Exception, e:
        #     syserr("Cannot plot simple model for %s: %s\n" % (mod_chrom, str(e)))

        if options.verbose:
            syserr(" - %s (complex model)\n" % mod_chrom)
        try:
            com_fig, com_ax = plot_chomosome(mod_chrom,
                                             chromosome_lengths,
                                             complex_model,
                                             get_target_sites(snoRNAs, mod_chrom_alt),
                                             0.0)
            sns.despine()
            com_fig.savefig(os.path.join(options.dir,
                                         "probabilities_complex_model_%s.pdf" % mod_chrom))
            com_fig.savefig(os.path.join(options.dir,
                                         "probabilities_complex_model_%s.png" % mod_chrom),
                            dpi=300)
            com_fig, com_ax = plot_chomosome(mod_chrom,
                                             chromosome_lengths,
                                             complex_model,
                                             get_target_sites(snoRNAs, mod_chrom_alt),
                                             0.5)
            sns.despine()
            com_fig.savefig(os.path.join(options.dir,
                                         "probabilities_complex_model_%s_05.pdf" % mod_chrom))
            com_fig.savefig(os.path.join(options.dir,
                                         "probabilities_complex_model_%s_05.png" % mod_chrom),
                            dpi=300)
        except Exception, e:
            syserr("Cannot plot complex model for %s: %s\n" % (mod_chrom, str(e)))


    # TODO ideal situation is when all the chromosomes matches
    if options.verbose:
        syserr("Plotting coverage\n")
    for mod_chrom_alt, mod_chrom in chromosome_mapping.iteritems():
        if options.verbose:
            syserr(" - %s\n" % mod_chrom)
        fig, ax = pl.subplots()
        profile = get_profile(coverage,
                              mod_chrom,
                              0,
                              chromosome_lengths[mod_chrom],
                              '+')
        ax.plot(profile, color=c2, label='Hybrid Coverage')
        ax.plot(get_target_sites(snoRNAs, mod_chrom_alt),
                [[-1*max(profile)/25.0]]*len(get_target_sites(snoRNAs, mod_chrom_alt)),
                linestyle='none',
                marker='o',
                markeredgecolor='none',
                color='ForestGreen',
                markersize=6,
                label='Methylation')
        ax.set_ylim(bottom=int(np.floor(-1*max(profile)/13.0)))
        ax.set_xlabel("Position")
        ax.set_title(mod_chrom)
        pl.legend()
        sns.despine()
        fig.savefig(os.path.join(options.dir, "raw_read_profile_%s.pdf" % mod_chrom))
        fig.savefig(os.path.join(options.dir, "raw_read_profile_%s.png" % mod_chrom), dpi=300)


def get_predictor_paramters(data, chroms, snoRNAs):
    model_data = {}
    for threshold in np.arange(0.0, 1.0, 0.01):
        ndf = data[(data["Probability"] < threshold)]
        df = data[(data["Probability"] >= threshold)]

        TP = 0.0
        FP = 0.0
        TN = 0.0
        FN = 0.0
        for mod_chrom_alt, mod_chrom in chroms.iteritems():
            # Calculate statistics
            if mod_chrom not in ["RNA28S", "RNA18S"]:
                continue
            mods = get_target_sites(snoRNAs, mod_chrom_alt)
            lTP, lTN, lFP, lFN, known_sites = get_statistics(mod_chrom, df, ndf, mods)
            TP += lTP
            FP += lFP
            TN += lTN
            FN += lFN

        try:
            TPR = TP / (TP + FN)  # sensitivity
        except ZeroDivisionError:
            TPR = 0.0
        try:
            TNR = TN / (FP + TN)  # specificity
        except ZeroDivisionError:
            TNR = 0.0
        try:
            PPV = TP / (TP + FP)  # positive predictive value - precision
        except ZeroDivisionError:
            PPV = 0.0
        try:
            NPV = TN / (TN + FN)  # negative predictive value
        except ZeroDivisionError:
            NPV = 0.0
        try:
            FPR = FP / (FP + TN)  # fall-out
        except ZeroDivisionError:
            FPR = 0.0
        try:
            FNR = FN / (FN + TP)  # false negative rate
        except ZeroDivisionError:
            FNR = 0.0
        try:
            FDR = FP / (TP + FP)  # false discovery rate
        except ZeroDivisionError:
            FDR = 0.0
        try:
            accuracy = (TP + TN) / (TP + FP + FN + TN)  # accuracy
        except ZeroDivisionError:
            accuracy = 0.0

        model_data[threshold] = {"TPR (Sensitivity)": TPR,
                                 "TNR (Specificity)": TNR,
                                 "PPV (Precision)": PPV,
                                 "NPV": NPV,
                                 "FPR (Fall-out)": FPR,
                                 "FNR": FNR,
                                 "TotalNumberOfFoundSites": (TP + FP),
                                 "Accuracy": accuracy,
                                 "FDR": FDR
                                 }
    return pd.DataFrame(model_data).transpose()


def get_statistics(chrom, data, neg_data, mod):
    known_sites = set(mod)
    found_sites = set(data[data.chrom == chrom]["Modification"].tolist())
    not_found_sites = set(neg_data[neg_data.chrom == chrom]["Modification"].tolist())
    true_positives = float(len(known_sites & found_sites))
    true_negatives = float(len(not_found_sites - known_sites))
    false_positives = float(len(found_sites - known_sites))
    false_negatives = float(len(not_found_sites & known_sites))
    return (true_positives,
            true_negatives,
            false_positives,
            false_negatives,
            len(known_sites))


def plot_predictor_parameters(data_to_plot, title=""):
    columns_to_plot = ['TPR (Sensitivity)',
                       'TNR (Specificity)',
                       'PPV (Precision)',
                       'FPR (Fall-out)',
                       'Accuracy']
    fig, ax = pl.subplots()
    data_to_plot[columns_to_plot].plot(lw=5, ax=ax)
    ax.set_xlim((0.0, 1.4))
    ax.set_ylim((0.0, 1.0))
    ax.set_xlabel('Threshold', fontsize=30)
    ax.set_ylabel('Value', fontsize=30)
    ax.set_title(title, fontsize=25)
    sns.despine()
    legend = pl.legend(fontsize=20)
    pl.tight_layout()
    return fig, ax


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


def get_coverage(path_to_file, chromosomes="auto", stranded=True, storage="step"):
    coverage = HTSeq.GenomicArray(chroms=chromosomes,
                                  stranded=stranded,
                                  typecode="d",
                                  storage=storage)
    # print "I am in get_coverage"
    with open(path_to_file, 'r') as f:
        chromosomes = set()
        for row in csv.reader(f, delimiter='\t'):
            chromosomes.add(row[0])
            coverage[HTSeq.GenomicInterval(chrom=row[0],
                                           start=int(row[1]),
                                           end=int(row[2]),
                                           strand=row[5])] += 1.0
    return coverage, chromosomes


def get_chromosome_lengths(path_to_chromosome_directory, glob_regular_expression):
    chromosome_lengths = {}
    for chrom in glob.glob(os.path.join(path_to_chromosome_directory,
                                        glob_regular_expression)):
        with open(chrom) as f:
            for rec in SeqIO.parse(f, 'fasta'):
                # print " - %s\n" % rec.id
                chromosome_lengths[str(rec.id)] = len(rec.seq)
    return chromosome_lengths


def get_profile(coverage, chrom, start, end, strand):
        window = HTSeq.GenomicInterval(chrom,
                                       start,
                                       end,
                                       strand)
        wincvg = np.fromiter(coverage[window],
                             dtype='i',
                             count=end-start)
        return wincvg


def plot_chomosome(chrom, chrom_lengths, data, mod_sites, thresh=0.7):
    fig, ax = pl.subplots()
    barcolor = 'RoyalBlue'
    metcolor = "FireBrick"
    xdata_prob = data[(data.chrom == chrom) & (data.Probability >= thresh)]["Modification"].tolist()
    ydata_prob = data[(data.chrom == chrom) & (data.Probability >= thresh)]["Probability"].tolist()
    ydata_true = [-0.1]*len(mod_sites)
    ax.bar(xdata_prob,
           ydata_prob,
           color=barcolor,
           edgecolor=barcolor,
           width=3,
           label='Found Modification Site')
    ax.bar(mod_sites,
           ydata_true,
           color=metcolor,
           edgecolor=metcolor,
           width=3,
           label='True Methylation Site')
    ax.axhline(0, color='black', ls='solid', lw=1)
    pl.legend()
    pl.ylim(-0.2, 1.2)
    pl.xlim(left=0, right=chrom_lengths[chrom])
    txt_height = 0.01*(pl.ylim()[1] - pl.ylim()[0])
    txt_width = 0.001*(pl.xlim()[1] - pl.xlim()[0])
    for x, y in zip(mod_sites, ydata_true):
        ax.text(x - txt_width, y*1.2, '%d' % int(x), rotation=90, fontsize=6)
    for x, y in zip(xdata_prob, ydata_prob):
        ax.text(x - txt_width, y*1.05, '%d' % int(x), rotation=90, fontsize=6)
    return fig, ax


def get_target_sites(snors, chrom):
    sites = []
    for snor_id, snor in snors.iteritems():
        if snor.modified_sites:
            # print snor.modified_sites, chrom
            if chrom in snor.modified_sites:
                for site, nuc in snor.modified_sites[chrom]:
                    sites.append(site)
    return sites


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
