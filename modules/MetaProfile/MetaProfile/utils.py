import os
import sys
import numpy as np
from pandas import DataFrame
import matplotlib.gridspec as gridspec
import pylab as pl
import seaborn as sns
import pandas as pd
from Bio import SeqIO, Seq, motifs
from Bio.Alphabet import IUPAC
import HTSeq

def create_nucleotide_frequency_from_fasta(path_to_fasta, sequence_length=None, verbose=False):
    """Create a pandas.DataFrame with the frequency of each nucleotide
    on each position. Sequences are converted to uppercase DNA and saved
    as IUPAC.ambiguous_dna, however only the results for A, C, T and G are
    returned.

    :param path_to_fasta: path to the fasta file
    :sequence_length: assumed length of the input sequences. If not provided length of the first
                      sequence will be chosen. It will be check and any sequece of different 
                      length will be ignored
    :returns: 2-tuple: pandas.DataFrame, Bio.motifs.Motif

    """
    seqs = []
    number_of_ignored = 0
    seq_counter = 0
    if not sequence_length:
        for rec in SeqIO.parse(path_to_fasta, 'fasta'):
            sequence_length = len(rec.seq)
            break
        sys.stderr.write("Presumed sequence length was not provided. First encounterd sequence length" +
                         " will be used i.e. %i\n" % sequence_length)
    for rec in SeqIO.parse(path_to_fasta, 'fasta'):
        seq_counter += 1
        if len(rec.seq) != sequence_length:
            if verbose:
                sys.stderr.write("%s has wrong length (%i)." % (rec.id, len(rec.seq)) +
                                 " It will be ignored.\n")
            number_of_ignored += 1
        else:
            seqs.append(Seq.Seq(str(rec.seq).upper().replace("U", "T"), IUPAC.ambiguous_dna))
    motifs_obj = motifs.create(seqs, IUPAC.ambiguous_dna)
    frequency_df = DataFrame(motifs_obj.pwm)[["A", "C", "T", "G"]]
    sys.stderr.write("%i sequences out of %i was ignored because of the length issue.\n" % (number_of_ignored, seq_counter))
    return frequency_df, motifs_obj


def plot_line_profile(profile, x=None, ax=None, smooth_by=1, **kwargs):
    if ax is None:
        fig, ax = pl.subplots()
    if profile.ndim == 1:
        if x is None:
            x = range(-(len(profile)//2), len(profile)//2 + len(profile)%2)
        assert len(profile) == len(x), "x vector and profile length are different"
        smoothed = pd.rolling_mean(profile, window=smooth_by, center=True)
        ax.plot(x, smoothed, **kwargs)
        return ax
    elif profile.ndim == 2:
        # if x is None:
        #     x = range(-(profile.shape[1]//2), profile.shape[1]//2 + profile.shape[1]%2)
        # assert profile.shape[1] == len(x), "x vector and profile length are different"
        smoothed = pd.rolling_mean(profile, window=smooth_by, center=True, axis=1)
        sns.tsplot(smoothed, ax=ax, **kwargs)
        offset = profile.shape[1]//2
        locs, labels = pl.xticks()
        t = pl.xticks(locs, map(lambda y: "%g" % y, locs - offset))
        return ax
    else:
        raise Exception("Input profile has wrong dimentionality")


def clean_axis(ax):
    """
    Remove ticks and tick labels

    Parameters
    ----------
    ax : matplotlib.axes
        Axes to modify
    """
    ax.get_xaxis().set_ticks([])
    ax.get_yaxis().set_ticks([])


def plot_heatmap_profile(df, ax=None, xlabel='', ylabel='', title='', percentiles=(1, 99), sort_by=None, **imshow_kws):
    if ax is None:
        fig, ax = pl.subplots()
    vmin = np.percentile(df, percentiles[0])
    vmax = np.percentile(df, percentiles[1])

    # heatmap
    if sort_by is not None:
        assert df.shape[0] == sort_by.shape[0], "Sort values and profile has different shape"
        ind = sort_by.argsort()
    else:
        ind = np.arange(df.shape[0])
        assert df.shape[0] == len(ind)

    axi = ax.imshow(df[ind, :], interpolation='none', aspect='auto',
                    origin='lower', vmin=vmin, vmax=vmax, **imshow_kws)
    ax.yaxis.set_label_position('right')
    ax.set_ylabel(ylabel)
    ax.set_xlabel(xlabel)
    clean_axis(ax)

    ax.set_title(title, fontsize=23)
    ax.set_ylabel(ylabel)
    ax.set_xlabel(xlabel)
    ax.heatmap_axis = axi
    ax.vmin = vmin
    ax.vmax = vmax
    return ax

def create_grid(figsize=(8, 12), strip=False, height_ratios=(1, 4), width_ratios=(1, 4), subplot_params=None):
    if subplot_params is None:
        subplot_params = {}
    fig = pl.figure(figsize=figsize)
    gs = gridspec.GridSpec(
        len(height_ratios),
        len(width_ratios),
        height_ratios=height_ratios,
        width_ratios=width_ratios,
        wspace=0.1, hspace=0.1,
        **subplot_params)
    fig.array_axes = plt.subplot(gs[1, 0:])
    fig.line_axes = pl.subplot(gs[0, 0:], sharex=fig.array_axes)
    fig.gs = gs
    return fig

def infer_extension(path):
    """Infer extension of a file based on its path

    :param path: @todo
    :returns: @todo

    """
    filename, extension = os.path.splitext(path)
    if extension == '':
        raise NoExtensionException("File has no extension, please provide file type manually")
    else:
        return extension[1:]
