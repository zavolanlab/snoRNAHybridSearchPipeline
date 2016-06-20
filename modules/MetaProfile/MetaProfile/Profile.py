import os
import sys
import pylab as pl
from numpy import fromiter, asarray, nan_to_num
from .utils import plot_line_profile, plot_heatmap_profile

class WrongFiletypeException(Exception): pass
class NoSuchAFileException(Exception): pass
class NoSuchAProfileException(Exception): pass
class NoSuchAProfileTypeException(Exception): pass

class Profile(object):

    """Profile of the ..."""

    def __init__(self, signal, window):
        """Initialize profile object

        :param signals: @todo
        :param windows: @todo

        """
        #TODO assert for signal and windows type
        self.signal = signal
        self.window = window
        self.name = "%s on %s" % (signal.name, window.name)
        self.profile_raw = None
        self.profile_normalized_to_gene = None
        self.profile_normalized_to_library = None
        self.__create_profile()

    def __create_profile(self):
        """Prepare profile with given signal and windows. For each window a pseudocount
        is added.
        """
        sys.stderr.write("Calculating profile for %s on %s\n" % (self.signal.name, self.window.name))
        profile = []
        window_lengths = set()
        for window in self.window.windows:
            window_lengths.add(window.length)
            assert len(window_lengths) == 1, "Windows has at least two different lengths"
            try:
                wincvg = fromiter(self.signal.coverage[window], dtype='d', count=window.length)
            except IndexError:
                sys.stderr.write("Wrong window: %s\n" % str(window))
                continue
            if window.strand == "+" and self.signal.stranded:
                profile.append(wincvg + self.window.pseudocount)
            elif window.strand == "-" and self.signal.stranded:
                profile.append(wincvg[::-1] + self.window.pseudocount)
            else:
                profile.append(wincvg + self.window.pseudocount)

        self.profile_raw = asarray(profile)

    def normalize_to_library(self):
        """Normalize profile to the library size
        """
        self.profile_normalized_to_library = asarray(map(lambda i: i/float(self.signal.library_size), self.profile_raw))

    def normalize_to_gene(self):
        """Normalize profile to the library size and to the gene count i.e.
        normalize each window to the total count in it by dividing each
        entry by the sum of the counts.
        """
        adjusted_profile = asarray(map(lambda i: i/float(i.sum()), self.profile_raw))
        self.profile_normalized_to_gene = nan_to_num(adjusted_profile)

    def get_profile_by_type(self, profile_type):
        if profile_type == "raw":
            return self.profile_raw
        elif profile_type == "normalized_to_gene":
            if self.profile_normalized_to_gene is not None:
                return self.profile_normalized_to_gene
            else:
                raise NoSuchAProfileTypeException("Profile is not normalized by default. "
                                            "First run Profile.normalize_to_gene to normalize it")
        elif profile_type == "normalized_to_library":
            if self.profile_normalized_to_gene is not None:
                return self.profile_normalized_to_library
            else:
                raise NoSuchAProfileTypeException("Profile is not normalized by default. "
                                            "First run Profile.normalize_to_library to normalize it")
        else:
            raise NoSuchAProfileTypeException("No profile type %s. " % profile_type + \
                    "Available types are: raw, normalized_to_gene and normalized_to_library")

    def get_aggregated_profile(self, profile_type, metric='mean', axis=0):
        """Get the profile that is aggregated i.e. summarized along
        all windows.
        :param profile_type: @todo
        :param metric: @todo
        :returns: @todo

        """
        profile_to_aggregate = self.get_profile_by_type(profile_type=profile_type)
        return self.__agregate(profile_to_aggregate, metric, axis)

    def plot_line(self, profile_type, x=None, ax=None, smooth_by=1, aggregate=False, metric="mean", **kwargs):
        if ax is None:
            fig, ax = pl.subplots()
        if aggregate:
            profile_to_plot = self.get_aggregated_profile(profile_type=profile_type, metric=metric)
        else:
            profile_to_plot = self.get_profile_by_type(profile_type=profile_type)

        return plot_line_profile(profile_to_plot, x=x, ax=ax, smooth_by=smooth_by, **kwargs)

    def plot_heatmap(self, profile_type, ax=None, xlabel='', ylabel='', title='', percentiles=(1, 99), sort_by=None, **imshow_kws):
        if ax is None:
            fig, ax = pl.subplots()
        profile_to_plot = self.get_profile_by_type(profile_type=profile_type)
        return plot_heatmap_profile(profile_to_plot, ax=ax, xlabel=xlabel, ylabel=ylabel, title=title,
                                    percentiles=percentiles, sort_by=sort_by, **imshow_kws)

    def __agregate(self, profile, metric, axis=0):
        if metric.upper() == "MEAN" or metric.upper() == "AVERAGE":
            return profile.mean(axis=axis)
        elif metric.upper() == "SUM":
            return profile.sum(axis=axis)
        else:
            raise NoSuchAMetricException("%s is not valid metric. Available are sum or mean.")

    def __repr__(self):
        return "Profile of %s on %s" % (self.signal.name, self.window.name)
