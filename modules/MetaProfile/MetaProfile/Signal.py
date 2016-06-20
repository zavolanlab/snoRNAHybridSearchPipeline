import os
import pybedtools
import multiprocessing
from sys import stderr
import HTSeq
from .utils import infer_extension
from numpy import fromiter, zeros

from bx.bbi.bigwig_file import BigWigFile

class WrongFiletypeException(Exception): pass
class NotSuchAFileException(Exception): pass
class NoSignalsException(Exception): pass


def get_signals(signals_list):
    """Get Signals based on the provided list of dicts. Use multiprocessing.

    :param signals_list: list of dicts -- a list with dictionaries of signals to take into account [{'argname': arg_value}]
    :param number_of_processes: int -- number of parallel processes
    :returns: list of Signals

    Example:
        >>> from MetaPy import Signal
        >>> signals_list = ({'name': "input1",
        >>>                  'filepath': "/path/to/input1.bed"},
        >>>                 {'name': "input2",
        >>>                  'filepath': "/path/to/input2.bed"})
        >>> signals = Signal.get_signals(signals_list)

    """
    results = []
    for params in signals_list:
        results.append(Signal(**params))
    return results


class Signal(object):

    """Genomic signal representation as HTSeq.GenomicArray"""

    def __init__(self, name, filepath, filetype=None, stranded=True, func=None):
        """
        Create genomic Signal

        :param name: str -- name of the signal
        :param filepath: str -- path to file with signals
        :param filetype: str -- type of the file (can be one of the following: BED, GTF, GFF, BAM, SAM or OTHER),
                                if the OTHER is specified a func parameter has to be provided
        :param stranded: bool -- specify whether the input is stranded
        :param func: function -- iterator function that takes filepath as argument and in each iteration it returns an
                                 object that has iv attribute that is an HTSeq.GenomicInterval
        """
        self.name = name
        if not os.path.isfile(filepath):
            raise NotSuchAFileException("No file named %s" % filepath)
        if filetype is None:
            filetype = infer_extension(filepath)
        if filetype.upper() not in ["BED", "GFF", "GTF", "BAM", "SAM", "OTHER", "BW", "BIGWIG", "BG", "BEDGRAPH"]:
            raise WrongFiletypeException("Wrong file type, allowed: %s" % "(BED, GTF, GFF, BAM, SAM, BW, BIGWIG, BEDGRAPH, BG, OTHER)")
        if filetype.upper() == "OTHER":
            if func is not None:
                raise ValueError("If the file type is OTHER you have to provide func parameter!")
            if not hasattr(func, '__call__'):
                raise TypeError("The func argument shall be a function!")
        self.filepath = filepath
        self.filetype = filetype
        self.stranded = stranded
        self.__create_genomic_signals(stranded=stranded, func=func)


    def __create_genomic_signals(self, stranded=True, func=None, use_wrappers=True):
        """Prepares coverage as a HTSeq.GenomicArray

        :param filepath: path to file
        :param filetype: type of the file (can be bed etc.)
        """
        stderr.write("Creating %s signal. It may take few minutes...\n" % self.name)
        self.coverage = HTSeq.GenomicArray("auto", stranded=stranded, typecode="d")
        self.library_size = 0
        if self.filetype.upper() == "BED":
            if use_wrappers:
                self.coverage = BedWrapper(self.filepath)
            else:
                for line in HTSeq.BED_Reader(self.filepath):
                    self.coverage[line.iv] += 1
                    self.library_size += 1
        elif self.filetype.upper() == "GFF" or self.filetype.upper() == "GTF":
            for line in HTSeq.GFF_Reader(self.filepath):
                self.coverage[line.iv] += 1
                self.library_size += 1
        elif self.filetype.upper() == "SAM":
            for line in HTSeq.SAM_Reader(self.filepath):
                self.coverage[line.iv] += 1
                self.library_size += 1
        elif self.filetype.upper() == "BAM":
            if use_wrappers:
                raise NotImplementedError("Bam wrapper is not yet implemented!")
                self.coverage = BamWrapper(self.filetype)
            for line in HTSeq.BAM_Reader(self.filepath):
                self.coverage[line.iv] += 1
                self.library_size += 1
        elif (self.filetype.upper() == "BG") or (self.filetype.upper() == "BEDGRAPH"):
            raise NotImplementedError("BedGraph is not yet implemented!")
        elif (self.filetype.upper() == "BW") or (self.filetype.upper() == "BIGWIG"):
            self.coverage = BigWigWrapper(self.filepath)
        elif self.filetype.upper() == "OTHER":
            for line in func(self.filepath):
                self.coverage[line.iv] += 1
                self.library_size += 1
        else:
            assert False, "I should not be here!"

    def __repr__(self):
        return "<Signal: %s>" % self.name


class BigWigWrapper(object):

    """A wrapper for bx-python BigWig file"""

    def __init__(self, filepath):
        self.bw = BigWigFile(open(filepath))

    def __getitem__(self, iv):
        return self.bw.get_as_array(iv.chrom, iv.start, iv.end)


class BedWrapper(object):

    """Wrapper around bed file using pybedtools"""

    def __init__(self, filepath):
        """

        :param filepath: @todo

        """
        self.filepath = filepath
        self.bed = pybedtools.BedTool(self.filepath)
        self.__prepare_bed()

    def __prepare_bed(self):
        if not self.bed._tabixed():
            self.bed = self.bed.sort().tabix(in_place=False, force=False, is_sorted=True)

    def __getitem__(self, iv):
        coverage = zeros(iv.length, dtype='i')
        if iv.strand == "+" or iv.strand == "-":
            coverage = HTSeq.GenomicArray("auto", stranded=True, typecode="d")
            for line in self.bed.tabix_intervals("%s:%i-%i[%s]" % (iv.chrom, iv.start, iv.end, iv.strand)):
                # coverage[line.start - iv.start: iv.end - line.end] += 1
                coverage[HTSeq.GenomicInterval(line.chrom, line.start, line.end, line.strand)] += 1
            return fromiter(coverage[iv], dtype='i', count=iv.length)
            # return coverage
        else:
            raise Exception
            coverage = HTSeq.GenomicArray("auto", stranded=False, typecode="d")
            for line in self.bed.tabix_intervals("%s:%i-%i" % (iv.chrom, iv.start, iv.end)):
                coverage[HTSeq.GenomicInterval(line.chrom, line.start, line.end, line.strand)] += 1
            return fromiter(coverage[iv], dtype='i', count=iv.length)

