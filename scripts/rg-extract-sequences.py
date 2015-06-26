#!/usr/bin/env python
"""
Given bed file extract sequences according to chromosome and strand and save it as additional
column in input file or fasta
"""

__date_ = "2014-09-24"
__author__ = "Rafal Gumienny"
__email__ = "r.gumienny@unibas.ch"
__license__ = "GPL"

# imports
import os
import sys
import csv
import time
import errno
import string
from Bio import SeqIO
from collections import defaultdict
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
                    default=sys.stdin,
                    help="Input file in Bed format.")
parser.add_argument("--output",
                    dest="output",
                    default=sys.stdout,
                    help="Output file in Bed format.")
parser.add_argument("--format",
                    dest="format",
                    choices=("bed", "fasta"),
                    default="bed",
                    help="Output format, defaults to bed")
parser.add_argument("--sequence-length",
                    dest="sequence_length",
                    type=int,
                    help="Final length of sequence to extract independently of coordinates.")
parser.add_argument("--genome-dir",
                    dest="genome_dir",
                    required=True,
                    help="Directory where the fasta sequences with all the chromosomes are stored")
parser.add_argument("--window-left",
                    dest="window_left",
                    type=int,
                    default=0,
                    help="Add nucleotides to the left (upstream). This option does not work if sequence-length is specified, defaults to 0")
parser.add_argument("--window-right",
                    dest="window_right",
                    type=int,
                    default=0,
                    help="Add nucleotides to the right (downstream). This option does not work if sequence-length is specified, defaults to 0")

try:
    options = parser.parse_args()
except Exception, e:
    parser.print_help()

# redefine a functions for writing to stdout and stderr to save some writting
syserr = sys.stderr.write
sysout = sys.stdout.write

class NoSuchAFileException(IOError): pass


def main():
    """Main logic of the script"""
    data = defaultdict(list)
    with smart_open(options.input) as infile:
        for row in csv.reader(infile, delimiter="\t"):
            data[row[0]].append(row)

    # make table for making reversed complement of the sequence
    # ACGTNRYWSMKBDHVacgtnrywsmkbdhv -> TGCANYRWSKMVDHBtgcanyrwskmvdhb
    translation_table = string.maketrans('ACGTNRYWSMKBDHVacgtnrywsmkbdhv', 'TGCANYRWSKMVDHBtgcanyrwskmvdhb')

    # iterate through data
    with smart_open(options.output, "w") as outfile:
        for chromosome, rows in data.iteritems():
            if options.verbose:
                syserr("Extracting sequences from %s\n" % chromosome)
            path_to_chromosome = os.path.join(options.genome_dir, chromosome + ".fa")
            try:
                with open(path_to_chromosome) as chromosome_file_handle:
                    chromosome_sequence = str(SeqIO.parse(chromosome_file_handle, 'fasta').next().seq)
            except IOError, e:
                syserr("No such a file: %s\n" % path_to_chromosome)
                continue
            for row in rows:
                old_beg, old_end = (int(row[1]), int(row[2]))
                if options.sequence_length:
                    target_sequence, beg, end = get_equal_length_sequence(chromosome_sequence,
                                                                          old_beg,
                                                                          old_end,
                                                                          options.sequence_length)
                    if row[5] == "-":
                        seq = target_sequence.translate(translation_table)[::-1]
                        if options.format == "bed":
                            outfile.write("%s\t%s\n" % ("\t".join(row[:1] + [str(beg),
                                                                             str(end)] + row[3:]),
                                                        seq))
                        elif options.format == "fasta":
                            header = "%s_%s|%s|%s|%i|%i|%i|%i" % (row[3],
                                                                  row[4],
                                                                  row[0],
                                                                  row[5],
                                                                  old_beg + 1, # 1-based
                                                                  old_end,
                                                                  old_beg - beg,
                                                                  end - old_end)
                            outfile.write(">%s\n%s\n" % (header, seq))
                        else:
                            raise Exception("No such a format: %s" % options.format)
                    else:
                        if options.format == "bed":
                            outfile.write("%s\t%s\n" % ("\t".join(row[:1] + [str(beg),
                                                                             str(end)] + row[3:]),
                                                        target_sequence))
                        elif options.format == "fasta":
                            header = "%s_%s|%s|%s|%i|%i|%i|%i" % (row[3],
                                                                  row[4],
                                                                  row[0],
                                                                  row[5],
                                                                  old_beg + 1, # 1-based
                                                                  old_end,
                                                                  old_beg - beg,
                                                                  end - old_end)
                            outfile.write(">%s\n%s\n" % (header, target_sequence))
                        else:
                            raise Exception("No such a format: %s" % options.format)
                else:
                    if row[5] == "-":
                        new_beg = old_beg - options.window_right
                        new_end = old_end + options.window_left
                        target_sequence = get_sequence(chromosome_sequence, new_beg, new_end)
                        seq = target_sequence.translate(translation_table)[::-1]
                        if options.format == "bed":
                            outfile.write("%s\t%s\n" % ("\t".join(row), seq))
                        elif options.format == "fasta":
                            header = "%s_%s|%s|%s|%i|%i|%i|%i" % (row[3],
                                                                  row[4],
                                                                  row[0],
                                                                  row[5],
                                                                  old_beg + 1, # 1-based
                                                                  old_end,
                                                                  options.window_left,
                                                                  options.window_right)
                            outfile.write(">%s\n%s\n" % (header, seq))
                        else:
                            raise Exception("No such a format: %s" % options.format)
                    else:
                        new_beg = old_beg - options.window_left
                        new_end = old_end + options.window_right
                        target_sequence = get_sequence(chromosome_sequence, new_beg, new_end)
                        if options.format == "bed":
                            outfile.write("%s\t%s\n" % ("\t".join(row), target_sequence))
                        elif options.format == "fasta":
                            header = "%s_%s|%s|%s|%i|%i|%i|%i" % (row[3],
                                                                  row[4],
                                                                  row[0],
                                                                  row[5],
                                                                  old_beg + 1, # 1-based
                                                                  old_end,
                                                                  options.window_left,
                                                                  options.window_right)
                            outfile.write(">%s\n%s\n" % (header, target_sequence))


def get_equal_length_sequence(sequence, start, end, final_length):
    assert start < end
    new_start, new_end = get_indices(final_length, start, end)
    new_start = new_start if new_start >= 0 else 0
    new_end = new_end if new_start > 0 else final_length
    if new_end > len(sequence):
        new_end = len(sequence)
        new_start = len(sequence) - final_length
    try:
        assert len(sequence[new_start: new_end]) == final_length
    except AssertionError:
        import ipdb; ipdb.set_trace()
    return sequence[new_start: new_end], new_start, new_end

def get_sequence(sequence, start, end):
    assert start < end
    return sequence[start: end]

def get_indices(desired_length, beg, end):
    """Calculate upper and lower index based on what size we want and what are the coordinates"""
    lower_index = beg - ((desired_length - (end - beg))//2 + (desired_length - (end - beg))%2)
    upper_index = end + ((desired_length - (end - beg))//2)
    return lower_index, upper_index


@contextmanager
def smart_open(filepath, mode='r'):
    if mode == 'r':
        if filepath is not sys.stdin:
            fh = open(filepath, 'r')
        else:
            fh = filepath
        try:
            yield fh
        except IOError as e:
            if e.errno == errno.EPIPE:
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
        except IOError as e:
            if e.errno == errno.EPIPE:
                pass
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
            syserr("############## Started script on %s ##############\n" % start_date)
        main()
        if options.verbose:
            syserr("### Successfully finished in %i seconds, on %s ###\n" % (time.time() - start_time, time.strftime("%d-%m-%Y at %H:%M:%S")))
    except KeyboardInterrupt:
        syserr("Interrupted by user after %i seconds!\n" % (time.time() - start_time))
        sys.exit(-1)
