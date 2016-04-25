#!/usr/bin/env python
"""
Clean RNAs from modomics
"""

__date__ = "2016-04-12"
__author__ = "Rafal Gumienny"
__email__ = "r.gumienny@unibas.ch"
__license__ = "GPL"

# imports
import sys
import csv
import time
import pandas as pd
from Bio import SeqIO
from contextlib import contextmanager
from argparse import ArgumentParser, RawTextHelpFormatter

parser = ArgumentParser(description=__doc__, formatter_class=RawTextHelpFormatter)
parser.add_argument("-v",
                    "--verbose",
                    dest="verbose",
                    action="store_true",
                    default=False,
                    help="Be loud!")
parser.add_argument("--RNAs",
                    dest="RNAs",
                    required=True,
                    help="RNA sequences in fasta format")
parser.add_argument("--modifications",
                    dest="modifications",
                    required=True,
                    help="File from modomics with modifications")

try:
    options = parser.parse_args()
except Exception, e:
    parser.print_help()

# redefine a functions for writing to stdout and stderr to save some writting
syserr = sys.stderr.write
sysout = sys.stdout.write


def main():
    """Main logic of the script"""
    modifications = pd.read_table(options.modifications,  quoting=csv.QUOTE_NONE, encoding='utf-8')
    modified_nucs = modifications.rnamods_abbrev.tolist()

    rnas = {}
    with open(options.RNAs) as rnas_file:
        for i, rec in enumerate(SeqIO.parse(rnas_file, 'fasta')):
            name = str(rec.description).decode('utf-8').replace(" ", '').split("|")[1] + \
                    translate_to_raw(str(rec.description).decode('utf-8').replace(" ", '').split("|")[2], modifications) + "_" + str(i + 1)
            rnas[name] = str(rec.seq).decode("utf-8").replace("-", "").replace("_", "")

    # for key, value in rnas.iteritems():
    #     print key, value, translate_to_raw(value, modifications)
    with open("tRNAs_raw.fa", 'w') as out:
        for name, seq in rnas.iteritems():
            for counter, letter in enumerate(seq, start=1):
                if letter in modified_nucs and letter not in ['A', 'C', 'G', 'U']:
                    original_base = modifications[modifications.rnamods_abbrev == letter].originating_base.iloc[0]
                    modification_name = modifications[modifications.rnamods_abbrev == letter].name.iloc[0]
                    modification_short_name = modifications[modifications.rnamods_abbrev == letter].short_name.iloc[0]
                    line = u"%s\t%i\t%s\t%s\t%s" % (name, counter, original_base, modification_name, modification_short_name)
                    print line.encode('utf-8')
            out.write(">%s\n%s\n" % (name, translate_to_raw(seq, modifications)))



def translate_to_raw(seq, modifications):
    modified_nucs = modifications.rnamods_abbrev.tolist()
    new_seq = ""
    for letter in seq:
        if letter in modified_nucs:
            nl = modifications[modifications.rnamods_abbrev == letter].originating_base.iloc[0]
            if pd.isnull(nl) or len(nl) > 1:
                nl = 'N'
            new_seq += nl
        else:
            new_seq += letter
    return new_seq


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
