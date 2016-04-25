#!/usr/bin/env python
"""

"""

__date__ = "2016-02-23"
__author__ = "Rafal Gumienny"
__email__ = "r.gumienny@unibas.ch"
__license__ = "GPL"

# imports
import sys
import time
import nwalign
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
parser.add_argument("--a",
                    dest="a",
                    required=True,
                    help="First set of sequences")
parser.add_argument("--b",
                    dest="b",
                    required=True,
                    help="Second set of sequences")


try:
    options = parser.parse_args()
except Exception, e:
    parser.print_help()

# redefine a functions for writing to stdout and stderr to save some writting
syserr = sys.stderr.write
sysout = sys.stdout.write


def main():
    """Main logic of the script"""
    data = {}
    for arec in SeqIO.parse(options.a, 'fasta'):
        data[str(arec.id)] = {}
        for brec in SeqIO.parse(options.b, 'fasta'):
            alignment = nwalign.global_align(str(arec.seq), str(brec.seq), gap_open=-10, gap_extend=-4, matrix='matrix.txt')
            data[str(arec.id)][str(brec.id)] = nwalign.score_alignment(alignment[0], alignment[1], gap_open=-10, gap_extend=-4, matrix='matrix.txt')

    df = pd.DataFrame(data)
    for idx, row in df.iterrows():
        pos = idx.split("|")
        # ENSMUSG00000099280_1|chr7|-|121040344|121040480|0|0
        pos_id = pos[1] + ":" + pos[3] + "-" + pos[4]
        if row.max() > 10:
            bed_line = "%s\t%s\t%s\t%s\t%s\t%s" % (pos[1],
                                                   pos[3],
                                                   pos[4],
                                                   pos[0],
                                                   str(row.max()),
                                                   pos[2]
                                                   )
            aligned_to = ",".join(set([i.split("-")[0] for i in row[row == row.max()].index]))
            print bed_line + "\t" + aligned_to
            # print "%s\t%s\t%i\t%s" % (idx, , row.max(), pos)

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
