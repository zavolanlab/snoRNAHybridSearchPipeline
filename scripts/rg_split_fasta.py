#!/usr/bin/env python
"""
Split fasta file into batches
"""

__date__ = "2015-05-22"
__author__ = "Rafal Gumienny"
__email__ = "r.gumienny@unibas.ch"
__license__ = "GPL"

# imports
import os
import sys
import time
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
parser.add_argument("--input",
                    dest="input",
                    default=sys.stdin,
                    help="Input file in fasta format. Defaults to sys.stdin.")
parser.add_argument("--output-dir",
                    dest="output_dir",
                    default="./",
                    help="Output directory for split files. Defaults to .")
parser.add_argument("--batch-size",
                    dest="batch_size",
                    type=int,
                    default=100,
                    help="Batch size to split, defaults to 100")
parser.add_argument("--prefix",
                    dest="prefix",
                    default="part_",
                    help="Prefix to file name , defaults to part_")
parser.add_argument("--suffix",
                    dest="suffix",
                    default="inputfasta",
                    help="Suffix (extension) to the file name , defaults to inputfasta")


# redefine a functions for writing to stdout and stderr to save some writting
syserr = sys.stderr.write
sysout = sys.stdout.write


def main(options):
    """Main logic of the script"""
    with smart_open(options.input) as f:
        fasta = SeqIO.parse(f, 'fasta')
        for i, batch in enumerate(batch_iterator(fasta, options.batch_size)):
                filename = os.path.join(options.output_dir,
                                        "%s%i.%s" % (options.prefix, i+1, options.suffix))
                handle = open(filename, 'w')
                SeqIO.write(batch, handle, "fasta")
                handle.close()

    if options.verbose:
        syserr("File was split into %i parts.\n" % (i + 1))


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

def batch_iterator(iterator, batch_size) :
    """Returns lists of length batch_size.

    This can be used on any iterator, for example to batch up
    SeqRecord objects from Bio.SeqIO.parse(...), or to batch
    Alignment objects from Bio.AlignIO.parse(...), or simply
    lines from a file handle.

    This is a generator function, and it returns lists of the
    entries from the supplied iterator.  Each list will have
    batch_size entries, although the final list may be shorter.

    Args:
        iterator (iterable): Usually a SeqIO.parse iterator
        batch_size (int): batch size

    Returns: iterable batch
    """
    entry = True #Make sure we loop once
    while entry :
        batch = []
        while len(batch) < batch_size :
            try :
                entry = iterator.next()
            except StopIteration :
                entry = None
            if entry is None :
                #End of file
                break
            batch.append(entry)
        if batch:
            yield batch

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
