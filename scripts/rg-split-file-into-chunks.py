#!/usr/bin/env python
"""
Split text file into files with desired number of lines
"""

__date__ = "2015-07-03"
__author__ = "Rafal Gumienny"
__email__ = "r.gumienny@unibas.ch"
__license__ = "GPL"

# imports
import sys
import os
import time
import itertools
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
                    required=True,
                    help="Input file in txt format. Defaults to sys.stdin.")
parser.add_argument("--lines",
                    dest="lines",
                    type=int,
                    required=True,
                    help="Number of lines in each file")
parser.add_argument("--prefix",
                    dest="prefix",
                    default="file_",
                    help="Prefix to the file, defaults to file_")
parser.add_argument("--dir",
                    dest="dir",
                    default="./",
                    help="Directory to put files, defaults to ./")
parser.add_argument("--suffix",
                    dest="suffix",
                    default=".part",
                    help="Suffix to the file, defaults to .part")

try:
    options = parser.parse_args()
except Exception, e:
    parser.print_help()

# redefine a functions for writing to stdout and stderr to save some writting
syserr = sys.stderr.write
sysout = sys.stdout.write


def main():
    """Main logic of the script"""
    with open(options.input) as f:
        groups = itertools.groupby(f, key=lambda x, line=itertools.count(): next(line)//options.lines)
        counter = 0
        for k, group in groups:
            counter += 1
            with open(os.path.join(options.dir, "%s%i%s" % (options.prefix,
                                                             counter,
                                                             options.suffix)), 'w') as out:
                for line in group:
                    out.write(line)


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
