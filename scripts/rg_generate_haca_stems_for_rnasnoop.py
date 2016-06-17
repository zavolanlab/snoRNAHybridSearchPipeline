#!/usr/bin/env python
"""
Generate fasta files for PLEXY from snoRNA input
"""

__date__ = "2015-02-20"
__author__ = "Rafal Gumienny"
__email__ = "r.gumienny@unibas.ch"
__license__ = "GPL"

# imports
import os
import sys
file_dir = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.join(file_dir, ".."))
import time
from modules.snoRNA import CD_snoRNA, WrongCDBoxPlacementException, IncompatibleStrandAndCoordsException
import pandas as pd
from argparse import ArgumentParser, RawTextHelpFormatter

class WrongTypeException(Exception): pass

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
                    help="Input file in tab format.")
parser.add_argument("--type",
                    dest="type",
                    required=True,
                    choices=("CD", "HACA"),
                    help="Type of snoRNA")
parser.add_argument("--dir",
                    dest="dir",
                    default="Plexy",
                    help="Directory to put output , defaults to Plexy")
parser.add_argument("--switch-boxes",
                    dest="switch_boxes",
                    action="store_true",
                    default=False,
                    help="If the CD box is located wrongly it will try to relabel it")


# redefine a functions for writing to stdout and stderr to save some writting
syserr = sys.stderr.write
sysout = sys.stdout.write


def main(options):
    """Main logic of the script"""
    snoRNAs = read_snoRNAs_from_table(options.input, type_of_snor=options.type, only_with_box=True)
    for s in snoRNAs:
        with open(os.path.join(options.dir, snor.snor_id + ".fa"), "w") as o:
            o.write(s.get_plexy_string())

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
