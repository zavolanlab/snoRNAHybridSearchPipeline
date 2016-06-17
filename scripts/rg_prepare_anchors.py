#!/usr/bin/env python
"""
Prepare anchor sequences from provided fasta
"""

__date_ = "2014-06-04"
__author__ = "Rafal Gumienny"
__email__ = "r.gumienny@unibas.ch"
__license__ = "GPL"

# imports
import sys
import time
import csv
from Bio import SeqIO
from collections import defaultdict
from argparse import ArgumentParser

parser = ArgumentParser(description=__doc__)
parser.add_argument("-v",
                    "--verbose",
                    dest="verbose",
                    action="store_true",
                    default=False,
                    help="Be loud!")
parser.add_argument("--fasta-to-anchor",
                    dest="fasta_to_anchor",
                    help="Fasta to anchor")
parser.add_argument("--anchor-length",
                    dest="anchor_length",
                    type=int,
                    default=12,
                    help="Anchor length, defaults to 12")
parser.add_argument("--output",
                    dest="output",
                    help="Output file name")
parser.add_argument("--expressed-snoRNAs",
                    dest="expressed_snoRNAs",
                    required=True,
                    help="A list with expressed snoRNAs with RPKMs in form of: snoR_ID\tRPKM")



# redefine a functions for writing to stdout and stderr to save some writting
syserr = sys.stderr.write
sysout = sys.stdout.write


def main(options):
    """Main logic of the script"""
    with open(options.expressed_snoRNAs) as esnors:
        snornas = [row[0] for row in csv.reader(esnors, delimiter="\t")]
    anchors = defaultdict(list)
    snor_counter = 0
    expressed_counter = 0
    with open(options.fasta_to_anchor) as f:
        for rec in SeqIO.parse(f, 'fasta'):
            snor_counter += 1
            if str(rec.id) not in snornas:
                syserr(" - %s is not expressed.\n" % rec.id)
                continue
            expressed_counter += 1
            for seqpart, winpos_str, winpos in slide_windows(str(rec.seq).upper().replace("U", "T"), options.anchor_length, 1):
                anchor_name = "%s" % (str(rec.id))
                anchors[seqpart].append(anchor_name)
    with open(options.output, 'w') as out:
        for anchor, names in anchors.iteritems():
            out.write("%s\t%s\n" % (anchor, ",".join(names)))
    syserr("%i snoRNAs was written out of %i provided" % (expressed_counter, snor_counter))


def slide_windows(seq, windowsize, slidesize):
    """
    Slide window on string in order to divide it into parts of equal windowsize length
    """
    for i in range(0, len(seq), slidesize):
        if len(seq[i:i+windowsize]) == windowsize:
            # yield a part of the sequence of length windowsize
            yield seq[i:i+windowsize], "("+str(i+1)+","+str(i+windowsize)+")", (i+1, i+windowsize)
        else:
            # if we reach end of the sequence return last fragment of size 50
            yield seq[-windowsize:], "("+str(len(seq)-windowsize)+","+str(len(seq))+")", (len(seq)-windowsize, len(seq))
            break

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
            syserr("############## Started script on %s ##############\n" % start_date)
        main(options)
        if options.verbose:
            syserr("### Successfully finished in %i seconds, on %s ###\n" % (time.time() - start_time, time.strftime("%d-%m-%Y at %H:%M:%S")))
    except KeyboardInterrupt:
        syserr("Interrupted by user after %i seconds!\n" % (time.time() - start_time))
        sys.exit(-1)
