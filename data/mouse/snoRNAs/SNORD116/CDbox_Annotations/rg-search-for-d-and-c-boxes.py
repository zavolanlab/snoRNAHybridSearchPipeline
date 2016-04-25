#!/usr/bin/env python
"""
Search for D and C boxes in snoRNAs
"""

__date_ = "2014-09-09"
__author__ = "Rafal Gumienny"
__email__ = "r.gumienny@unibas.ch"
__license__ = "GPL"

# imports
import os
import sys
import time
import glob
from Bio import motifs, SeqIO, Seq, Alphabet
from argparse import ArgumentParser, RawTextHelpFormatter

parser = ArgumentParser(description=__doc__, formatter_class=RawTextHelpFormatter)
parser.add_argument("-v",
                    "--verbose",
                    dest="verbose",
                    action="store_true",
                    default=False,
                    help="Be loud!")
parser.add_argument("--seqs",
                    dest="seqs",
                    required=True,
                    help="Sequences to search in")

try:
    options = parser.parse_args()
except Exception, e:
    parser.print_help()

# redefine a functions for writing to stdout and stderr to save some writting
syserr = sys.stderr.write
sysout = sys.stdout.write


def main():
    """Main logic of the script"""
    with open("Cbox.fa") as f:
        instances_cbox = [rec.seq for rec in SeqIO.parse(f, 'fasta', Alphabet.IUPAC.unambiguous_dna)]
    with open("Dbox.fa") as f:
        instances_dbox = [rec.seq for rec in SeqIO.parse(f, 'fasta', Alphabet.IUPAC.unambiguous_dna)]
    with open("Cprime.fa") as f:
        instances_cpri = [rec.seq for rec in SeqIO.parse(f, 'fasta', Alphabet.IUPAC.unambiguous_dna)]
    with open("Dprime.fa") as f:
        instances_dpri = [rec.seq for rec in SeqIO.parse(f, 'fasta', Alphabet.IUPAC.unambiguous_dna)]
    with open(options.seqs) as f:
        sequences = {str(rec.id): rec.seq for rec in SeqIO.parse(f, 'fasta', Alphabet.IUPAC.unambiguous_dna)}


    cbox = motifs.create(instances_cbox)
    dbox = motifs.create(instances_dbox)
    cpri = motifs.create(instances_cpri)
    dpri = motifs.create(instances_dpri)

    # cbox.weblogo("Cbox.pdf", format='pdf', version='3.3')
    # dbox.weblogo("Dbox.pdf", format='pdf', version='3.3')
    # cpri.weblogo("Dpri.pdf", format='pdf', version='3.3')
    # dpri.weblogo("Cpri.pdf", format='pdf', version='3.3')

    snornas = [os.path.basename(f).split(".")[0] for f in glob.glob("../Annotated/*.fa")]
    print snornas

    for seqid, seq in sequences.iteritems():
        if seqid in snornas:
            continue
        print ">%s" % seqid
        print seq

        positions = []
        starts = []
        for position, score in cbox.pssm.search(seq, threshold=7.0, both=False):
            # print("Cbox position %d: score = %5.3f" % (position, score))
            positions.extend(range(position, position + cbox.length))
            starts.append(position + 1)
        text = ""
        for i in range(len(seq)):
            if i in positions:
                text += "C"
            else:
                text += "."
        print text, starts

        positions = []
        starts = []
        for position, score in dbox.pssm.search(seq, threshold=7.0, both=False):
            # print("Dbox position %d: score = %5.3f" % (position, score))
            positions.extend(range(position, position + dbox.length))
            starts.append(position + 1)
        text = ""
        for i in range(len(seq)):
            if i in positions:
                text += "D"
            else:
                text += "."
        print text, starts

        positions = []
        starts = []
        for position, score in cpri.pssm.search(seq, threshold=5.0, both=False):
            # print("C' position %d: score = %5.3f" % (position, score))
            positions.extend(range(position, position + cpri.length))
            starts.append(position + 1)
        text = ""
        for i in range(len(seq)):
            if i in positions:
                text += "c"
            else:
                text += "."
        print text, starts

        positions = []
        starts = []
        for position, score in dpri.pssm.search(seq, threshold=4.0, both=False):
            # print("D' position %d: score = %5.3f" % (position, score))
            positions.extend(range(position, position + dpri.length))
            starts.append(position + 1)
        text = ""
        for i in range(len(seq)):
            if i in positions:
                text += "d"
            else:
                text += "."
        print text, starts

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
