#!/usr/bin/env python
"""
RNA5-8S5|NR_003285.2    15  30  SNORD16 1   +
RNA5-8S5|NR_003285.2    86  105 SNORD16 1   +
RNA28S5|NR_003287.2 1563    1582    SNORD56B    1   +

SNORD50A|chr7|+|57640816|57640830|20|20 SNORD50A    TCATGCTTTGTGTTGTGAAGACCGCCTGGGACTACCGGGCAGGGTGTAGTAGGCA
SNORD50A|chr7|+|68527467|68527482|20|20 SNORD50A    ACTGAAGAAATTCAGTGAAATGCGGGTAAACGGCGGGAGTAACTATGACTCTCTTA
SNORD50A|chr7|+|68527638|68527654|20|20 SNORD50A    AATCAGCGGGGAAAGAAGACCCTGTTGAGTTTGACTCTAGTCTGGCATGGTGAAGAG
"""

__date_ = "2014-08-26"
__author__ = "Rafal Gumienny"
__email__ = "r.gumienny@unibas.ch"
__license__ = "GPL"

# imports
import re
import os
import sys
import csv
import glob
import time
import random
import subprocess
import pandas as pd
from Bio import SeqIO
from ushuffle import shuffle
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
                    help="Input file in tab format.")
parser.add_argument("--output",
                    dest="output",
                    required=True,
                    help="Output file in tab format.")
parser.add_argument("--snoRNA-paths",
                    dest="snoRNA_paths",
                    default="./Plexy/",
                    help="Path to snoRNAs with Plexy , defaults to ./Plexy/")
parser.add_argument("--RNAduplex-bin",
                    dest="RNAduplex_bin",
                    default="RNAduplex",
                    help="Path to RNAduplex binary, defaults to RNAcofold")


# redefine a functions for writing to stdout and stderr to save some writting
syserr = sys.stderr.write
sysout = sys.stdout.write


def main(options):
    nan = "NaN"
    """Main logic of the script"""
    snorna_paths = options.snoRNA_paths
    snorna_paths_glob = glob.glob(os.path.join(snorna_paths, "*.fa"))

    counter = 0
    with open(options.input) as infile, open(options.output, 'w') as outfile:
        for row in csv.reader(infile, delimiter="\t"):
            chrom, start, end, seqid, count, strand, snorna, sequence = row[:8]
            # syserr("%i\n" % counter)
            snorna = seqid.split(":")[-2]
            counter += 1
            with open(os.path.join(snorna_paths, snorna + '.fa')) as snorna_in:
                snorna_sequence = str(SeqIO.parse(snorna_in, 'fasta').next().seq)
            with open(random.choice(snorna_paths_glob)) as snorna_in:
                snorna_sequence_random = str(SeqIO.parse(snorna_in, 'fasta').next().seq)

            shuffled_target = shuffle(sequence, len(sequence), 2)
            shuffled_snorna = shuffle(snorna_sequence, len(snorna_sequence), 2)

            proc = subprocess.Popen('printf "%s\n%s" | %s' % (snorna_sequence, sequence, options.RNAduplex_bin),
                                    stdout=subprocess.PIPE, shell=True)
            proc_shuf = subprocess.Popen('printf "%s\n%s" | %s' % (snorna_sequence_random, sequence, options.RNAduplex_bin),
                                    stdout=subprocess.PIPE, shell=True)
            proc_shuf_tar = subprocess.Popen('printf "%s\n%s" | %s' % (snorna_sequence, shuffled_target, options.RNAduplex_bin),
                                    stdout=subprocess.PIPE, shell=True)
            sout, serr = proc.communicate()
            sout_shuf, serr_shuf = proc_shuf.communicate()
            sout_shuf_tar, serr_shuf_tar = proc_shuf_tar.communicate()

            # .(.((.(((((((((((((((((((....(((((&.)))).)....)))).))))))))))))))))).).  17,50  :  35,70  (-29.20)

            score_shuf = sout_shuf.split()[-1][1:-1]
            score_shuf_tar = sout_shuf_tar.split()[-1][1:-1]
            score = sout.split()[-1][1:-1]
            part1 = sout.split()[1]
            part2 = sout.split()[3]
            structure = get_complete_structure(snorna_sequence, sout.split()[0].split("&")[0], part1)
            structure_random = get_complete_structure(snorna_sequence_random, sout_shuf.split()[0].split("&")[0], sout_shuf.split()[1])
            outfile.write("%s\t%s\n" % ("\t".join(row), "\t".join([str(i) for i in [len(snorna_sequence), get_GC_frac(snorna_sequence), score, score_shuf, score_shuf_tar, structure, part1, part2, structure_random]])))

def get_complete_structure(snoRNA_sequence, partial_structure, positions, index=None):
    start, stop = [int(i) for i in positions.split(',')] # 1-based
    struc = ("#"*len(snoRNA_sequence[: start - 1])
             + partial_structure
             + "#"*len(snoRNA_sequence[stop:]))
    assert len(struc) == len(snoRNA_sequence), "%s, %i, %i, %s" % (positions, len(snoRNA_sequence), len(struc),
                                                                  str(index))
    return struc

def get_GC_frac(seq):
    seq = seq.upper().replace("U", "T")
    Gcon = seq.count("G")
    Ccon = seq.count("C")
    return (Gcon + Ccon)/float(len(seq))

def get_modification_position(rel_pos, start, end, strand):
    if strand == "-":
        return end - rel_pos
    else:
        return start + rel_pos


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
