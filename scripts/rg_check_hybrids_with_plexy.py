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
parser.add_argument("--plexy-tmp",
                    dest="plexy_tmp",
                    default="temp/",
                    help="Plexy temporary directory , defaults to temp/")
parser.add_argument("--plexy-bin",
                    dest="plexy_bin",
                    required=True,
                    help="Path to PLEXY binary")

try:
    options = parser.parse_args()
except Exception, e:
    parser.print_help()

# redefine a functions for writing to stdout and stderr to save some writting
syserr = sys.stderr.write
sysout = sys.stdout.write


def main():
    nan = "NaN"
    """Main logic of the script"""
    snorna_paths = options.snoRNA_paths
    snorna_paths_glob = glob.glob(os.path.join(snorna_paths, "*.fa"))

    counter = 0
    input_name = options.input + ".plexyinput"
    with open(options.input) as infile, open(options.output, 'w') as outfile:
        for chrom, start, end, seqid, count, strand, sequence in csv.reader(infile, delimiter="\t"):
            # syserr("%i\n" % counter)
            snorna = seqid.split(":")[-2]
            row = "\t".join([chrom, start, end, seqid, count, strand, snorna, sequence])
            counter += 1
            with open(input_name, 'w') as infasta:
                infasta.write(">%s\n%s\n" % ("%s|%s|%s|%s|%s|%s" % (chrom, snorna, seqid, start, end, sequence), sequence))
            proc = subprocess.Popen("%s -f %s.fa -o %s -T %s -e 0.00" % (options.plexy_bin,
                                                                         os.path.join(snorna_paths, snorna),
                                                                         options.plexy_tmp,
                                                                         input_name),
                                    stdout=subprocess.PIPE, shell=True)
            proc_shuf = subprocess.Popen("%s -f %s -o %s -T %s -e 0.00" % (options.plexy_bin,
                                                                           random.choice(snorna_paths_glob),
                                                                           options.plexy_tmp,
                                                                           input_name),
                                         stdout=subprocess.PIPE, shell=True)
            sout, serr = proc.communicate()
            sout_shuf, serr_shuf = proc_shuf.communicate()
            try:
                myres_shuf = [re.split("\s*", i.rstrip()) for i in sout_shuf.split("\n") if not i.startswith("#")]
                myres_shuf = pd.DataFrame(myres_shuf).dropna()
                myres_shuf[3] = myres_shuf[3].astype(float)
                score_shuf = myres_shuf.sort(3).iloc[0][3]
            except KeyError:
                score_shuf = 0.0
            try:
                myres = [re.split("\s*", i.rstrip()) for i in sout.split("\n") if not i.startswith("#")]
                myres = pd.DataFrame(myres).dropna()
                myres[3] = myres[3].astype(float)
                results = myres.sort(3).iloc[0].tolist()
                modification_position = get_modification_position(int(results[2].split("-")[-1]), int(start), int(end), strand)
                outfile.write("%s\t%s\n" % (row,
                                     "\t".join([str(i) for i in [results[3], results[4], results[5], results[1], modification_position, score_shuf]])))
                # print "\t".join([str(i) for i in myres.sort(3).iloc[0].tolist()] + [str(score_shuf)])
            except KeyError:
                outfile.write("%s\t%s\n" % (row,
                                     "\t".join([str(i) for i in ["0.0", nan, nan, nan, nan, score_shuf]])))
                # print "\t".join([snorna, nan, seqid, "0.0", nan, nan, nan, str(score_shuf)])
            if options.verbose:
                syserr("%s\n" % serr)

def get_modification_position(rel_pos, start, end, strand):
    if strand == "-":
        return end - rel_pos
    else:
        return start + rel_pos


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
