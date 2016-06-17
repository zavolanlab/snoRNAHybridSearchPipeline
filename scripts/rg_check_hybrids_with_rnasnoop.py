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
parser.add_argument("--rnasnoop",
                    dest="rnasnoop",
                    default="RNAsnoop",
                    help="Path to RNAsnoop binary, defaults to RNAsnoop")
parser.add_argument("--snoRNA-paths",
                    dest="snorna_paths",
                    # default="/import/bc2/home/zavolan/gumiennr/Pipelines/Pipelines/pipeline_snoRNASearch/data/HACA/snornas_stems",
                    required=True,
                    help="Path to snoRNAs stems")


# redefine a functions for writing to stdout and stderr to save some writting
syserr = sys.stderr.write
sysout = sys.stdout.write


def main(options):
    nan = "NaN"
    """Main logic of the script"""
    snorna_paths_glob = glob.glob(os.path.join(options.snorna_paths, "*.fa"))

    counter = 0
    input_name = options.input + ".rnasnoopinput"
    with open(options.input) as infile, open(options.output, 'w') as outfile:
        for chrom, start, end, seqid, count, strand, sequence in csv.reader(infile,
                                                                            delimiter="\t"):
            # syserr("%i\n" % counter)
            snorna = seqid.split(":")[-2]
            row = "\t".join([chrom, start, end, seqid, count, strand, snorna, sequence])
            counter += 1
            with open(input_name, 'w') as infasta:
                infasta.write(">%s\n%s\n" % ("%s" % (seqid),
                                             sequence))
            command_correct = "%s -s %s.fa -t %s" % (options.rnasnoop,
                                                     os.path.join(options.snorna_paths,
                                                                  snorna),
                                                     input_name)
            command_shuffled = "%s -s %s -t %s" % (options.rnasnoop,
                                                      random.choice(snorna_paths_glob),
                                                      input_name)
            # print "stem1", command_correct_stem1
            proc_stem1 = subprocess.Popen(command_correct,
                                          stdout=subprocess.PIPE, shell=True)
            sout_stem1, serr_stem1 = proc_stem1.communicate()
            # print "shuffled", command_shuffled
            proc_shuf = subprocess.Popen(command_shuffled,
                                         stdout=subprocess.PIPE, shell=True)
            sout_shuf, serr_shuf = proc_shuf.communicate()

            #
            # Get score for shuffled snoRNA
            #
            try:
                myres_shuf = [re.split("\s*", i.rstrip())
                               for i in sout_shuf.split("\n") if not i.startswith(">")]
                score_shuf = float(myres_shuf[0][6][1:])
            except ValueError:
                score_shuf = 0.0
            #
            # Get score for snoRNA
            #
            try:
                myres_stem1 = [re.split("\s*", i.rstrip())
                                for i in sout_stem1.split("\n") if not i.startswith(">")]
                score_stem1 = float(myres_stem1[0][6][1:])
                seq_stem1 = myres_stem1[1][0]
                modification_position_stem1 = get_modification_position(int(myres_stem1[0][3]),
                                                                        int(start),
                                                                        int(end),
                                                                        strand)
            except ValueError:
                score_stem1 = nan
                seq_stem1 = nan
                modification_position_stem1 = nan
            #
            # Get the position from the stem loop that has the lowest energy
            #
            outfile.write("%s\t%s\n" % (row, "\t".join([str(i) for i in [score_stem1,
                                                                         myres_stem1[0][0],
                                                                         seq_stem1,
                                                                         snorna.split("_")[-1],
                                                                         modification_position_stem1,
                                                                         score_shuf]])))
            # if options.verbose:
            #     syserr("%s\n" % serr)
        #
        # Clean postscript files in the end
        #
        ps_files = glob.glob("./*.ps")
        for ps_file in ps_files:
            print 'rm -f "%s"' % ps_file
            os.system('rm -f "%s"' % ps_file)

# ['<<.|.<<.<.<<.<.<.<<<<&.>>>>.>.>.>>.>.>>.((((((((((........)))).))))))>>............', '2,22', ';', '5', ':', '3,50', '(-15.80', '=', '-4.30', '+', '-5.90', '+', '-9.70', '+', '0.00', '+', '4.1', ')']

# RNA5-8S5    57  96  2249554_1:SNORD38B:2    30  +   SNORD38B    AGAATTAATGTGAATTGCAGGACACATTGATCATCGACA -11.6   .(((((((.&.))))))). TGCAmGGACA&AGTTCTGCT    D'  76  -2.9
# RNA5-8S5    58  97  2249555_3:SNORD38B:2    34  +   SNORD38B    GAATTAATGTGAATTGCAGGACACATTGATCATCGACAC -11.6   .(((((((.&.))))))). TGCAmGGACA&AGTTCTGCT    D'  76  0.0
# RNA5-8S5    58  97  2249558_3:SNORD38B:14   34  +   SNORD38B    GAATTAATGTGAATTGCAGGACACATTGATCATCGACAC -11.6   .(((((((.&.))))))). TGCAmGGACA&AGTTCTGCT    D'  76  0.0
# RNA5-8S5    58  97  2249559_5:SNORD38B:N    34  +   SNORD38B    GAATTAATGTGAATTGCAGGACACATTGATCATCGACAC -11.6   .(((((((.&.))))))). TGCAmGGACA&AGTTCTGCT    D'  76  -4.2


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
