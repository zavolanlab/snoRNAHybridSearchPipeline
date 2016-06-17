#!/usr/bin/env python
"""
Convert output table from alignment search into fasta
"""

__date_ = "2014-06-06"
__author__ = "Rafal Gumienny"
__email__ = "r.gumienny@unibas.ch"
__license__ = "GPL"

# imports
import re
import sys
import csv
import time
import swalign
from argparse import ArgumentParser

parser = ArgumentParser(description=__doc__)
parser.add_argument("-v",
                    "--verbose",
                    dest="verbose",
                    action="store_true",
                    default=False,
                    help="Be loud!")
parser.add_argument("--input",
                    dest="input",
                    help="Input table")
parser.add_argument("--output",
                    dest="output",
                    help="Output fasta file")
parser.add_argument("--stats",
                    dest="stats",
                    help="")
parser.add_argument("--length",
                    dest="length",
                    type=int,
                    default=15,
                    help="Length of the target site to keep, defaults to 15")
parser.add_argument("--assign-score-threshold",
                    dest="assign_score_threshold",
                    action="store_true",
                    default=False,
                    help="")
parser.add_argument("--filter-ambiguous",
                    dest="filter_ambiguous",
                    action="store_true",
                    default=False,
                    help="Filter reads that can be assigned to more than one snoRNA")
parser.add_argument("--five-prime-adapter",
                    dest="five_prime_adapter",
                    help="Five prime adapter sequence used in experiment - will be used to remove reads that are similar")
parser.add_argument("--three-prime-adapter",
                    dest="three_prime_adapter",
                    help="Three prime adapter sequence used in experiment - will be used to remove reads that are similar")
parser.add_argument("--five-prime-adapter-threshold",
                    dest="five_prime_adapter_threshold",
                    type=float,
                    default=0.8,
                    help="Threshold of the identity to the 5' adapter, defaults to 0.8")
parser.add_argument("--three-prime-adapter-threshold",
                    dest="three_prime_adapter_threshold",
                    type=float,
                    default=0.8,
                    help="Threshold of the identity to the 3' adapter, defaults to 0.8")


# redefine a functions for writing to stdout and stderr to save some writting syserr = sys.stderr.write
sysout = sys.stdout.write
syserr = sys.stderr.write


def main(options):
    """Main logic of the script"""
    if options.assign_score_threshold:
        with open(options.stats) as s:
            stats = {}
            for row in csv.reader(s, delimiter="\t"):
                # print row
                stats[row[0]] = row[1]
            score_thersh = float(stats["Score threshold"])
    else:
        score_thersh = 0.0

    if options.five_prime_adapter or options.three_prime_adapter:
        scoring = swalign.NucleotideScoringMatrix(match=1,
                                                  mismatch=-1)
        sw = swalign.LocalAlignment(scoring,
                                    gap_penalty=-6,
                                    gap_extension_penalty=-4)

    with open(options.input) as f:
        with open(options.output, 'w') as out:
            if options.filter_ambiguous:
                for row in csv.reader(f, delimiter='\t'):
                    if row[6] == "NA":
                        continue
                    if len(row[6]) >= options.length and float(row[2]) >= score_thersh and row[7] == "Unique":
                        if options.five_prime_adapter is not None:
                            alignment_five = sw.align(options.five_prime_adapter, row[6])
                            score_five = alignment_five.score
                            normed_five = score_five/float(len(row[6]))
                        else:
                            normed_five = 0.0
                        if options.three_prime_adapter is not None:
                            alignment_three = sw.align(options.three_prime_adapter, row[6])
                            score_three = alignment_three.score
                            normed_three = score_three/float(len(row[6]))
                        else:
                            normed_three = 0.0

                        if normed_three < options.three_prime_adapter_threshold and normed_five < options.five_prime_adapter_threshold:
                            out.write(">%s:N\n%s\n" % (row[0], row[6]))
                            for match in re.finditer("C", row[6]):
                                out.write(">%s\n%s\n" % (row[0] + ":" + str(match.span()[0]),
                                                         row[6][:match.span()[0]] + "T" + row[6][match.span()[1]:]))
            else:
                counter_adapt = 0
                counter = 0
                for row in csv.reader(f, delimiter='\t'):
                    if row[6] == "NA":
                        continue
                    if len(row[6]) >= options.length and float(row[2]) >= score_thersh:
                        counter += 1
                        if options.five_prime_adapter is not None:
                            alignment_five = sw.align(options.five_prime_adapter, row[6])
                            score_five = alignment_five.score
                            normed_five = score_five/float(len(row[6]))
                        else:
                            normed_five = 0.0
                        if options.three_prime_adapter is not None:
                            alignment_three = sw.align(options.three_prime_adapter, row[6])
                            score_three = alignment_three.score
                            normed_three = score_three/float(len(row[6]))
                        else:
                            normed_three = 0.0

                        if normed_three < options.three_prime_adapter_threshold and normed_five < options.five_prime_adapter_threshold:
                            counter_adapt += 1
                            out.write(">%s:N\n%s\n" % (row[0], row[6]))
                            for match in re.finditer("C", row[6]):
                                out.write(">%s\n%s\n" % (row[0] + ":" + str(match.span()[0]),
                                                         row[6][:match.span()[0]] + "T" + row[6][match.span()[1]:]))

                if options.verbose:
                    syserr("%i seqs was removed out of %s because they were similar to the adapter\n" % (counter - counter_adapt, counter))


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

# hsa-let-7a_@ago2_kishore-hc_3834205_x1  hsa-let-7a-5p   44.000000   0   22  5p  TTATTACGTCTTCTTC
# hsa-let-7a_@ago2_kishore-hc_1451499_x1  hsa-let-7a-5p   40.000000   0   20  5p  ACGTATGTGGTGTTCTGT
# hsa-let-7a_@ago2_kishore-pc_4275663_x1  hsa-let-7a-5p   40.000000   0   20  5p  ACTCTTCTATACAACTTG

