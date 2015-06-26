#!/import/bc2/soft/app/Python/2.7.3/Linux/bin/python2.7

import optparse
import os
import csv
from Bio import SeqIO
from sys import exit
import sys
from os.path import abspath
import subprocess
import re
import string
from collections import defaultdict

# parse arguments and store them in the variables
argparser = optparse.OptionParser()
argparser.add_option('-v', '--verbose', dest='verbose', action="store_true",
                     default=False)
argparser.add_option('--genome-dir', dest='genome_dir',
                     action='store', help='fasta with mRNA sequences')
argparser.add_option('--out', '-o', default='output.tab', dest='out',
                     action='store', help='output table')
argparser.add_option('--coords', default='coords.tab', dest='coords',
                     action='store', help='file with target\
        sites positions, miRNA and target gene ID')
argparser.add_option('--contextLen', type=int, dest='contextLen', default=50,
                     action='store', help='length of the \
        context sequence serounding binding site')

arguments, args = argparser.parse_args()
verbose = arguments.verbose

# redefine a functions for writing to stdout and stderr to save some writting
syserr = sys.stderr.write
sysout = sys.stdout.write

def calculate_flanks_composition(lower_seq, upper_seq):
    extr_seq = lower_seq + upper_seq
    extr_seq = extr_seq.upper().replace('U', 'T')
    G_count = float(extr_seq.count('G'))
    A_count = float(extr_seq.count('A'))
    C_count = float(extr_seq.count('C'))
    T_count = float(extr_seq.count('T'))
    g_perc = str(G_count / float(len(extr_seq)))
    a_perc = str(A_count / float(len(extr_seq)))
    c_perc = str(C_count / float(len(extr_seq)))
    t_perc = str(T_count / float(len(extr_seq)))
    return g_perc, a_perc, c_perc, t_perc


def get_sequence(sequence, start, end):
    assert start < end
    return sequence[start: end]


# Read coords file into the table: [geneID, miR name, begining position,
# end position]

coords = defaultdict(list)
with open(arguments.coords) as infile:
    for row in csv.reader(infile, delimiter="\t"):
        coords[row[0]].append(row)



# Open output file and write first lines
try:
    outfile = open(arguments.out, 'w')
except IOError:
    print "Connot open output file %s" % (arguments.out)
outfile.write('#siteID\tflanksG\tflanksA\tflanksC\tflanksU\n')

# Iterate through the binding coordinates to calculate their score
clen = float(len(coords))
translation_table = string.maketrans('ACGTNRYWSMKBDHVacgtnrywsmkbdhv', 'TGCANYRWSKMVDHBtgcanyrwskmvdhb')
for chromosome, rows in coords.iteritems():
    # here we assume that the coordinates are given in 0-based for start and
    # 1-based for end
    try:
        syserr("Extracting sequences from %s\n" % chromosome)
        path_to_chromosome = os.path.join(arguments.genome_dir, chromosome + ".fa")
        with open(path_to_chromosome) as chromosome_file_handle:
            chromosome_sequence = str(SeqIO.parse(chromosome_file_handle, 'fasta').next().seq)
        for row in rows:
            # print row
            lower_sequence = get_sequence(chromosome_sequence, int(row[2]) - arguments.contextLen, int(row[2]))
            upper_sequence = get_sequence(chromosome_sequence, int(row[3]), int(row[3]) + arguments.contextLen)
            # print lower_sequence, len(lower_sequence)
            # print upper_sequence, len(upper_sequence)
            if row[4] == "-":
                lower_sequence = lower_sequence.translate(translation_table)[::-1]
                upper_sequence = upper_sequence.translate(translation_table)[::-1]
                if int(row[2]) - arguments.contextLen >= 0 and int(row[3]) <= len(chromosome_sequence):
                    gcont, acont, ccont, tcont = calculate_flanks_composition(lower_sequence, upper_sequence)
                    outtext = '%s,%s,%s,%s\t%s\t%s\t%s\t%s\n' % (row[0],
                                                                   row[1],
                                                                   row[2],
                                                                   row[3],
                                                                   gcont,
                                                                   acont,
                                                                   ccont,
                                                                   tcont)
                    outfile.write(outtext)
                else:
                    sys.stderr.write("Flanks out of boarders\n")
                    outtext = '%s,%s,%s,%s\t%s\t%s\t%s\t%s\n' % (row[0],
                                                                   row[1],
                                                                   row[2],
                                                                   row[3],
                                                                   'NA',
                                                                   'NA',
                                                                   'NA',
                                                                   'NA')
                    outfile.write(outtext)
            else:
                if int(row[2]) - arguments.contextLen >= 0 and int(row[3]) <= len(chromosome_sequence):
                    gcont, acont, ccont, tcont = calculate_flanks_composition(lower_sequence, upper_sequence)
                    outtext = '%s,%s,%s,%s\t%s\t%s\t%s\t%s\n' % (row[0],
                                                                   row[1],
                                                                   row[2],
                                                                   row[3],
                                                                   gcont,
                                                                   acont,
                                                                   ccont,
                                                                   tcont)
                    outfile.write(outtext)
                else:
                    sys.stderr.write("Flanks out of boarders\n")
                    outtext = '%s,%s,%s,%s\t%s\t%s\t%s\t%s\n' % (row[0],
                                                                   row[1],
                                                                   row[2],
                                                                   row[3],
                                                                   'NA',
                                                                   'NA',
                                                                   'NA',
                                                                   'NA')
                    outfile.write(outtext)
    except IOError:
        syserr("Cannot read from %s file\n" % (path_to_chromosome))

outfile.close()
