import os
import re
import sys
import csv
import string
import optparse
import subprocess
from sys import exit
from Bio import SeqIO
from os.path import abspath
from collections import defaultdict

# parse arguments and store them in the variables
argparser = optparse.OptionParser()
argparser.add_option('-v', '--verbose', dest='verbose', action="store_true",
                     default=False)
argparser.add_option('--genome-dir', dest='genome_dir',
                     action='store', help='fasta with mRNA sequences')
argparser.add_option('--out', '-o', default='output.tab', dest='out',
                     action='store', help='output table')
argparser.add_option('--energies', '-e',
                     default='/import/bc2/home/zavolan/gumiennr/Data/ContrafoldEnergies/results.txt',
                     dest='energies', action='store', help='file with energies')
argparser.add_option('--coords', default='coords.tab', dest='coords',
                     action='store', help='file with target\
        sites positions, miRNA and target gene ID')
argparser.add_option(
    '--contextLen_L', type=int, dest='contextLen_L', default=0,
    action='store', help='length of the \
        context sequence downstream binding site to be unwinded')
argparser.add_option(
    '--contextLen_U', type=int, dest='contextLen_U', default=0,
    action='store', help='length of the \
        context sequence upstream binding site to be unwinded')
argparser.add_option('--context', type=int, dest='context', default=50,
                     action='store', help='length of the \
        context of the seed to be checked')

arguments, args = argparser.parse_args()
verbose = arguments.verbose
# redefine a functions for writing to stdout and stderr to save some writting
syserr = sys.stderr.write
sysout = sys.stdout.write

# Define runContrafold function to run MIRZA and obtaine score


def runContrafold(mRNAseq, mRNAID, unw_start, unw_end):
    # TODO add description of the function
    """"""
    Contrabin = 'contrafold'
    #Contrabin = '/import/bc2/home/zavolan/gumiennr/Soft/contrafold/src/contrafold'
    # open all files and write in fasta format
    try:
        mrnainput = open('contrafold_mrna_input.bpseq', 'w')
    except IOError:
        print 'Cannot write files for Contrafold input'
        exit()

    # get the path for the file
    mrna_path = abspath(mrnainput.name)

    # write sequence to the file
    bpseq = []
    try:
        for i in range(1, len(mRNAseq) + 1):
            if i >= unw_start + 1 and i <= unw_end:
                bpseq.append('%i\t%s\t%i' % (i, mRNAseq[i - 1], 0))
            else:
                bpseq.append('%i\t%s\t%i' % (i, mRNAseq[i - 1], -1))
        mrnainput.write("\n".join(bpseq))
        mrnainput.close()
    except IOError:
        'Cannot write into the MIRZA input files'
        exit()
    #print "\n".join(bpseq)
    # Run Contrafold with these files
    print "Calculating accessibility for %s" % mRNAID
    contra_command_unwind = Contrabin + ' predict ' + \
        mrna_path + ' --constraints --partition'
    contrarun_unwind = subprocess.Popen(
        contra_command_unwind, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout_unwind, stderr_unwind = contrarun_unwind.communicate()

    contra_command_wind = Contrabin + ' predict ' + mrna_path + ' --partition'
    contrarun_wind = subprocess.Popen(
        contra_command_wind, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout_wind, stderr_wind = contrarun_wind.communicate()

    try:
        val_wind = float(re.search(': (\d+\.?\d*)', stdout_wind).group(1))
        val_unwind = float(re.search(': (\d+\.?\d*)', stdout_unwind).group(1))
        assert val_wind >= val_unwind
        #print "Energy: ", val_wind, val_wind
        #print "OUT: ",stdout_unwind
        #print stderr_unwind
        #print val_unwind
        return val_unwind - val_wind, stdout_unwind
    except Exception, e:
        ratio = "NA"
        print "Cannot find value among contrafold results for: mRNA %s" % (mRNAID)
        print stdout_unwind
        print stderr_unwind
        print str(e)
        print "\n #################################################################### \n"
        return ratio, stdout_unwind


def clean():
    """Remove all files that were used during calculations"""
    rmrnafile = subprocess.Popen(
        'rm ' + 'contrafold_mrna_input.bpseq', shell=True)
    rmrnafile.communicate()


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
outfile.write('#siteID\tContraScore\n')

# Iterate through the binding coordinates to calculate their score
if verbose == True:
    print "Calculating Contrafold accessibility..."
count = 0
corlen = len(coords)


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
            coorLen = int(row[3]) - int(row[2])
            lowerIndex = arguments.context - arguments.contextLen_L
            upperIndex = arguments.context + arguments.contextLen_U + coorLen
            sequence = get_sequence(chromosome_sequence, int(row[2]) - arguments.context, int(row[3]) + arguments.context)
            # print upper_sequence, len(upper_sequence)
            if row[4] == "-":
                sequence = sequence.translate(translation_table)[::-1]
                if len(sequence) >= upperIndex and  lowerIndex >= 0 and len(sequence) == 2 * arguments.context + coorLen:
                    result, unwind_out = runContrafold(sequence, row[0], lowerIndex, upperIndex)

                    outtext = '%s,%s,%s,%s\t%s\n' % (row[0], row[1], row[2],
                                                     row[3], result)
                    outfile.write(outtext)
                else:
                    outtext = '%s,%s,%s,%s\t%s\n' % (row[0], row[1], row[2],
                                                     row[3], 'NA')
                    outfile.write(outtext)
            else:
                print row[0], row[1], row[2], len(sequence)
                if len(sequence) >= upperIndex and  lowerIndex >= 0 and len(sequence) == 2 * arguments.context + coorLen:
                    result, unwind_out = runContrafold(sequence, row[0], lowerIndex, upperIndex)

                    outtext = '%s,%s,%s,%s\t%s\n' % (row[0], row[1], row[2],
                                                     row[3], result)
                    outfile.write(outtext)
                else:
                    outtext = '%s,%s,%s,%s\t%s\n' % (row[0], row[1], row[2],
                                                     row[3], 'NA')
                    outfile.write(outtext)
    except IOError:
        syserr("Cannot read from %s file\n" % (path_to_chromosome))



outfile.close()
if verbose == True:
    print "Cleaning"
clean()
if verbose == True:
    print "Done!"
