#!/usr/bin/env python
"""
Convert Hadi annotations to input for Plexy
"""

__date_ = "2014-09-08"
__author__ = "Rafal Gumienny"
__email__ = "r.gumienny@unibas.ch"
__license__ = "GPL"

# imports
import re
import sys
import time
from itertools import groupby
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
                    help="Input file in Tab format.")
parser.add_argument("--output",
                    dest="output",
                    help="Output file in fasta format.")

try:
    options = parser.parse_args()
except Exception, e:
    parser.print_help()

# redefine a functions for writing to stdout and stderr to save some writting
syserr = sys.stderr.write
sysout = sys.stdout.write


def main():
    """Main logic of the script"""
    # The input looks like that:
    #
    # >snoRNA:SNORD116_ENSMUSG00000093869_92|chr7|-|59781781|59781873|0|0
    # GGATCTATGATGATTCCCAGTCAAACATTCCTTGGAAAAGCTGAACAAAATGAGTGAAAACTCTGTACCGCCACTCTCATCTGAACTGAGGTC
    # ......CCCCCCC...........................dddd.........ccccccc.........................DDDD....
    #
    # The output should look like that:
    # chr12 snoStrip    CD-snoRNA   109649617   109649687   "snoBoard:CD_10-1,family:SNORD114,boxes:ATTATGA_39_ATGA_28_ATGATGA_8_CTGA_62,seq:TGGATCAATGATGACTGTGGGTGCTGTATGAGTCGTGTATTATGACTATGCGTCTGAGAGTCTGAGGTTCA"
    with open(options.output, 'w') as outfile:
        i = 0
        for header, data in iter_input(options.input):
            i += 1
            seq, boxes = data
            print seq, boxes
            Cbox = re.search("C+", boxes)
            Dbox = re.search("D+", boxes)
            Cpri = re.search("c+", boxes)
            Dpri = re.search("d+", boxes)
            try:
                annotation = "%s_%i_%s_%i_%s_%i_%s_%i" % (seq[Cpri.start(): Cpri.end()],
                                                          Cpri.start() + 1,
                                                          seq[Dpri.start(): Dpri.end()],
                                                          Dpri.start() + 1,
                                                          seq[Cbox.start(): Cbox.end()],
                                                          Cbox.start() + 1,
                                                          seq[Dbox.start(): Dbox.end()],
                                                          Dbox.start() + 1)
            except AttributeError:
                annotation = "_%s_%i_%s_%i" % (seq[Cbox.start(): Cbox.end()],
                                               Cbox.start() + 1,
                                               seq[Dbox.start(): Dbox.end()],
                                               Dbox.start() + 1)

            name, chrom, strand, start, end, win1, win2 = header.split("|")
            name = name.split(":")[-1]

            gtfline = ('{chrom}\tENSEMBL\tCD-snoRNA\t{start}\t{end}\t.\t{strand}\t.\t'
                       '"snoBoard:{name},family:SNORD116,boxes:{boxes},seq:{seq}"\n').format(chrom=chrom,
                                                                                                  start=start,
                                                                                                  end=end,
                                                                                             strand=strand,
                                                                                                  name=name,
                                                                                                  boxes=annotation,
                                                                                                  seq=seq)
            outfile.write(gtfline)
            # with open("./Annotated/" + nhead + ".fa", 'w') as individual:
            #     individual.write(">%s_%s_%s\n%s\n" % (nhead[:5], nhead[5:], annotation, seq))


def iter_input(path_to_input):
    """Iterate input

    Args:
        path_to_input (@todo): @todo

    Returns: @todo

    """
    with open(path_to_input) as infile:
        sequences = (x[1] for x in groupby(infile, lambda line: line[0]==">"))
        for header in sequences:
            fasta_header = header.next()[1:].strip()
            data = [s.strip() for s in sequences.next()]
            yield fasta_header, data


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
