#!/usr/bin/env python
"""
Take SuppD1 and SuppD5 and compile snoRNA table for pipeline input. The output columns are:

1. Chromosome   Chromosome on which gene is located
2. Start    0-based start coordinate
3. End  1-based end coordinate
4. snoRNA ID    ID for snoRNA
5. Type General type of snoRNA (C/D or H/ACA)
6. Strand   Strand on which gene is located
7. Sequence Sequence of snoRNA
8. D-boxes  1-based coordinates snoRNA D-boxes in form of D-box_start..D-box_end or D'-box_start..D'-box_end,D-box_start..D-box_end eg. 60..63 or 30..33,60..63. It must be 4 nucleotides long.
9. C-boxes  Same as above (not that prime boxes are reversed in positions) with exception of the length constraint
10. H-box   1-based position of the H box in H/ACA box snoRNAs in form of start..end
11. ACA-box 1-based position of the ACA box in H/ACA box snoRNA in form of start..end
12. Alias   Alias/Alternative name for snoRNA
13. Gene name   Name of the snoRNA gene
14. Accession   Accession (NCBI, RefSeq etc) for snoRNA
15. Modification sites  Coma separated list of modified sites in form of ModifiedRNA:NucleotidePosition eg. 28S:U46,18S:G52
16. Host gene   Host gene name if any
17. Host ID Host locus ID
18. Organization    Organisation of the gene eg. intronic
19. Organism    Species
20. Note    Additional info about snoRNA
21. snoRNA name    The official name for snoRNA (if available)
22. snoRNA family    The family of snoRNA (if available)
23. snoRNA precise type    Precise type of the snoRNA as opposed to general type
"""

__date__ = "2015-10-31"
__author__ = "Rafal Gumienny"
__email__ = "r.gumienny@unibas.ch"
__license__ = "GPL"

# imports
import re
import sys
import csv
import time
from contextlib import contextmanager
import pandas as pd
from argparse import ArgumentParser, RawTextHelpFormatter

parser = ArgumentParser(description=__doc__, formatter_class=RawTextHelpFormatter)
parser.add_argument("-v",
                    "--verbose",
                    dest="verbose",
                    action="store_true",
                    default=False,
                    help="Be loud!")
parser.add_argument("--supp-D5",
                    dest="supp_D5",
                    required=True,
                    help="Supplementary data 5 with snoRNAs as GFF file from atlas paper")
parser.add_argument("--organism",
                    dest="organism",
                    required=True,
                    help="Organism of the snoRNAs")
parser.add_argument("--output",
                    dest="output",
                    default="snoRNAs_table.tab",
                    help="The name of the output table , defaults to snoRNAs_table.tab")


try:
    options = parser.parse_args()
except Exception, e:
    parser.print_help()

# redefine a functions for writing to stdout and stderr to save some writting
syserr = sys.stderr.write
sysout = sys.stdout.write


def main():
    """Main logic of the script"""
    snor_type_map = {"ALUACA-snoRNA": "HACA",
                     "CD-SCARNA-snoRNA": "CD",
                     "CD-snoRNA": "CD",
                     "HACA-SCARNA-snoRNA": "HACA",
                     "HACA-snoRNA": "HACA",
                     "Hybrid-snoRNA": "Hybrid",
                     "sno-lncRNA-snoRNA": "Hybrid",
                     "SNORD-like-snoRNA": "CD",
                     "Tandem-CD-snoRNA": "CD",
                     "Tandem-HACA-snoRNA": "HACA"
                     }
    # families, modification_sites = parse_supp_D1(options.supp_D1)
    with open(options.output, 'w') as outfile:
        for (chrom, start, end, strand, snor_type, source, attrs) in iter_snoRNAs(options.supp_D5):
            if attrs['boxes'] == 'NA':
                dboxes = 'NA'
                cboxes = 'NA'
                hbox = 'NA'
                acabox = 'NA'
            elif snor_type_map[snor_type] == 'CD':
                try:
                    dboxes, cboxes = get_box_string(attrs['boxes'], 'CD')
                    hbox = 'NA'
                    acabox = 'NA'
                except ValueError:
                    print "Line with ID %s is wrong." % attrs["ID"]
            elif snor_type_map[snor_type] == 'HACA':
                try:
                    hbox, acabox = get_box_string(attrs['boxes'], 'HACA')
                    dboxes = 'NA'
                    cboxes = 'NA'
                except ValueError:
                    print "Line with ID %s is wrong." % attrs["ID"]
            else:
                print "Unknown type: %s" % snor_type_map[snor_type]
            output_tuple = (chrom,
                            start,
                            end,
                            attrs['snoBoard'],
                            snor_type_map[snor_type],
                            strand, attrs['seq'],
                            dboxes,
                            cboxes,
                            hbox,
                            acabox,
                            'NA',
                            'NA',
                            'NA',
                            'NA',
                            'NA',
                            'NA',
                            'NA',
                            options.organism,
                            'NA',
                            'NA',
                            attrs['family'],
                            snor_type)
            outfile.write("\t".join(map(str, output_tuple)) + "\n")


def parse_supp_D1(path):
    df = pd.read_csv(path)
    #df[['ID', 'Name', 'RNA family name', 'Reported1 (R) ', 'Reported2 (R)']]
    df.index = df['ID']
    mods = {}
    for snor_id, row in df.iterrows():
        # print snor_id
        tmp_mods = []
        rep1 = row['Reported1 (R) ']
        rep2 = row['Reported2 (R)']
        if not pd.isnull(rep1):
            if rep1 != '-':
                for m in rep1.split(";"):
                    if len(m.split("-")) == 2:
                        tmp_mods.append(m.split("-")[0] + ":N" + m.split("-")[1])
                    elif len(m.split("-")) == 3:
                        # this is for situations like this: 5-8S-75
                        tmp_mods.append(".".join(m.split("-")[:2]) + ":N" + m.split("-")[2])
                    else:
                        syserr("The modifications site is strange for %s: %s\n" % (snor_id, m))
        if not pd.isnull(rep2):
            if rep2 != '-':
                for m in rep2.split(";"):
                    if len(m.split("-")) == 2:
                        tmp_mods.append(m.split("-")[0] + ":N" + m.split("-")[1])
                    elif len(m.split("-")) == 3:
                        # this is for situations like this: 5-8S-75
                        tmp_mods.append(".".join(m.split("-")[:2]) + ":N" + m.split("-")[2])
                    else:
                        syserr("The modifications site is strange for %s: %s\n" % (snor_id, m))
        if len(tmp_mods) > 0:
            mods[snor_id] = ",".join(tmp_mods)
        else:
            mods[snor_id] = 'NA'

    fams = {}
    for snor_id, row in df.iterrows():
        if not pd.isnull(row['RNA family name']):
            fams[snor_id] = row['RNA family name']
        else:
            fams[snor_id] = 'NA'

    return fams, mods


def iter_snoRNAs(input_file):
    with open(input_file) as infile:
        for row in csv.reader(infile, delimiter='\t'):
            chrom = row[0] if not row[0].startswith('chr') else row[0][3:]
            source = row[1]
            snor_type = row[2]
            start = int(row[3]) - 1
            end = int(row[4])
            strand = row[6]
            attrs = {atr.split(":")[0]: atr.split(":")[1] for atr in row[8].replace('"','').split(",")}
            yield (chrom, start, end, strand, snor_type, source, attrs)


def get_box_string(box, snor_type):
    if snor_type == 'CD':
        boxes = re.split("[_]+", box)
        if len(boxes) == 4:
            dbox = "%s..%s" % (boxes[-1], int(boxes[-1]) + 3)
            if "NNNN" in boxes[-4]:
                cbox = 'NA'
            else:
                cbox = "%s..%s" % (boxes[-3], int(boxes[-3]) + 5)
            return dbox, cbox
        elif len(boxes) == 8:
            dbox = "%s..%s" % (boxes[-1], int(boxes[-1]) + 3)

            if "NNNN" in boxes[-4]:
                cbox = 'NA'
            else:
                cbox = "%s..%s" % (boxes[-3], int(boxes[-3]) + 5)

            if "NNNN" in boxes[-6]:
                dpbox = 'NA'
            else:
                dpbox = "%s..%s" % (boxes[-5], int(boxes[-5]) + 3)

            if "NNNN" in boxes[-8]:
                cpbox = 'NA'
            else:
                cpbox = "%s..%s" % (boxes[-7], int(boxes[-7]) + 5)

            if dpbox == 'NA':
                dboxes = "%s" % (dbox)
            else:
                dboxes = "%s,%s" % (dpbox, dbox)

            if cpbox == 'NA':
                cboxes = "%s" % (cbox)
            else:
                cboxes = "%s,%s" % (cpbox, cbox)

            return dboxes, cboxes

    elif snor_type == 'HACA':
        boxes = re.split("[_]+", box)
        acabox = '%s..%i' % (boxes[-1], int(boxes[-1]) + 2)
        hbox = '%s..%i' % (boxes[-3], int(boxes[-3]) + 5)
        return hbox, acabox
    else:
        print "Unknown snoRNA type %s" % snor_type
        return 'NA', 'NA'

# this function is also defined in utils but I put it here to avoid
# unnecessary module import that might not be available everywhere as
# it is my own module
@contextmanager
def smart_open(filepath, mode='r'):
    """Open file intelligently depending on the source

    :param filepath: can be both path to file or sys.stdin or sys.stdout
    :param mode: mode can be read "r" or write "w". Defaults to "r"
    :yield: context manager for file handle

    """
    if mode == 'r':
        if filepath is not sys.stdin:
            fh = open(filepath, 'r')
        else:
            fh = filepath
        try:
            yield fh
        except IOError as e:
            if fh is not sys.stdin:
                fh.close()
            elif e.errno == errno.EPIPE:
                pass
        finally:
            if fh is not sys.stdin:
                fh.close()
    elif mode == 'w':
        if filepath is not sys.stdout:
            fh = open(filepath, 'w')
        else:
            fh = filepath
        try:
            yield fh
        finally:
            if fh is not sys.stdout:
                fh.close()
    else:
        raise NoSuchModeException("No mode %s for file" % mode)


if __name__ == '__main__':
    try:
        if options.verbose:
            start_time = time.time()
            start_date = time.strftime("%d-%m-%Y at %H:%M:%S")
            syserr("############## Started script on %s ##############\n" %
                   start_date)
        main()
        if options.verbose:
            syserr("### Successfully finished in %i seconds, on %s ###\n" %
                   (time.time() - start_time,
                    time.strftime("%d-%m-%Y at %H:%M:%S")))
    except KeyboardInterrupt:
        syserr("Interrupted by user after %i seconds!\n" %
               (time.time() - start_time))
        sys.exit(-1)
