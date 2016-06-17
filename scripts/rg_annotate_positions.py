#!/usr/bin/env python
"""
Annotate found snoRNA target sites
"""

__date__ = "2015-02-24"
__author__ = "Rafal Gumienny"
__email__ = "r.gumienny@unibas.ch"
__license__ = "GPL"

# imports
import sys
import csv
import time
import HTSeq
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
                    help="Output file in tab format.")
parser.add_argument("--regions",
                    dest="regions",
                    required=True,
                    help="GFF file with annotations for different gene regions eg. UTRs")
parser.add_argument("--genes",
                    dest="genes",
                    required=True,
                    help="Positions of all genes in GFF format")
parser.add_argument("--snoRNAs",
                    dest="snoRNAs",
                    help="GFF file with annotations for snoRNAs in the same format as genes file")
parser.add_argument("--repeats",
                    dest="repeats",
                    required=True,
                    help="GTF file with annotations for repeats in the format from rmsk table in UCSC")

try:
    options = parser.parse_args()
except Exception, e:
    parser.print_help()

# redefine a functions for writing to stdout and stderr to save some writting
syserr = sys.stderr.write
sysout = sys.stdout.write

class GeneInfo:

    def __init__(self, gene_id, gene_name, gene_type):
        self.gene_id = gene_id
        self.gene_name = gene_name
        self.gene_type = gene_type

    def __repr__(self):
        return "%s:%s" % (self.gene_type,
                          self.gene_id)


def main():
    """Main logic of the script"""
    genes = HTSeq.GenomicArrayOfSets("auto", stranded=True)
    intron_exon = HTSeq.GenomicArrayOfSets("auto", stranded=True)
    mRNA_part = HTSeq.GenomicArrayOfSets("auto", stranded=True)
    if options.verbose:
        syserr("Reading annotations\n")

    with open(options.genes) as g_inf:
        for rec in HTSeq.GFF_Reader(g_inf):
            try:
                genes[rec.iv] += GeneInfo(rec.name,
                                         rec.attr["gene_name"],
                                         rec.attr["gene_biotype"])
            except KeyError:
                genes[rec.iv] += GeneInfo(rec.name,
                                         rec.attr["gene_id"],
                                         rec.attr["gene_biotype"])
    if options.snoRNAs:
        with open(options.snoRNAs) as g_inf:
            for rec in HTSeq.GFF_Reader(g_inf):
                try:
                    genes[rec.iv] += GeneInfo(rec.name,
                                             rec.attr["gene_name"],
                                             rec.attr["gene_biotype"])
                except KeyError:
                    genes[rec.iv] += GeneInfo(rec.name,
                                             rec.attr["gene_id"],
                                             rec.attr["gene_biotype"])
    if options.repeats:
        with open(options.repeats) as g_inf:
            for rec in HTSeq.GFF_Reader(g_inf):
                rep_type = 'repeat'
                rec.iv.chrom = rec.iv.chrom[3:] if rec.iv.chrom.startswith("chr") else rec.iv.chrom
                genes[rec.iv] += GeneInfo(rec.attr['gene_id'],
                                         rec.attr["transcript_id"],
                                         rep_type)

    id_to_name = {}
    for rec in HTSeq.GFF_Reader(options.regions):
        if rec.type == "gene":
            id_to_name[rec.type] = rec.attr["Name"]
        elif rec.type == 'intron' or rec.type == 'exon':
            intron_exon[rec.iv] += rec.type
        elif rec.type == "CDS" or rec.type == "three_prime_UTR" or rec.type == "five_prime_UTR":
            mRNA_part[rec.iv] += rec.type
        else:
            pass
    if options.verbose:
        syserr("Annotating file\n")

    with open(options.input) as infile, open(options.output, "w") as outfile:
        for row in csv.reader(infile, delimiter="\t"):
            chrom = row[0][3:] if row[0].startswith("chr") else row[0]
            iv = HTSeq.GenomicInterval(chrom, int(row[1]), int(row[2]), row[5])
            annot_gene_types = ",".join(set([i.gene_type for g in genes[iv] for i in g]))
            annot_gene_ids = ",".join(set([i.gene_id for g in genes[iv] for i in g]))
            annot_gene_names = ",".join(set([i.gene_name for g in genes[iv] for i in g]))
            annot_int_ex = ",".join(set([str(i) for g in intron_exon[iv] for i in g]))
            annot_mRNA_part = ",".join(set([str(i) for g in mRNA_part[iv] for i in g]))
            outfile.write("%s\t%s\t%s\t%s\t%s\t%s\n" % ("\t".join(row),
                                                        annot_gene_types if annot_gene_types else 'NA',
                                                        annot_gene_ids if annot_gene_ids else 'NA',
                                                        annot_gene_names if annot_gene_names else 'NA',
                                                        annot_int_ex if annot_int_ex else 'NA',
                                                        annot_mRNA_part if annot_mRNA_part else 'NA'))

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
