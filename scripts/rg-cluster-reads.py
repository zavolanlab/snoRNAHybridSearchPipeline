#!/usr/bin/env python
"""
Cluster reads into more convinient bed file
"""

__date_ = "2014-06-30"
__author__ = "Rafal Gumienny"
__email__ = "r.gumienny@unibas.ch"
__license__ = "GPL"

# imports
import csv
import sys
import time
import collections
import pylab as pl
from operator import itemgetter
from optparse import OptionParser
from argparse import ArgumentParser
from Bio import SeqIO
from bx.intervals.cluster import ClusterTree

parser = ArgumentParser(description=__doc__)
parser.add_argument("-v",
                    "--verbose",
                    dest="verbose",
                    action="store_true",
                    default=False,
                    help="Be loud!")
parser.add_argument("--input",
                    dest="input",
                    required=True,
                    help="Input file in special asmbed format or in bed format")
parser.add_argument("--bed",
                    dest="bed",
                    action="store_true",
                    default=False,
                    help="Specifies if the input file is in bed format")
parser.add_argument("--cluster-size",
                    dest="cluster_size",
                    type=int,
                    default=1,
                    help="Number of reads necessary for a group to be considered a cluster. "
                     "eg. 2 returns all groups with 2 or more overlapping reads, defaults to 1")
parser.add_argument("--overlap",
                    dest="overlap",
                    type=int,
                    default=-1,
                    help="Distance in basepairs for two reads to be in the same cluster. "
                     "For instance 20 would group all reads with 20bp of each other. Negative "
                     "number means overlap eg. -10 - read must overlap at leas 10 basepairs, defaults to -1")
parser.add_argument("--expand-cluster",
                    dest="expand_cluster",
                    type=int,
                    default=0,
                    help="Expand cluster in both directions, defaults to 0")
parser.add_argument("--expand-read",
                    dest="expand_read",
                    type=int,
                    default=15,
                    help="Expand read in both directions (some alternative to expand cluseter), defaults to 15")
parser.add_argument("--output",
                    dest="output",
                    default="output.bed",
                    help="Output file in bed format , defaults to output.bed")
parser.add_argument("--asmbed",
                    dest="asmbed",
                    action="store_true",
                    default=False,
                    help="Write in asmbed format for fasta extraction")
parser.add_argument("--rRNAs",
                    dest="rRNAs",
                    required=False,
                    help="rRNAs to add in the end of the clusters")
parser.add_argument("--tRNAs",
                    dest="tRNAs",
                    required=False,
                    help="tRNAs to add in the end of the clusters")
parser.add_argument("--snRNAs",
                    dest="snRNAs",
                    required=False,
                    help="snRNAs to add in the end of the clusters")
filters = parser.add_mutually_exclusive_group(required=False)
filters.add_argument("--filter-by",
                     dest="filter_by",
                     help="Keep only read with these tags in read_ids.\nInput is coma separated list of tags")
filters.add_argument("--filter-except",
                     dest="filter_except",
                     help="Keep read except with these tags in read_ids.\nInput is coma separated list of tags")

try:
    options = parser.parse_args()
except Exception, e:
    parser.print_help()

# redefine a functions for writing to stdout and stderr to save some writting
syserr = sys.stderr.write
sysout = sys.stdout.write


def main():
    """Main logic of the script"""
    if not options.bed:
        cluster_trees, read_map = build_cluster_trees(parse_input(options.input),
                                                        options.overlap,
                                                        options.cluster_size)
    else:
        cluster_trees, read_map = build_cluster_trees(parse_bed_input(options.input, options.filter_by, options.filter_except, expand_read=options.expand_read),
                                                            options.overlap,
                                                            options.cluster_size)
    i = 0
    single_clusters = 0
    if options.verbose:
        syserr("Saving output\n")
    with open(options.output, 'w') as o:
        for contig, cluster_tree in cluster_trees.items():
            for start, end, read_ids in cluster_tree.getregions():
                cluster_name = "Cluster"
                if len(read_ids) == 1:
                    single_clusters += 1
                i += 1
                chrom, strand = contig.split(":")
                if options.asmbed:
                    o.write("%s\t%s\t%s\t%i\t%i\n" %("%s_%i_%i_%i" % (cluster_name, i, end - start, len(read_ids)),
                                                     chrom,
                                                     strand,
                                                     start - options.expand_cluster + 1,  # asmbed is 1-based format
                                                     end + options.expand_cluster))
                else:
                    o.write("%s\t%i\t%i\t%s\t%i\t%s\n" %(chrom,
                                                         start - options.expand_cluster,
                                                         end + options.expand_cluster,
                                                         "%s_%i_%i" % (cluster_name, i, end - start),
                                                         len(read_ids),
                                                         strand))
        if options.rRNAs:
            with open(options.rRNAs) as f:
                for rec in SeqIO.parse(f, 'fasta'):
                    i += 1
                    if options.asmbed:
                        o.write("Cluster_%i_%i_%i\t%s\t%s\t%i\t%i\n" % (i,
                                                                    len(rec.seq),
                                                                    1,
                                                                    str(rec.id).split("|")[0],
                                                                    "+",
                                                                    1,
                                                                    len(rec.seq)))
                    else:
                        o.write("%s\t%i\t%i\t%s\t%i\t%s\n" %(str(rec.id).split("|")[0],
                                                             0,
                                                             len(rec.seq),
                                                             "Cluster_%s_%i" % (i, len(rec.seq)),
                                                             1,
                                                             "+"))
        if options.tRNAs:
            with open(options.tRNAs) as f:
                for rec in SeqIO.parse(f, 'fasta'):
                    i += 1
                    if options.asmbed:
                        o.write("Cluster_%i_%i_%i\t%s\t%s\t%i\t%i\n" % (i,
                                                                    len(rec.seq),
                                                                    1,
                                                                    str(rec.id).split("|")[0],
                                                                    "+",
                                                                    1,
                                                                    len(rec.seq)))
                    else:
                        o.write("%s\t%i\t%i\t%s\t%i\t%s\n" %(str(rec.id).split("|")[0],
                                                             0,
                                                             len(rec.seq),
                                                             "Cluster_%s_%i" % (i, len(rec.seq)),
                                                             1,
                                                             "+"))
        if options.snRNAs:
            with open(options.snRNAs) as f:
                for rec in SeqIO.parse(f, 'fasta'):
                    i += 1
                    if options.asmbed:
                        o.write("Cluster_%i_%i_%i\t%s\t%s\t%i\t%i\n" % (i,
                                                                    len(rec.seq),
                                                                    1,
                                                                    str(rec.id).split("|")[0],
                                                                    "+",
                                                                    1,
                                                                    len(rec.seq)))
                    else:
                        o.write("%s\t%i\t%i\t%s\t%i\t%s\n" %(str(rec.id).split("|")[0],
                                                             0,
                                                             len(rec.seq),
                                                             "Cluster_%s_%i" % (i, len(rec.seq)),
                                                             1,
                                                             "+"))
    if options.verbose:
        syserr("%i clusters were written (%i of them are single-read clusters)\n" % (i, single_clusters))


def parse_input(bed_path):
    """
    Parse asmbed file and yield read_id, chrom, strand, start, end
    Eg. input:
    seq33   chr2    +   26676380    26676398    1   mRNA
    seq46   chr5    +   134193246   134193263   2   mRNA
    seq56   chr19   +   37468373    37468391    1   mRNA
    seq57   chr19   +   30281781    30281798    1   repeat
    """
    with open(bed_path) as bed:
        for row in csv.reader(bed, delimiter='\t'):
            # returns: read_id, chrom, strand, start, end
            yield (row[0],
                   row[1],
                   row[2],
                   int(row[3]) - 1,  # input is 1-based
                   int(row[4]))

def parse_bed_input(bed_path, filter_by=None, filter_except=None, expand_read=15):
    """
    Parse bed file and yield read_id, chrom, strand, start, end
    Eg. input:
    seq33   chr2    +   26676380    26676398    1   mRNA
    seq46   chr5    +   134193246   134193263   2   mRNA
    seq56   chr19   +   37468373    37468391    1   mRNA
    seq57   chr19   +   30281781    30281798    1   repeat
    """
    if not filter_by and not filter_except:
        with open(bed_path) as bed:
            for row in csv.reader(bed, delimiter='\t'):
                # returns: read_id, chrom, strand, start, end
                yield (row[3],
                       row[0],
                       row[5],
                       int(row[1]) - expand_read,
                       int(row[2]) + expand_read)
    elif filter_by:
        with open(bed_path) as bed:
            filter_by_list = filter_by.split(",")
            for row in csv.reader(bed, delimiter='\t'):
                # returns: read_id, chrom, strand, start, end
                if any([i in row[3] for i in filter_by_list]):
                    yield (row[3],
                           row[0],
                           row[5],
                           int(row[1]) - expand_read,
                           int(row[2]) + expand_read)
    elif filter_except:
        with open(bed_path) as bed:
            filter_except_list = filter_except.split(",")
            for row in csv.reader(bed, delimiter='\t'):
                # returns: read_id, chrom, strand, start, end
                if not any([i in row[3] for i in filter_except_list]):
                    yield (row[3],
                           row[0],
                           row[5],
                           int(row[1]) - expand_read,
                           int(row[2]) + expand_read)


def build_cluster_trees(bed_generator, cluser_distance, read_count):
    """
    arguments to ClusterTree are:
    - Distance in basepairs for two reads to be in the same cluster;
      for instance 20 would group all reads with 20bp of each other
    - Number of reads necessary for a group to be considered a cluster;
      2 returns all groups with 2 or more overlapping reads
    """
    if options.verbose:
        syserr("Making ClusterTree\n")
    cluster_trees = collections.defaultdict(lambda: ClusterTree(cluser_distance, read_count))
    i = 0
    read_ids_mapping = {}
    for read_id, match_id, strand, start, end in bed_generator:
        cluster_trees["%s:%s" % (match_id, strand)].insert(start, end, i)
        read_ids_mapping[i] = read_id
        i += 1
    return cluster_trees, read_ids_mapping


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
