#!/usr/bin/env python
"""
Cluster results according to the position of the hit and miRNA
The input file looks like that:

chr6    99846856    99846871    2628039_1-Unique-1:hsa-miR-129-3p:8   30  -
chr3    30733346    30733368    2630171_1-Unique-1:hsa-miR-93:N       36  +
chr17   3627403     3627417     2632714_1-Unique-1:hsa-miR-186:N      28  +
chr17   3627403     3627417     2639898_1-Unique-1:hsa-miR-16:N       28  +
"""

__date_ = "2014-08-11"
__author__ = "Rafal Gumienny"
__email__ = "r.gumienny@unibas.ch"
__license__ = "GPL"

# imports
import sys
import csv
import time
import numpy as np
import pandas as pd
from operator import itemgetter
from argparse import ArgumentParser
from collections import defaultdict
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
                    help="Input table file in bed like format")
parser.add_argument("--output",
                    dest="output",
                    default="output.tab",
                    help="Output table , defaults to output.tab")
parser.add_argument("--cluster-size",
                    dest="cluster_size",
                    type=int,
                    default=1,
                    help="Number of reads necessary for a group to be considered a cluster. "
                     "eg. 2 returns all groups with 2 or more overlapping reads, defaults to 1")
parser.add_argument("--overlap",
                    dest="overlap",
                    type=int,
                    default=-40,
                    help="Distance in basepairs for two reads to be in the same cluster. "
                     "For instance 20 would group all reads with 20bp of each other. Negative "
                     "number means overlap eg. -10 - read must overlap at leas 10 basepairs, defaults to -1")


# redefine a functions for writing to stdout and stderr to save some writting
syserr = sys.stderr.write
sysout = sys.stdout.write


def main(options):
    """Main logic of the script"""
    cluster_trees, read_map = build_cluster_trees(parse_bed_input(options.input),
                                                  options.overlap,
                                                  options.cluster_size)
    i = 0
    single_clusters = 0
    if options.verbose:
        syserr("Saving output\n")
    with open(options.output, 'w') as o:
        for contig, cluster_tree in cluster_trees.items():
            for start, end, read_ids in cluster_tree.getregions():
                chrom, strand, snoRNA = contig.split(":")
                row = read_map[read_ids[0]] # just take one representative
                o.write("%s\t%i\t%i\t%s\t%s\t%s\t%s\n" % (chrom,
                                                          start,
                                                          end,
                                                          snoRNA,
                                                          len(read_ids),
                                                          strand,
                                                          "\t".join(row[8:]),
                                                          ))


def parse_bed_input(bed_path):
    """
    """
    with open(bed_path) as bed:
        for row in csv.reader(bed, delimiter='\t'):
            # returns: read_id, chrom, strand, start, end
            yield row


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
    cluster_trees = defaultdict(lambda: ClusterTree(cluser_distance, read_count))
    i = 0
    read_ids_mapping = {}
    for row in bed_generator:
        cluster_trees["%s:%s:%s" % (row[0], row[5], row[6])].insert(int(row[1]), int(row[2]), i)
        read_ids_mapping[i] = row
        i += 1
    return cluster_trees, read_ids_mapping

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
