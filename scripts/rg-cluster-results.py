#!/usr/bin/env python
"""
Cluster results according to the position of the hit
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
                    default=-1,
                    help="Distance in basepairs for two reads to be in the same cluster. "
                     "For instance 20 would group all reads with 20bp of each other. Negative "
                     "number means overlap eg. -10 - read must overlap at leas 10 basepairs, defaults to -1")

try:
    options = parser.parse_args()
except Exception, e:
    parser.print_help()

# redefine a functions for writing to stdout and stderr to save some writting
syserr = sys.stderr.write
sysout = sys.stdout.write


def main():
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
                chrom, strand, mirna = contig.split(":")
                o.write("%s\t%i\t%i\t%s\t%i\t%s\n" %(chrom,
                                                     start,
                                                     end,
                                                     mirna,
                                                     len(read_ids),
                                                     strand))


def parse_bed_input(bed_path):
    """
    """
    with open(bed_path) as bed:
        for row in csv.reader(bed, delimiter='\t'):
            # returns: read_id, chrom, strand, start, end
            yield (row[3],
                   row[0],
                   row[5],
                   int(row[1]),  # input is 1-based
                   int(row[2]),
                   row[6])


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
    for read_id, match_id, strand, start, end, mirna in bed_generator:
        cluster_trees["%s:%s:%s" % (match_id, strand, mirna)].insert(start, end, i)
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
