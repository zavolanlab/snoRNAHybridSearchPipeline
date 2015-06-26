#!/usr/bin/env python
"""
For each read in the file check if there is an anchor sequence
and if this is the case make local alignment (SW) for each associated
sequence. As a sequence in the read take only that with the best score.
"""

__date_ = "2014-06-05"
__author__ = "Rafal Gumienny"
__email__ = "r.gumienny@unibas.ch"
__license__ = "GPL"

# imports
import sys
import csv
import time
import swalign
from ushuffle import shuffle
from Bio import SeqIO
from collections import OrderedDict, defaultdict
from argparse import ArgumentParser
parser = ArgumentParser(description=__doc__)
parser.add_argument("-v",
                    "--verbose",
                    dest="verbose",
                    action="store_true",
                    default=False,
                    help="Be loud!")
parser.add_argument("--anchors",
                    dest="anchors",
                    help="File with anchors (tab-separated)")
parser.add_argument("--anchor-sequences",
                    dest="anchor_sequences",
                    help="Sequences from which anchors were generated")
parser.add_argument("--reads",
                    dest="reads",
                    help="File with reads")
parser.add_argument("--match",
                    dest="match",
                    type=int,
                    default=2,
                    help="Match score, defaults to 2")
parser.add_argument("--mismatch",
                    dest="mismatch",
                    type=int,
                    default=-5,
                    help="Mismatch penalty, defaults to -5")
parser.add_argument("--gap-open",
                    dest="gap_open",
                    type=int,
                    default=-6,
                    help="Open gap penalty, defaults to -6")
parser.add_argument("--gap-extend",
                    dest="gap_extend",
                    type=int,
                    default=-4,
                    help="Gap extension penalty, defaults to -4")
parser.add_argument("--output",
                    dest="output",
                    help="Output table")
parser.add_argument("--RNase-T1",
                    dest="RNase_T1",
                    action="store_true",
                    default=False,
                    help="Indicates if in the experiment RNase T1 was used")


try:
    options = parser.parse_args()
except Exception, e:
    parser.print_help()

# redefine a functions for writing to stdout and stderr to save some writting
syserr = sys.stderr.write
sysout = sys.stdout.write


def main():
    """Main logic of the script"""
    anchors = read_anchors(options.anchors)
    anchor_seqs = read_fasta_to_dict(options.anchor_sequences)
    all_reads = read_fasta_to_dict_and_shuffle(options.reads)
    #
    # prepare swalign parameters
    #
    scoring = swalign.NucleotideScoringMatrix(match=options.match,
                                              mismatch=options.mismatch)
    sw = swalign.LocalAlignment(scoring,
                                gap_penalty=options.gap_open,
                                gap_extension_penalty=options.gap_extend)

    #
    # Iterate reads to find matches
    #
    # counter = 0
    with open(options.output, 'w') as out:
        # out.write("#read_id\tname\tscore\tbeg\tend\twhich_end\ttarget\tqbeg\tqend\tisambig\tmatches:mismatches:identity")
        for read_id, reads in all_reads.iteritems():
            read, shuffled_read = reads
            alignments = get_alignments(anchors, read, anchor_seqs, sw)
            shuffled_alignments = get_alignments(anchors, shuffled_read, anchor_seqs, sw)
            if len(alignments.keys()) != 0:
                outtexts = get_best_alignment(alignments, read_id, read)
            else:
                outtexts = ["%s\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA" % (read_id)]
            if len(shuffled_alignments.keys()) != 0:
                shuffled_text = get_best_shuffled_text(shuffled_alignments, read)
                for i, text in enumerate(outtexts):
                    outtexts[i] += shuffled_text
            else:
                for i, text in enumerate(outtexts):
                    outtexts[i] += "\tNA\tNA\tNA\tNA\n"
            for text in outtexts:
                out.write(text)

def get_alignments(anchors, read, anchor_seqs, sw):
    alignments = defaultdict(list)
    for anchor, anchor_names in anchors.iteritems():
        if anchor in read:
            for anchor_name in anchor_names:
                tmp_alignment = sw.align(read, anchor_seqs[anchor_name])
                alignments[tmp_alignment.score].append((anchor_name,
                                                                tmp_alignment))
    return alignments

def get_best_alignment(alignments, read_id, read):
    max_score = max([a for a in alignments])
    my_alignments = alignments[max_score]
    ambiguous = False
    ambiguous_count = len(set([al[0] for al in my_alignments]))
    if ambiguous_count != 1:
        # raise Exception("Ambiguous Alignment!")
        ambiguous = True
    names = []
    texts = []
    for name, alignment in my_alignments:
        if name in names:
            continue
        names.append(name)
        # name, alignment = my_alignments[0]
        score = alignment.score
        beg = alignment.r_pos
        end = alignment.r_end
        five_prime_end = read[:beg]
        three_prime_end = read[end:]
        mirseq = read[beg: end]
        which_end = '3p' if len(five_prime_end) > len(three_prime_end) else '5p'
        if mirseq[-2].upper() == 'G' and options.RNase_T1:
            read_rest = read[:beg] if len(five_prime_end) > len(three_prime_end) else read[end - 1:]
        if mirseq[-3].upper() == 'G' and options.RNase_T1:
            read_rest = read[:beg] if len(five_prime_end) > len(three_prime_end) else read[end - 2:]
        else:
            read_rest = five_prime_end if len(five_prime_end) > len(three_prime_end) else three_prime_end
        #read_id name score beg end which_end target qbeg qend isambig matches:mismatches:identity
        ambiguous_string = "Ambiguous-%i" % ambiguous_count if ambiguous else "Unique-1"
        read_name = read_id + "-" + ambiguous_string + ":" + name
        text = "%s\t%s\t%f\t%i\t%i\t%s\t%s\t%i\t%i\t%s\t%s" % (read_name,
                                                               name,
                                                               score,
                                                               beg,
                                                               end,
                                                               which_end,
                                                               read_rest,
                                                               alignment.q_pos,
                                                               alignment.q_end,
                                                               "Ambiguous" if ambiguous else "Unique",
                                                               "%i:%i:%.4f" % (alignment.matches,
                                                                               alignment.mismatches,
                                                                               alignment.identity))
        texts.append(text)
    return texts

def get_best_shuffled_text(alignments, read):
    max_score = max([a for a in alignments])
    my_alignments = alignments[max_score]
    name, alignment = my_alignments[0]
    five_prime_end = read[:alignment.r_pos]
    three_prime_end = read[alignment.r_end:]
    read_rest = five_prime_end if len(five_prime_end) > len(three_prime_end) else three_prime_end
    #shuf_score shuf_target_len shuf_qbeg shuf_qend
    return "\t%s\t%i\t%i\t%i\t%i\n" % (name,
                                       alignment.score,
                                       len(read_rest),
                                       alignment.q_pos,
                                       alignment.q_end)

def read_anchors(path):
    """Read anchors into dictionary

    Args:
        path (str): path to file

    Returns: dictionary

    """
    out_dict = {}
    with open(path) as f:
        for row in csv.reader(f, delimiter="\t"):
            out_dict[row[0]] = row[1].split(",")
    return out_dict


def read_fasta_to_dict_and_shuffle(path_to_file):
    """Read fasta file into dictionary

    Args:
        path_to_file (str): path to file

    Returns: dictionary

    """
    if options.verbose:
        syserr("Reading sequences from %s \n" % (path_to_file))
    try:
        seq_obj = open(path_to_file, 'Ur')
        seqs = OrderedDict()
        for seq in SeqIO.parse(seq_obj, 'fasta'):
            read = str(seq.seq).upper().replace("U", "T")
            seqs[str(seq.id)] = (read, shuffle(read, len(read), 2)) # second read is the shuffled version of the first
    except IOError:
        raise IOError('Cannot read from %s' % (path_to_file))

    return seqs

def read_fasta_to_dict(path_to_file):
    """Read fasta file into dictionary

    Args:
        path_to_file (str): path to file

    Returns: dictionary

    """
    if options.verbose:
        syserr("Reading sequences from %s \n" % (path_to_file))
    try:
        seq_obj = open(path_to_file, 'Ur')
        seqs = OrderedDict()
        for seq in SeqIO.parse(seq_obj, 'fasta'):
            read = str(seq.seq).upper().replace("U", "T")
            seqs[str(seq.id)] = read
    except IOError:
        raise IOError('Cannot read from %s' % (path_to_file))

    return seqs

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
