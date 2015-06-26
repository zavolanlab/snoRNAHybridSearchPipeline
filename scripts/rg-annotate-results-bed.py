"""
rg-annotate-bed.py
@Author:      Rafal Gumienny (gumiennr@unibas.ch)
@Created:     12-Dec-12
@Description: Annotate bed file with another bed file containing annotations
@Usage:       python rg-annotate-bed.py -h
"""


import sys
import csv
import StringIO
#import logging
#import logging.handlers
from bx.intervals.intersection import Intersecter, Interval
from optparse import OptionParser


####################################
#   DEFINE OPTIONS OF THE SCRIPT   #
####################################

class MyParser(OptionParser):
    """Change parsing of the epilog"""
    def format_epilog(self, formatter):
        return self.epilog

# you might want to add eppilog. It will be formatted as you write it!
epilog = """
        ########################## FILE DESCRIPTION ###################################################

        BED FILE FOR WITH ANNOTATION EXAMPLE
            1   24740163    24740215    miRNA:ENST00000003583 0   -
            1   24727808    24727946    miRNA:ENST00000003583 0   -
            1   24710391    24710493    miRNA:ENST00000003583 0   -

        fields: chr start end annot_type:annot_name num strand"]

        INPUT BED FILE EXAMPLE
            1   24685109    24687340    ENST00000003583 0   -
            1   24687531    24696163    ENST00000003583 0   -
            1   24696329    24700191    ENST00000003583 0   -


        ########################## FILE DESCRIPTION ###################################################
        """
parser = MyParser(version="%prog 1.0", usage="\n\n    %progs [options]", epilog=epilog)
#
#   general options: verbosity / logging
#
parser.add_option("-v",
                  "--verbose",
                  dest="verbose",
                  action="count",
                  default=0,
                  help="Print more verbose messages for each additional verbose level.")
#
#   script options
#
parser.add_option("--input",
                  dest="input",
                  metavar="FILE",
                  type="string",
                  help="a bed file that you want to annotate")
parser.add_option("--output",
                  dest="output",
                  metavar="FILE",
                  type="string",
                  default="output.tab",
                  help="an output table with annotations")
parser.add_option("--annotations",
                  dest="annotations",
                  metavar="FILE",
                  type="string",
                  help="a bed file with annotations")
parser.add_option("--fraction",
                  dest="fraction",
                  metavar="FLOAT",
                  type=float,
                  default=0.1,
                  help="Fraction of read that must overlap the feature to be accepted")
parser.add_option("--placeholder", dest="placeholder",
                  metavar="STRING",
                  type="string",
                  default=".",
                  help="A placeholder for empty annotations")
parser.add_option("--un_stranded",
                  action="store_true",
                  default=False,
                  help="Pass if your protocol is un-stranded")
parser.add_option("--filter-by",
                  dest="filter_by",
                  help="Filter by these (coma separated) list of annotation types")

# in the end add default values for options
for option in parser.option_list:
    if option.default != ("NO", "DEFAULT"):
        option.help += (" " if option.help else "") + "[default: " + str(option.default) + "]"


####################################
#   FUNCTIONS                      #
####################################


def check_mandatory_options(options, mandatory_options, helpstr):
    """
    Check if specified mandatory options have been defined.
    Parameter mandatory_options is a list of strings that
    correspond to the "dest" parameter in the options.
    """
    missing_options = []
    for o in mandatory_options:
        if not getattr(options, o):
            missing_options.append("--" + o)

    if not len(missing_options):
        return

    raise Exception("Missing mandatory parameter%s: %s.\n\n%s\n\n" %
                    ("s" if len(missing_options) > 1 else "",
                     ", ".join(missing_options),
                     helpstr))


def index_bed(genes, un_stranded):
    """ """
    chroms = set([i['chrom'] for i in genes])
    if not un_stranded:
        bed_quicksect = {chr + ":" + str: Intersecter()
                         for str in ['+', '-']
                         for chr in set(chroms)}
    else:
        bed_quicksect = {chr: Intersecter()
                         for chr in set(chroms)}
    for annotation in genes:
        chrom = annotation['chrom']
        strand = annotation['strand']
        contig = chrom if un_stranded else chrom + ":" + strand
        bed_quicksect[contig].insert_interval(Interval(annotation['beg'],
                                                       annotation['end'],
                                                       chrom=chrom,
                                                       strand=strand,
                                                       value={"annotation_type": annotation['annotation_type'],
                                                              "annotation_name": annotation['annotation_name'],
                                                              }))
    return bed_quicksect


def read_annotation_bed_file(filepath):
    """Read annotation bed file into list of dictionaries"""
    output_list = []
    with open(filepath) as annotations:
        for line in csv.reader(annotations, delimiter="\t"):
            if not line[0].startswith("#"):
                chrom = "chr" + line[0] if "chr" not in line[0] else line[0]
                strand = line[5]
                beg = int(line[1])
                end = int(line[2])
                annotation_type = line[3].split(":")[0]
                annotation_name = line[3].split(":")[1]
                output_dict = {"chrom": chrom,
                               "strand": strand,
                               "beg": beg,
                               "end": end,
                               "annotation_type": annotation_type,
                               "annotation_name": annotation_name}
                output_list.append(output_dict)
    return output_list


def get_overlap(a, b):
    """Return overlap between two intervals (0-based, 1-based)"""
    return max(0, min(a[1], b[1]) - max(a[0], b[0]))


def filter_results(start, end, results, min_moverlap):
    """Filter the results according to the minimal overlap needed"""
    new_results = []
    for gene in results:
        if get_overlap((start, end), (gene.start, gene.end)) >= min_moverlap:
            new_results.append(gene)
    return new_results


if __name__ == '__main__':
    options, remaining_args = parser.parse_args()
    f = StringIO.StringIO()
    parser.print_help(f)
    helpstr = f.getvalue()
    mandatory_options = ["input",
                         "annotations"]
    check_mandatory_options(options, mandatory_options, helpstr)
    #
    # read the bed file with annotations
    #
    sys.stderr.write("Reading files\n")
    annotation_bed = read_annotation_bed_file(options.annotations)
    #
    # get intervaltree of the provided annotations
    #
    sys.stderr.write("Making index files\n")
    annotation_quicksect = index_bed(annotation_bed, options.un_stranded)
    #
    # scan bed file in order to get annotations
    #
    sys.stderr.write("Reading done!\n")
    with open(options.input) as bedfile, open(options.output, 'w') as outfile:
        reader = csv.reader(bedfile, delimiter="\t")
        for row in reader:
            try:
                #sys.stderr.write("\t".join(row) + "\n")
                chrom = row[0] if row[0] != 'chrM' else 'chrMT'
                strand = row[5]
                posstart = int(row[1])
                posend = int(row[2])
                min_moverlap = round((posend - posstart) * options.fraction)
                contig = chrom if options.un_stranded else chrom + ":" + strand
                #
                # scan IntervalTrees to find overlaps in genes, mirnas, exons, introns and repeats
                #
                annotation_scan = filter_results(posstart, posend, annotation_quicksect[contig].find(posstart, posend), min_moverlap)
                #
                # collect annotations into strings
                #
                annotation_types = set([i.value['annotation_type'] for i in annotation_scan])
                # annotation_type = ";".join(annotation_types)
                # annotation_name = ";".join(set([i.value['annotation_name'] for i in annotation_scan]))
                #
                # print line with annotations
                #
                # additional_columns = "\t".join(row[6:])
                for annotation in annotation_scan:
                    annotation_type = annotation.value['annotation_type']
                    annotation_name = annotation.value['annotation_name']
                    if options.filter_by:
                        if any(i in annotation_types for i in options.filter_by.split(",")):
                            outtext = "%s\t%i\t%i\t%s\t%s\t%s\t%s\t%s" % (chrom,
                                                                              posstart,
                                                                              posend,
                                                                              row[3],
                                                                              row[4],
                                                                              strand,
                                                                              annotation_type if annotation_type != "" else options.placeholder,
                                                                              annotation_name if annotation_name != "" else options.placeholder,)
                            outfile.write(outtext + "\n")
                    else:
                        outtext = "%s\t%i\t%i\t%s\t%s\t%s\t%s\t%s" % (chrom,
                                                                          posstart,
                                                                          posend,
                                                                          row[3],
                                                                          row[4],
                                                                          strand,
                                                                          annotation_type if annotation_type != "" else options.placeholder,
                                                                          annotation_name if annotation_name != "" else options.placeholder,)
                        outfile.write(outtext + "\n")
            except KeyError, e:
                # sys.stderr.write(str(e) + '\n')
                # continue
                pass
