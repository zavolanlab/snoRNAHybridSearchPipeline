#!/usr/bin/env python
"""

"""

__date_ = "2014-10-06"
__author__ = "Rafal Gumienny"
__email__ = "r.gumienny@unibas.ch"
__license__ = "GPL"

# imports
import sys
import time
import numpy as np
import pandas as pd
from collections import Counter
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
parser.add_argument("--mapping",
                    dest="mapping",
                    default="/import/bc2/home/zavolan/gumiennr/Pipelines/Pipelines/pipeline_snoRNASearch/data/Annotations/transcript_2_gene_mapping.txt.clean",
                    help="Mapping from ENSEMBL transcript to gene,\ndefaults to /import/bc2/home/zavolan/gumiennr/Pipelines/Pipelines/pipeline_snoRNASearch/data/Annotations/transcript_2_gene_mapping.txt.clean")

try:
    options = parser.parse_args()
except Exception, e:
    parser.print_help()

# redefine a functions for writing to stdout and stderr to save some writting
syserr = sys.stderr.write
sysout = sys.stdout.write


def main():
    """Main logic of the script"""
    mapping = pd.read_table(options.mapping)
    mapping.index = mapping['Ensembl Transcript ID']
    df = pd.read_table(options.input, header=None)
    gene_mapping = {}
    for key, value in mapping['Ensembl Gene ID'].to_dict().iteritems():
            gene_mapping[key] = value
    description_mapping = {}
    for key, value in mapping['WikiGene Description'].to_dict().iteritems():
            description_mapping[key] = value
    transcript_count_mapping = {}
    for key, value in mapping['Transcript count'].to_dict().iteritems():
            transcript_count_mapping[key] = value
    name_mapping = {}
    for key, value in mapping['WikiGene Name'].to_dict().iteritems():
            name_mapping[key] = value
    entrez_mapping = {}
    for key, value in mapping['EntrezGene ID'].to_dict().iteritems():
            entrez_mapping[key] = value


    df['EnsemblGeneID'] = [gene_mapping[i] if i in gene_mapping else np.nan for i in df[7]]
    df['EntrezGene ID'] = [entrez_mapping[i] if i in entrez_mapping else np.nan for i in df[7]]
    df['WikiName'] = [name_mapping[i] if i in name_mapping else np.nan for i in df[7]]
    df['Description'] = [description_mapping[i] if i in description_mapping else np.nan for i in df[7]]
    df['TranscriptCount'] = [transcript_count_mapping[i] if i in transcript_count_mapping else np.nan for i in df[7]]
    ndf = df.groupby([0, 1, 2, 3, 4, 5]).agg({7: max,
                                              6: lambda x: Counter(x).most_common()[0][0],
                                              'EnsemblGeneID': lambda x: x[0],
                                              'EntrezGene ID': lambda x: x[0],
                                              'Description': lambda x: x[0],
                                              'WikiName': lambda x: x[0],
                                              'TranscriptCount': max}).reset_index()
    ndf.sort(4, ascending=False).to_csv(options.output, sep='\t', header=None, index=None, na_rep='NaN')

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
