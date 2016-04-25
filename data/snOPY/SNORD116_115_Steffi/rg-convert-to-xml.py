#!/usr/bin/env python
"""
Convert snoRNA defined as a PLEXY input into snOPY XML file
"""

__date__ = "2015-10-27"
__author__ = "Rafal Gumienny"
__email__ = "r.gumienny@unibas.ch"
__license__ = "GPL"

# imports
import os
import sys
import time
from jinja2 import Template
from Bio import SeqIO
from contextlib import contextmanager
from argparse import ArgumentParser, RawTextHelpFormatter

parser = ArgumentParser(description=__doc__, formatter_class=RawTextHelpFormatter)
parser.add_argument("-v",
                    "--verbose",
                    dest="verbose",
                    action="store_true",
                    default=False,
                    help="Be loud!")
parser.add_argument("--template",
                    dest="template",
                    default="./template.xml",
                    help="Jinja2 template for xml, defaults to ./template.xml")
parser.add_argument("--snoRNAs",
                    dest="snoRNAs",
                    required=True,
                    help="Fasta file with snoRNAs")
parser.add_argument("--output-dir",
                    dest="output_dir",
                    default="output",
                    help="Output directory, defaults to ./")

try:
    options = parser.parse_args()
except Exception, e:
    parser.print_help()

# redefine a functions for writing to stdout and stderr to save some writting
syserr = sys.stderr.write
sysout = sys.stdout.write


def main():
    """Main logic of the script"""
    with open(options.template) as tmpl:
        template = Template(tmpl.read())
    with open(options.snoRNAs) as snors:
        for rec in SeqIO.parse(snors, 'fasta'):
            header = str(rec.id)
            sequence = str(rec.seq)
            infos = header.split("_")
            snorna_id = infos[0]
            if "sapiens" in infos[1]:
                organism = "Homo_sapiens"
            elif "musculus" in infos[1]:
                organism = "Mus_musculus"
            else:
                raise Exception("Unknown organism %s" % infos[1])
            chrom = infos[2][1:-1]
            chrom = chrom[3:] if chrom.startswith("chr") else chrom
            coords = infos[3]
            strand = infos[4]
            dpos = infos[-1]
            cpos = infos[-3]
            dppos = infos[-5]
            cppos = infos[-7]
            try:
                xml = template.render(snoopy_id="\"%s_%s\"" % (organism, snorna_id),
                                      data_quality="",
                                      open_status="",
                                      last_update=time.strftime("%d-%m-%Y"),
                                      blind_note="",
                                      organism=organism,
                                      gene_name=snorna_id,
                                      family_name="",
                                      alias="",
                                      family_base="",
                                      chrom_or_contig=chrom,
                                      genome_position=coords.replace(",", ".."),
                                      strand=strand,
                                      box_c="%i..%i,%i..%i" % (int(cppos), int(cppos) + 6, int(cpos), int(cpos) + 6),
                                      box_d="%i..%i,%i..%i" % (int(dppos), int(dppos) + 3, int(dpos), int(dpos) + 3),
                                      box_h="",
                                      box_aca="",
                                      host_id="",
                                      host_gene="",
                                      host_intron="",
                                      host_position="",
                                      organization="",
                                      mod_type="C/D",
                                      mod_rna="Unknonwn",
                                      mod_site="",
                                      duplex_region="",
                                      accession="",
                                      evidence="",
                                      note="",
                                      sequence=sequence)
            except ValueError:
                xml = template.render(snoopy_id="\"%s_%s\"" % (organism, snorna_id),
                                      data_quality="",
                                      open_status="",
                                      last_update=time.strftime("%d-%m-%Y"),
                                      blind_note="",
                                      organism=organism,
                                      gene_name=snorna_id,
                                      family_name="",
                                      alias="",
                                      family_base="",
                                      chrom_or_contig=chrom,
                                      genome_position=coords.replace(",", ".."),
                                      strand=strand,
                                      box_c="%i..%i" % (int(cpos), int(cpos) + 6),
                                      box_d="%i..%i" % (int(dpos), int(dpos) + 3),
                                      box_h="",
                                      box_aca="",
                                      host_id="",
                                      host_gene="",
                                      host_intron="",
                                      host_position="",
                                      organization="",
                                      mod_type="C/D",
                                      mod_rna="Unknonwn",
                                      mod_site="",
                                      duplex_region="",
                                      accession="",
                                      evidence="",
                                      note="",
                                      sequence=sequence)
            with open(os.path.join(options.output_dir, "%s%s.xml" % (organism, snorna_id)), 'w') as out:
                out.write(xml)

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
