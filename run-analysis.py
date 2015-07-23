#This script generates the jobs of the pipeline that will do the analysis on the splitted files.
import glob
import os
import sys
from configobj import ConfigObj
from jinja2 import Template
from argparse import ArgumentParser, RawTextHelpFormatter

parser = ArgumentParser(description=__doc__, formatter_class=RawTextHelpFormatter)
parser.add_argument("-v",
                    "--verbose",
                    dest="verbose",
                    action="store_true",
                    default=False,
                    help="Be loud!")
parser.add_argument("--config",
                    dest="config",
                    required=True,
                    help="Config file")
parser.add_argument("--group-id",
                    dest="group_id",
                    required=True,
                    help="Group Id")
parser.add_argument("--input-dir",
                    dest="input_dir",
                    required=True,
                    help="Input and output directory")
parser.add_argument("--working-dir",
                    dest="working_dir",
                    required=True,
                    help="Working directory of the pipeline. Required because this file is launched from ~")

try:
    options = parser.parse_args()
except Exception, e:
    parser.print_help()

settings = ConfigObj(options.config).dict()

working_directory = options.working_dir
index_directory = os.path.join(working_directory, 'index')
output_directory = options.input_dir
pip_dir = os.path.dirname(os.path.abspath(__file__))
jobberDir = settings['general']['jobber_path']
sys.path.append(jobberDir)
from jobber import JobClient

jobber = JobClient.Jobber()

#We can call "extendGroup" if we want to create jobs into an already existing group. Don't forget to call "endGroup" and "launch"
#after you're done
jobber.extendGroup(options.group_id)

#We assume that each file to be analyzed ends with .seqs. Its important to always distinguish input files from any intermediary
#files in case we need to restart the jobs. We should make the jobs unique to prevent duplication of jobs in case this script
#is run multiple times

template = settings['general']['template']
with open(template) as tmpl:
    template = Template(tmpl.read())

files_to_run = {}
search_anchors_group = jobber.startGroup({'name': 'SearchAnchors'})
for f in glob.glob(options.input_dir + "/*.inputfasta"):
    input_name = os.path.splitext(f)[0]

    #
    # Search Anchors in unmapped reads
    #
    search_anchors_settings = settings['tasks']['SearchAnchors']
    search_anchors_script = 'scripts/rg-search-anchor-and-make-alignments.py'
    search_anchors_command = """python {script} \\
                                    --anchors {anchors} \\
                                    --anchor-sequences {anchor_sequences} \\
                                    --reads {reads} \\
                                    --match {match} \\
                                    --mismatch {mismatch} \\
                                    --gap-open {gap_open} \\
                                    --gap-extend {gap_extend} \\
                                    --output {output} \\
                                    -v
                  """
    #
    # If there is template use it for command
    #
    if settings['general'].get('executer', 'drmaa') == 'drmaa':
        #
        # Copy files by default to the tmp directory
        #
        copy_dir = "$TMPDIR"
        copy_files = {f: 'reads.fa',
                      os.path.join(working_directory, "anchors.tab"): "anchors.tab",
                      os.path.join(working_directory, "snoRNAs.fa"): 'anchor_sequences.fa'}
        moveback = {'output': input_name + ".anchorsearch"}

        search_anchors_command_rendered = template.render(modules=search_anchors_settings.get('modules', None),
                                           command=search_anchors_command,
                                           copy=copy_files,
                                           moveback=moveback,
                                           copydir=copy_dir)
        search_anchors_command = str(search_anchors_command_rendered).format(**{'script': os.path.join(pip_dir, search_anchors_script),
                                         'reads': 'reads.fa',
                                         'anchors': 'anchors.tab',
                                         'anchor_sequences': 'anchor_sequences.fa',
                                         'output': 'output',
                                         'match': search_anchors_settings.get('match', '2'),
                                         'mismatch': search_anchors_settings.get('mismatch', '-5'),
                                         'gap_open': search_anchors_settings.get('gap_open', '-6'),
                                         'gap_extend': search_anchors_settings.get('gap_extended', '-4'),
                                       })
    else:
        search_anchors_command = str(search_anchors_command).format(**{'script': os.path.join(pip_dir, search_anchors_script),
                                       'reads': f,
                                       'anchors': os.path.join(working_directory, "anchors.tab"),
                                       'anchor_sequences': os.path.join(working_directory, "snoRNAs.fa"),
                                       'output': input_name + ".anchorsearch",
                                       'match': search_anchors_settings.get('match', '2'),
                                       'mismatch': search_anchors_settings.get('mismatch', '-5'),
                                       'gap_open': search_anchors_settings.get('gap_open', '-6'),
                                       'gap_extend': search_anchors_settings.get('gap_extended', '-4'),
                                       })
    search_anchors_id = jobber.job(search_anchors_command,
                               {'name': 'AnchorSearch',
                                'uniqueId': True,
                                'options': [('q', search_anchors_settings.get('queue', 'short.q')),
                                            ('l', "membycore=%s" % search_anchors_settings.get('mem_req', '2G'))]
                                })
    files_to_run[input_name] = search_anchors_id

jobber.endGroup()

#
# Merge for statistics
#

make_statistics_group = jobber.startGroup({'name': "MakeStatistics",
                                           'dependencies': [search_anchors_group]})
merge_for_statistics_command = "cat {output_dir}/*.anchorsearch > {cwd}/search_anchors_for_statistics".format(
                                    output_dir=output_directory,
                                    cwd=working_directory)
merge_for_statistics_id = jobber.job(merge_for_statistics_command, {'name': "MergeForStats"})

#
# Make statistics
#


make_statistics_command = "python %s --input %s --output %s --length 15 --fpr 1.0 --dir %s" % (
            os.path.join(pip_dir, 'scripts/rg-make-stats-for-search.py'),
            os.path.join(working_directory, "search_anchors_for_statistics"),
            os.path.join(working_directory, "search_anchors.stats"),
            os.path.join(working_directory, "Plots")
              )

make_statistics_id = jobber.job(make_statistics_command, {'name': "MakeStats",
                                                          'options': [('l', "membycore=10G")],
                                                          'dependencies': [merge_for_statistics_id]})

jobber.endGroup()

#
# Convert to fasta
#
convert_to_fasta_group = jobber.startGroup({'name': 'ConvertToFasta',
                                            'dependencies': [make_statistics_group]})
convert_dependancies = {}
for input_name, search_anchor_id in files_to_run.iteritems():
    convert_to_fasta_settings = settings['tasks']['ConvertToFasta']
    convert_to_fasta_script = 'scripts/rg-convert-tab-to-fasta.py'
    convert_to_fasta_command = """python {script} \\
                    --input {input} \\
                    --output {output} \\
                    --length {length} \\
                    --assign-score-threshold \\
                    --stats {stats} \\
                    -v
              """

    if settings['general'].get('executer', 'drmaa') == 'drmaa':
        #
        # Copy files by default to the tmp directory
        #
        copy_dir = "$TMPDIR"
        copy_files = {input_name + ".anchorsearch": 'input.anchorsearch',
                      os.path.join(working_directory, "search_anchors.stats"): 'input.stats',
                      }
        moveback = {'output': input_name + ".anchorfasta"}

        convert_to_fasta_command_rendered = template.render(modules=convert_to_fasta_settings.get('modules', None),
                                                           command=convert_to_fasta_command,
                                                           copy=copy_files,
                                                           moveback=moveback,
                                                           copydir=copy_dir)
        convert_to_fasta_command = str(convert_to_fasta_command_rendered).format(**{'script': os.path.join(pip_dir, convert_to_fasta_script),
                            'output': 'output',
                            'input': "input.anchorsearch",
                            'stats': "input.stats",
                            'length': convert_to_fasta_settings.get('length', "15"),
                            })
    else:
        convert_to_fasta_command = str(convert_to_fasta_command).format(**{'script': os.path.join(pip_dir, convert_to_fasta_script),
                            'output': input_name + ".anchorfasta",
                            'input': input_name + ".anchorsearch",
                            'stats': os.path.join(working_directory, "search_anchors.stats"),
                            'length': convert_to_fasta_settings.get('length', "15"),
                            })

    convert_to_fasta_id = jobber.job(convert_to_fasta_command, {
                                      'name': 'ConvertToFasta',
                                      'dependencies': [search_anchor_id],
                                       'options': [('q', convert_to_fasta_settings.get('queue', 'short.q')),
                                                   ('l', "membycore=%s" % convert_to_fasta_settings.get('mem_req', '2G'))],
                                      'uniqueId': True})
    convert_dependancies[input_name] = convert_to_fasta_id

jobber.endGroup()


#
# Map reads
#
map_reads_group = jobber.startGroup({'name': 'MapReads'})
map_dependancies = {}
for input_name, convert_id in convert_dependancies.iteritems():
    map_reads_settings = settings['tasks']['MapReads']
    map_reads_script = 'bowtie2'
    map_reads_command = """{script} -x {index} -f -D100 -L 13 -i C,1 --local -k 10 -U {input} -S {output}"""

    if settings['general'].get('executer', 'drmaa') == 'drmaa':
        #
        # Copy files by default to the tmp directory
        #
        copy_dir = "$TMPDIR"
        copy_files = {input_name + ".anchorfasta": 'input.anchorfasta',
                      index_directory: "index"}
        moveback = {'output.sam': input_name + ".sam"}

        map_reads_command_rendered = template.render(modules=map_reads_settings.get('modules', None),
                                                           command=map_reads_command,
                                                           copy=copy_files,
                                                           moveback=moveback,
                                                           copydir=copy_dir)
        map_reads_command = str(map_reads_command_rendered).format(**{'script': map_reads_script,
                            'output': 'output.sam',
                            'input': "input.anchorfasta",
                            'index': "./index/bowtie_index",
                            })
    else:
        map_reads_command = str(map_reads_command).format(**{'script': map_reads_script,
                            'output': input_name + ".sam",
                            'input': input_name + ".anchorfasta",
                            'index': os.path.join(index_directory, "bowtie_index"),
                            })

    map_reads_id = jobber.job(map_reads_command, {
                                      'name': 'MapReads',
                                      'dependencies': [convert_id],
                                       'options': [('q', map_reads_settings.get('queue', 'short.q')),
                                                   ('l', "membycore=%s" % map_reads_settings.get('mem_req', '2G'))],
                                      'uniqueId': True})
    map_dependancies[input_name] = map_reads_id

jobber.endGroup()


# #
# # Filter mapped reads
# #
# filter_mapped_group = jobber.startGroup({'name': 'FilterMapped'})
# filter_map_dependancies = {}
# for input_name, map_id in map_dependancies.iteritems():
#     filter_mapped_settings = settings['tasks']['FilterMapped']
#     filter_mapped_script = 'grep'
#     filter_mapped_command = """{script} -v XS:i: {input} > {output}"""

#     if settings['general'].get('executer', 'drmaa') == 'drmaa':
#         #
#         # Copy files by default to the tmp directory
#         #
#         copy_dir = "$TMPDIR"
#         copy_files = {input_name + ".sam": 'input.sam'}
#         moveback = {'output': input_name + ".samfiltered"}

#         filter_mapped_command_rendered = template.render(modules=filter_mapped_settings.get('modules', None),
#                                                            command=filter_mapped_command,
#                                                            copy=copy_files,
#                                                            moveback=moveback,
#                                                            copydir=copy_dir)
#         filter_mapped_command = str(filter_mapped_command_rendered).format(**{'script': filter_mapped_script,
#                             'output': 'output',
#                             'input': "input.sam",
#                             })
#     else:
#         filter_mapped_command = str(filter_mapped_command).format(**{'script': filter_mapped_script,
#                             'output': input_name + ".samfiltered",
#                             'input': input_name + ".sam",
#                             })

#     filter_mapped_id = jobber.job(filter_mapped_command, {
#                                       'name': 'FilterMapped',
#                                       'dependencies': [map_id],
#                                        'options': [('q', filter_mapped_settings.get('queue', 'short.q')),
#                                                    ('l', "membycore=%s" % filter_mapped_settings.get('mem_req', '2G'))],
#                                       'uniqueId': True})
#     filter_map_dependancies[input_name] = filter_mapped_id

# jobber.endGroup()


# #
# # ConvertToBed
# #
# convert_mapped_to_bed_group = jobber.startGroup({'name': 'ConvertToBed'})
# convert_to_bed_dependancies = {}
# # for input_name, filter_map_id in filter_map_dependancies.iteritems():
# for input_name, map_id in map_dependancies.iteritems():
#     convert_mapped_to_bed_settings = settings['tasks']['ConvertToBed']
#     convert_mapped_to_bed_script = 'samtools view'
#     convert_mapped_to_bed_command = """{script} -S {input} -b -u | bamToBed -tag AS | grep -P \"\\t\\+\" > {output}"""

#     if settings['general'].get('executer', 'drmaa') == 'drmaa':
#         #
#         # Copy files by default to the tmp directory
#         #
#         copy_dir = "$TMPDIR"
#         copy_files = {input_name + ".samfiltered": 'input.samfiltered'}
#         moveback = {'output': input_name + ".mappedbed"}

#         convert_mapped_to_bed_command_rendered = template.render(modules=convert_mapped_to_bed_settings.get('modules', None),
#                                                            command=convert_mapped_to_bed_command,
#                                                            copy=copy_files,
#                                                            moveback=moveback,
#                                                            copydir=copy_dir)
#         convert_mapped_to_bed_command = str(convert_mapped_to_bed_command_rendered).format(**{'script': convert_mapped_to_bed_script,
#                             'output': 'output',
#                             'input': "input.samfiltered",
#                             })
#     else:
#         convert_mapped_to_bed_command = str(convert_mapped_to_bed_command).format(**{'script': convert_mapped_to_bed_script,
#                             'output': input_name + ".mappedbed",
#                             'input': input_name + ".samfiltered",
#                             })

#     convert_mapped_to_bed_id = jobber.job(convert_mapped_to_bed_command, {
#                                       'name': 'ConvertToBed',
#                                       'dependencies': [filter_map_id],
#                                        'options': [('q', convert_mapped_to_bed_settings.get('queue', 'short.q')),
#                                                    ('l', "membycore=%s" % convert_mapped_to_bed_settings.get('mem_req', '2G'))],
#                                       'uniqueId': True})
#     convert_to_bed_dependancies[input_name] = convert_mapped_to_bed_id

# jobber.endGroup()

#
# ConvertToBed
#
convert_mapped_to_bed_group = jobber.startGroup({'name': 'ConvertToBed'})
convert_to_bed_dependancies = {}
for input_name, map_id in map_dependancies.iteritems():
    convert_mapped_to_bed_settings = settings['tasks']['ConvertToBed']
    convert_mapped_to_bed_script = 'samtools view'
    convert_mapped_to_bed_command = """{script} -S {input} -b -u | bamToBed -tag AS | grep -P \"\\t\\+\" > {output}"""

    if settings['general'].get('executer', 'drmaa') == 'drmaa':
        #
        # Copy files by default to the tmp directory
        #
        copy_dir = "$TMPDIR"
        copy_files = {input_name + ".sam": 'input.sam'}
        moveback = {'output': input_name + ".mappedbed"}

        convert_mapped_to_bed_command_rendered = template.render(modules=convert_mapped_to_bed_settings.get('modules', None),
                                                           command=convert_mapped_to_bed_command,
                                                           copy=copy_files,
                                                           moveback=moveback,
                                                           copydir=copy_dir)
        convert_mapped_to_bed_command = str(convert_mapped_to_bed_command_rendered).format(**{'script': convert_mapped_to_bed_script,
                            'output': 'output',
                            'input': "input.sam",
                            })
    else:
        convert_mapped_to_bed_command = str(convert_mapped_to_bed_command).format(**{'script': convert_mapped_to_bed_script,
                            'output': input_name + ".mappedbed",
                            'input': input_name + ".sam",
                            })

    convert_mapped_to_bed_id = jobber.job(convert_mapped_to_bed_command, {
                                      'name': 'ConvertToBed',
                                      'dependencies': [map_id],
                                       'options': [('q', convert_mapped_to_bed_settings.get('queue', 'short.q')),
                                                   ('l', "membycore=%s" % convert_mapped_to_bed_settings.get('mem_req', '2G'))],
                                      'uniqueId': True})
    convert_to_bed_dependancies[input_name] = convert_mapped_to_bed_id

jobber.endGroup()

#
# Filter mapped bed
#
filter_bed_group = jobber.startGroup({'name': 'FilterBed'})
filter_bed_dependancies = {}
for input_name, convert_mapped_to_bed_id in convert_to_bed_dependancies.iteritems():
    filter_bed_settings = settings['tasks']['FilterBed']
    filter_bed_script = 'scripts/rg-filter-bed.py'
    filter_bed_command = """python {script} \\
                                --input {input} \\
                                --output {output} \\
                                -v
                      """

    if settings['general'].get('executer', 'drmaa') == 'drmaa':
        #
        # Copy files by default to the tmp directory
        #
        copy_dir = "$TMPDIR"
        copy_files = {input_name + ".mappedbed": 'input.mappedbed'}
        moveback = {'output': input_name + ".filteredbed"}

        filter_bed_command_rendered = template.render(modules=filter_bed_settings.get('modules', None),
                                                           command=filter_bed_command,
                                                           copy=copy_files,
                                                           moveback=moveback,
                                                           copydir=copy_dir)
        filter_bed_command = str(filter_bed_command_rendered).format(**{'script': os.path.join(pip_dir, filter_bed_script),
                            'output': 'output',
                            'input': "input.mappedbed",
                            })
    else:
        filter_bed_command = str(filter_bed_command).format(**{'script': os.path.join(pip_dir, filter_bed_script),
                            'output': input_name + ".filteredbed",
                            'input': input_name + ".mappedbed",
                            })

    filter_bed_id = jobber.job(filter_bed_command, {
                                      'name': 'FilterBed',
                                      'dependencies': [convert_mapped_to_bed_id],
                                       'options': [('q', filter_bed_settings.get('queue', 'short.q')),
                                                   ('l', "membycore=%s" % filter_bed_settings.get('mem_req', '2G'))],
                                      'uniqueId': True})
    filter_bed_dependancies[input_name] = filter_bed_id

jobber.endGroup()

#
# Calculate true positions on the chromosome
#
get_true_chrom_dependancies = {}
get_true_chrom_group = jobber.startGroup({'name': 'GetTrueChromosome'})
for input_name, filter_bed_id in filter_bed_dependancies.iteritems():
    get_true_chrom_settings = settings['tasks']['GetTrueChromosome']
    get_true_chrom_script = 'scripts/rg-get-true-chromosome-positions.py'
    get_true_chrom_command = """python {script} \\
                                    --input {input} \\
                                    --output {output} \\
                                    -v
                  """
    #
    # If there is template use it for command
    #
    if settings['general'].get('executer', 'drmaa') == 'drmaa':
        #
        # Copy files by default to the tmp directory
        #
        copy_dir = "$TMPDIR"
        copy_files = {input_name + ".filteredbed": 'input.filteredbed'}
        moveback = {'output': input_name + ".truechrombed"}

        get_true_chrom_command_rendered = template.render(modules=get_true_chrom_settings.get('modules', None),
                                           command=get_true_chrom_command,
                                           copy=copy_files,
                                           moveback=moveback,
                                           copydir=copy_dir)
        get_true_chrom_command = str(get_true_chrom_command_rendered).format(**{'script': os.path.join(pip_dir, get_true_chrom_script),
                                         'input': 'input.filteredbed',
                                         'output': 'output',
                                       })
    else:
        get_true_chrom_command = str(get_true_chrom_command).format(**{'script': os.path.join(pip_dir, get_true_chrom_script),
                                       'input': input_name + ".filteredbed",
                                       'output': input_name + ".truechrombed",
                                       })
    get_true_chrom_id = jobber.job(get_true_chrom_command,
                               {'name': 'GetTrueChromosome',
                                'uniqueId': True,
                                 'dependencies': [filter_bed_id],
                                'options': [('q', get_true_chrom_settings.get('queue', 'short.q')),
                                            ('l', "membycore=%s" % get_true_chrom_settings.get('mem_req', '2G'))]
                                })
    get_true_chrom_dependancies[input_name] = get_true_chrom_id

jobber.endGroup()


#
# Append sequence
#
append_sequence_dependancies = {}
append_sequence_group = jobber.startGroup({'name': 'AppendSequence'})
for input_name, get_true_chrom_id in get_true_chrom_dependancies.iteritems():
    append_sequence_settings = settings['tasks']['AppendSequence']
    append_sequence_script = 'scripts/rg-extract-sequences.py'
    append_sequence_command = """python {script} \\
                                    --input {input} \\
                                    --output {output} \\
                                    --genome-dir {genome} \\
                                    --sequence-length {sequence_length} \\
                                    -v
                  """
    #
    # If there is template use it for command
    #
    if settings['general'].get('executer', 'drmaa') == 'drmaa':
        #
        # Copy files by default to the tmp directory
        #
        copy_dir = "$TMPDIR"
        copy_files = {input_name + ".truechrombed": 'input.truechrombed'}
        moveback = {'output': input_name + ".sequencebed"}

        append_sequence_command_rendered = template.render(modules=append_sequence_settings.get('modules', None),
                                           command=append_sequence_command,
                                           copy=copy_files,
                                           moveback=moveback,
                                           copydir=copy_dir)
        append_sequence_command = str(append_sequence_command_rendered).format(**{'script': os.path.join(pip_dir, append_sequence_script),
                                         'input': 'input.truechrombed',
                                         'output': 'output',
                                         'genome': settings['general']['genome'],
                                         'sequence_length': append_sequence_settings.get('sequence_length', '50')
                                       })
    else:
        append_sequence_command = str(append_sequence_command).format(**{'script': os.path.join(pip_dir, append_sequence_script),
                                       'input': input_name + ".truechrombed",
                                       'output': input_name + ".sequencebed",
                                       'genome': settings['general']['genome'],
                                       'sequence_length': append_sequence_settings.get('sequence_length', '50')
                                       })
    append_sequence_id = jobber.job(append_sequence_command,
                               {'name': 'AppendSequence',
                                'uniqueId': True,
                                 'dependencies': [get_true_chrom_id],
                                'options': [('q', append_sequence_settings.get('queue', 'short.q')),
                                            ('l', "membycore=%s" % append_sequence_settings.get('mem_req', '2G'))]
                                })
    append_sequence_dependancies[input_name] = append_sequence_id

jobber.endGroup()


if settings["general"]['type'] == "CD":
    #
    # Append PLEXY
    #

    append_score_dependancies = {}
    append_score_group = jobber.startGroup({'name': 'AppendPLEXY'})
    for input_name, append_sequence_id in append_sequence_dependancies.iteritems():
        append_score_settings = settings['tasks']['AppendPLEXY']
        append_score_script = 'scripts/rg-check-hybrids-with-plexy.py'
        append_score_command = """python {script} \\
                                        --input {input} \\
                                        --output {output} \\
                                        --snoRNA-paths {snorna_path} \\
                                        --plexy-bin {plexy_bin} \\
                                        --plexy-tmp {plexy_tmp} \\
                                        -v
                      """
        #
        # If there is template use it for command
        #
        if settings['general'].get('executer', 'drmaa') == 'drmaa':
            #
            # Copy files by default to the tmp directory
            #
            copy_dir = "$TMPDIR"
            copy_files = {input_name + ".sequencebed": 'input.sequencebed'}
            moveback = {'output': input_name + ".scorebed"}

            append_score_command_rendered = template.render(modules=append_score_settings.get('modules', None),
                                               command=append_score_command,
                                               copy=copy_files,
                                               moveback=moveback,
                                               copydir=copy_dir)
            append_score_command = str(append_score_command_rendered).format(**{'script': os.path.join(pip_dir, append_score_script),
                                            'input': 'input.sequencebed',
                                            'output': 'output',
                                            'snorna_path': os.path.join(working_directory, "Input/"),
                                            'plexy_bin': settings['general']['PLEXY_bin'],
                                            'plexy_tmp': './',
                                           })
        else:
            append_score_command = str(append_score_command).format(**{'script': os.path.join(pip_dir, append_score_script),
                                           'input': input_name + ".sequencebed",
                                           'output': input_name + ".scorebed",
                                           'snorna_path': os.path.join(working_directory, "Input/"),
                                           'plexy_bin': settings['general']['PLEXY_bin'],
                                           'plexy_tmp': './tmp',
                                           })
        append_score_id = jobber.job(append_score_command,
                                   {'name': 'AppendPLEXY',
                                    'uniqueId': True,
                                    'dependencies': [append_sequence_id],
                                    'options': [('q', append_score_settings.get('queue', 'short.q')),
                                                ('l', "membycore=%s" % append_score_settings.get('mem_req', '2G'))]
                                    })
        append_score_dependancies[input_name] = append_score_id

    jobber.endGroup()
elif settings["general"]['type'] == "HACA":
    #
    # Append RNAsnoop
    #

    append_score_dependancies = {}
    append_score_group = jobber.startGroup({'name': 'AppendRNAsnoop'})
    for input_name, append_sequence_id in append_sequence_dependancies.iteritems():
        append_score_settings = settings['tasks']['AppendRNAsnoop']
        append_score_script = 'scripts/rg-check-hybrids-with-rnasnoop.py'
        append_score_command = """python {script} \\
                                        --input {input} \\
                                        --output {output} \\
                                        --snoRNA-paths {snorna_path} \\
                                        --rnasnoop {rnasnoop_bin} \\
                                        -v
                      """
        #
        # If there is template use it for command
        #
        if settings['general'].get('executer', 'drmaa') == 'drmaa':
            #
            # Copy files by default to the tmp directory
            #
            copy_dir = "$TMPDIR"
            copy_files = {input_name + ".sequencebed": 'input.sequencebed'}
            moveback = {'output': input_name + ".scorebed"}

            append_score_command_rendered = template.render(modules=append_score_settings.get('modules', None),
                                               command=append_score_command,
                                               copy=copy_files,
                                               moveback=moveback,
                                               copydir=copy_dir)
            append_score_command = str(append_score_command_rendered).format(**{'script': os.path.join(pip_dir, append_score_script),
                                            'input': 'input.sequencebed',
                                            'output': 'output',
                                            'snorna_path': os.path.join(working_directory, "Input/"),
                                            'rnasnoop_bin': settings['general']['RNAsnoop_bin'],
                                           })
        else:
            append_score_command = str(append_score_command).format(**{'script': os.path.join(pip_dir, append_score_script),
                                           'input': input_name + ".sequencebed",
                                           'output': input_name + ".scorebed",
                                           'snorna_path': os.path.join(working_directory, "Input/"),
                                           'rnasnoop_bin': settings['general']['RNAsnoop_bin'],
                                           })
        append_score_id = jobber.job(append_score_command,
                                   {'name': 'AppendRNAsnoop',
                                    'uniqueId': True,
                                    'dependencies': [append_sequence_id],
                                    'options': [('q', append_score_settings.get('queue', 'short.q')),
                                                ('l', "membycore=%s" % append_score_settings.get('mem_req', '2G'))]
                                    })
        append_score_dependancies[input_name] = append_score_id

    jobber.endGroup()
else:
    raise Exception("Unknown snoRNA type: %s" % settings["general"]['type'])


    # #
    # # Calculate seed matches
    # #
    # seed_count_settings = settings['tasks']['CalculateSeedMatches']
    # seed_count_script = 'scripts/rg-count-miRNA-seeds-and-filter-duplicates.py'
    # seed_count_command = """python {script} \\
    #                             --motifs {input} \\
    #                             --seqs {seqs} \\
    #                             --output {output} \\
    #                             --how {how} \\
    #                             --context {context} \\
    #                             --split-by "{split_by}" \\
    #                             --index-after-split {index_after_split} \\
    #                             -v
    #               """
    # #
    # # If there is template use it for command
    # #
    # if settings['general'].get('executer', 'drmaa') == 'drmaa':
    #     #
    #     # Copy files by default to the tmp directory
    #     #
    #     copy_dir = "$TMPDIR"
    #     copy_files = {f: 'mirnas.fa',
    #                   settings['general']['seqs']: 'seqs.fa'}
    #     moveback = {'output': input_name + ".seedcount"}

    #     seed_command_rendered = template.render(modules=seed_count_settings.get('modules', None),
    #                                        command=seed_count_command,
    #                                        copy=copy_files,
    #                                        moveback=moveback,
    #                                        copydir=copy_dir)
    #     seed_count_command = str(seed_command_rendered).format(**{'script': os.path.join(pip_dir, seed_count_script),
    #                                    'input': 'mirnas.fa',
    #                                    'seqs': 'seqs.fa',
    #                                    'output': 'output',
    #                                    'how': seed_count_settings.get('how', 'TargetScan'),
    #                                     'split_by': seed_count_settings.get('split_by', "NONE"),
    #                                    'index_after_split': seed_count_settings.get('index_after_split', 0),
    #                                    'context': seed_count_settings.get('context', 50)})
    # else:
    #     seed_count_command = str(seed_count_command).format(**{'script': os.path.join(pip_dir, seed_count_script),
    #                                    'input': f,
    #                                    'seqs': settings['general']['seqs'],
    #                                    'output': input_name + ".seedcount",
    #                                    'how': seed_count_settings.get('how', 'TargetScan'),
    #                                     'split_by': seed_count_settings.get('split_by', "NONE"),
    #                                    'index_after_split': seed_count_settings.get('index_after_split', 0),
    #                                    'context': seed_count_settings.get('context', 50)})
    # seed_count_id = jobber.job(seed_count_command,
    #                            {'name': 'SeedCount',
    #                             'uniqueId': True,
    #                             'options': [('q', seed_count_settings.get('queue', 'short.q')),
    #                                         ('l', "membycore=%s" % seed_count_settings.get('mem_req', '2G'))]
    #                             })



jobber.endGroup()
jobber.launch(options.group_id)
