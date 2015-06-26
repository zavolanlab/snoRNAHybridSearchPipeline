import os
import random
import sys
from errno import EEXIST

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
parser.add_argument("--name-suffix",
                    dest="name_suffix",
                    default="test_run",
                    help="Suffix to add to pipeline name in order to easily differentiate between different run, defaults to test_run")

try:
    options = parser.parse_args()
except Exception, e:
    parser.print_help()

settings = ConfigObj(options.config).dict()


def mkdir_p(path_to_dir):
    try:
        os.makedirs(path_to_dir)
    except OSError as e: # Python >2.5
        if e.errno == EEXIST and os.path.isdir(path_to_dir):
            sys.stderr.write("Directory %s already exists. Skipping.\n" % path_to_dir)
        else:
            raise e

working_directory = os.getcwd()
pipeline_directory = os.path.dirname(os.path.abspath(__file__))
output_directory = os.path.join(working_directory, "output")
index_directory = os.path.join(working_directory, "index")
plots_directory = os.path.join(working_directory, "Plots")
plexy_directory = os.path.join(working_directory, "Plexy")
mkdir_p(output_directory)
mkdir_p(index_directory)
mkdir_p(plots_directory)
mkdir_p(plexy_directory)
jobber_path = settings['general']['jobber_path']
sys.path.append(jobber_path)
from jobber import JobClient

jobber = JobClient.Jobber()

#Create a group for whole pipeline. The module "Python" will be inherited by all jobs that are in this group,
# so we don't need to define it for each job that calls a python script
pipeline_id = jobber.startGroup({'name': "snoRNAHybridSearch-%s" % options.name_suffix,
                                 'options': [['module', "Python"],
                                             ['module', "GCC"],
                                             ['module', "Bowtie2"],
                                             ['module', "OpenBLAS"],
                                             ['module', "BEDTools"],
                                             ['module', "SAMtools"],
                                             ],
                                 'executer': settings['general'].get('executer', 'drmaa')})

#First step is to split the file
split_command = "python %s --input %s --output-dir %s --batch-size %s -v" % (os.path.join(pipeline_directory, "scripts/rg-split-fasta.py"),
                                                          settings['general']['unmapped_reads'],
                                                          output_directory,
                                                          settings['general']['reads_per_file'],
                                                                             )
split_files_id = jobber.job(split_command, {'name': "SplitInput"})

#
# And the second is to generate the anchors - for that we need perform other tasks
#

# Generate file with snoRNAs
#
generate_fasta_command = "python %s --input %s --output %s --type %s -v --switch-box" % (
                                             os.path.join(pipeline_directory, "scripts/rg-generate-fasta.py"),
                                             settings['general']['snoRNAs'],
                                             os.path.join(working_directory, "snoRNAs.fa"),
                                             settings['general']['type'])
generate_fasta_id = jobber.job(generate_fasta_command, {'name': "GenerateFasta"})


#
# Generate plexy input files
generate_plexy_command = "python %s --input %s --dir %s --type %s -v --switch-box" % (
                                             os.path.join(pipeline_directory, "scripts/rg-generate-plexy.py"),
                                             settings['general']['snoRNAs'],
                                             plexy_directory,
                                             settings['general']['type'])
generate_plexy_id = jobber.job(generate_plexy_command, {'name': "GeneratePlexy"})

#
# Generate snoRNA bed

generate_snoRNA_bed_command = "python %s --input %s --output %s --type %s -v --switch-box" % (
                                             os.path.join(pipeline_directory, "scripts/rg-generate-snoRNA-bed.py"),
                                             settings['general']['snoRNAs'],
                                             os.path.join(working_directory, "snoRNAs.bed"),
                                             settings['general']['type'])
generate_snoRNA_bed_id = jobber.job(generate_snoRNA_bed_command, {'name': "GenerateSnoRNABed"})


#
# Annotate reads with snoRNAs


annotate_reads_with_snornas_tuple = (os.path.join(pipeline_directory, 'scripts/rg-annotate-bed.py'),
                                     settings['general']['bed_for_index'],
                                     os.path.join(working_directory, "mapped_reads_annotated_with_snornas.tab"),
                                     os.path.join(working_directory, "snoRNAs.bed"))

annotate_reads_with_snornas_command = "python %s --input %s --output %s --annotations %s --placeholder NaN --filter-by snoRNA" % annotate_reads_with_snornas_tuple


annotate_reads_with_snornas_id = jobber.job(annotate_reads_with_snornas_command, {'name': "AnnotateWithSnoRNAs",
                                        'dependencies': [generate_snoRNA_bed_id]})


#
# Calculate snoRNA expressions


calculate_snorna_expression_tuple = (os.path.join(pipeline_directory, 'scripts/rg-calculate-snoRNA-RPKM.py'),
                   os.path.join(working_directory, "mapped_reads_annotated_with_snornas.tab"),
                   os.path.join(working_directory, "snoRNAs.rpkm"),
                   settings['general']['bed_for_index'],
                   os.path.join(working_directory, "snoRNAs.fa"))

calculate_snorna_expression_command = "python %s --input %s --output %s --library %s --snoRNAs %s --quantile 0.0" % calculate_snorna_expression_tuple


calculate_snorna_expression_id = jobber.job(calculate_snorna_expression_command, {'name': "CalculateSnoRNAsExpression",
                                        'dependencies': [annotate_reads_with_snornas_id, generate_fasta_id]})


# Generate anchors
#
anchor_command = "python %s --fasta-to-anchor %s --anchor-length %s --output %s --expressed-snoRNAs %s -v" % (os.path.join(pipeline_directory, "scripts/rg-prepare-anchors.py"),
                                                          os.path.join(working_directory, "snoRNAs.fa"),
                                                          settings['general']['anchor_length'],
                                                          os.path.join(working_directory, "anchors.tab"),
                                                          os.path.join(working_directory, "snoRNAs.rpkm"),
                                                                             )
anchor_id = jobber.job(anchor_command, {'name': "PrepareAnchors",
                                        'dependencies': [generate_fasta_id, calculate_snorna_expression_id]})

#
# Build Bowtie index Group
#
build_index_id = jobber.startGroup({'name': "BuildIndex"})

#
# Cluster reads
#

cluster_reads_command = """python {script} \\
                            --input {input} \\
                            --output {output} \\
                            --overlap {overlap} \\
                            --expand-cluster {expand_cluster} \\
                            --expand-read {expand_read} \\
                            --cluster-size {cluster_size} \\
                            --bed \\
                            --rRNAs  {rRNAs} \\
                            --tRNAs  {tRNAs} \\
                            --snRNAs {snRNAs} \\
                            --filter-except snRNA,tRNA,rRNA,none,repeat
                        """.format(**{'script': os.path.join(pipeline_directory, 'scripts/rg-cluster-reads.py'),
                                      'input': settings['general']['bed_for_index'],
                                      'output': os.path.join(working_directory, "for_index.clustered"),
                                      'overlap': settings['tasks']['ClusterReads'].get('overlap', '1'),
                                      'expand_cluster': settings['tasks']['ClusterReads'].get('expand_cluster', '2'),
                                      'expand_read': settings['tasks']['ClusterReads'].get('expand_read', '15'),
                                      'cluster_size': settings['tasks']['ClusterReads'].get('cluster_size', '5'),
                                      'rRNAs':  settings['general']['rRNAs'],
                                      'tRNAs':  settings['general']['tRNAs'],
                                      'snRNAs': settings['general']['snRNAs'],
                                      })

cluster_reads_id = jobber.job(cluster_reads_command, {'name': "ClusterReads",
                                        'options': [('q', settings['tasks']['ClusterReads'].get('queue', 'short.q')),
                                                    ('l', "h_vmem=%s" % settings['tasks']['ClusterReads'].get('mem_req', '2G'))]
                                                      })

#
# Make fasta from clusters
#

make_fasta_tuple = (os.path.join(pipeline_directory, 'scripts/rg-extract-sequences.py'),
                   settings['general']['genome'],
                   os.path.join(working_directory, "for_index.clustered"),
                   os.path.join(working_directory, "for_index.fasta"))
make_fasta_command = "python %s --genome-dir %s --input %s --output %s --format fasta -v" % make_fasta_tuple

make_fasta_id = jobber.job(make_fasta_command, {'name': "MakeFasta",
                                                'dependencies': [cluster_reads_id],
                                                'options': [('q', settings['tasks']['MakeFastaFromClusters'].get('queue', 'short.q')),
                                                            ('l', "h_vmem=%s" % settings['tasks']['MakeFastaFromClusters'].get('mem_req', '2G'))]
                                                })


#
# Make Bowtie2 index
#

make_index_command =  "bowtie2-build %s %s/%s 2> /dev/null" % (os.path.join(working_directory, "for_index.fasta"),
                                                               index_directory,
                                                               "bowtie_index")
make_bowtie_index_id = jobber.job(make_index_command, {'name': "BuildBowtieIndex",
                                                'dependencies': [make_fasta_id],
                                                'options': [('q', settings['tasks']['BuildBowtieIndex'].get('queue', 'short.q')),
                                                            ('l', "h_vmem=%s" % settings['tasks']['BuildBowtieIndex'].get('mem_req', '2G'))]
                                                })

jobber.endGroup()
#######################################################################################

# We create a group where the jobs to analyse the splitted files will be put into
analyse_files_id = jobber.startGroup({'name': "Analysis",
                                      'dependencies': [split_files_id, anchor_id, build_index_id]})

#We call the script that will generate the jobs that will analyse the split files. We pass the id of the group
#and the folder where the script will find the splitted files.
analysis_tuple = (os.path.join(pipeline_directory, "run-analysis.py"),
                  output_directory,
                  analyse_files_id,
                  os.path.abspath(options.config),
                  working_directory)
analysis_command = "python %s --input-dir %s --group-id %s --config %s --working-dir %s -v" % analysis_tuple
jobber.job(analysis_command, {'name': "CreateAnalysisJobs"})


jobber.endGroup()

# We merge the files into our result file after analysis finishes
merge_results_command = "cat {output_dir}/*.plexybed > {cwd}/results_with_plexy.tab".format(output_dir=output_directory,
                                                                                                   cwd=working_directory)
merge_results_id = jobber.job(merge_results_command, {'name': "MergeResults",
                                 'dependencies': [analyse_files_id]})

# And we merge the raw read files into
merge_raw_results_command = "cat {output_dir}/*.truechrombed > {cwd}/raw_reads_results.tab".format(output_dir=output_directory,
                                                                                                   cwd=working_directory)
merge_raw_results_id = jobber.job(merge_raw_results_command, {'name': "MergeRawResults",
                                 'dependencies': [analyse_files_id]})


# # Cluster results
# #
# generate_fasta_command = "python %s --input %s --output %s" % (
#                                              os.path.join(pipeline_directory, "scripts/rg-generate-fasta.py"),
#                                                                os.path.join(working_directory, "results_with_plexy.tab"),
#                                                                os.path.join(working_directory, "results_with_plexy_clustered.tab")
#                                              )
# generate_fasta_id = jobber.job(generate_fasta_command, {'name': "GenerateFasta"})

# Add RPKM to results
#
add_rpkm_command = "python %s --input %s --output %s --rpkm %s --annotated-reads %s -v" % (
                                             os.path.join(pipeline_directory, "scripts/rg-add-rpkm-to-plexy.py"),
                                             os.path.join(working_directory, "results_with_plexy.tab"),
                                             os.path.join(working_directory, "results_with_plexy_and_rpkm.tab"),
                                             os.path.join(working_directory, "snoRNAs.rpkm"),
                                             os.path.join(working_directory, "mapped_reads_annotated_with_snornas.tab")
                                             )
add_rpkm_id = jobber.job(add_rpkm_command, {'name': "AddRPKM",
                                            'dependencies': [merge_results_id]})

# Cluster results
#
aggregate_by_site_command = "python %s --input %s --output %s" % (
                               os.path.join(pipeline_directory, "scripts/rg-aggregate-plexy-results.py"),
                                                 os.path.join(working_directory, "results_with_plexy_and_rpkm.tab"),
                                                 os.path.join(working_directory, "results_with_plexy_and_rpkm_aggregated.tab")
                               )
aggregate_by_site_id = jobber.job(aggregate_by_site_command, {'name': "AggregateBySite",
                                                              'dependencies': [add_rpkm_id]})
#
# Features for model
#

calculate_features_group_id = jobber.startGroup({'name': "FeaturesGroup",
                                                 'dependencies': [aggregate_by_site_id]})

template = settings['general']['template']
with open(template) as tmpl:
    template = Template(tmpl.read())

#
# Calculate accessibility of the sites

calculate_contrafold_settings = settings['tasks']['CalculateCONTRAFold']
calculate_contrafold_script = 'scripts/Contrafold.py'
calculate_contrafold_command = """python {script} \\
                                    --coords {coords} \\
                                    --out {output} \\
                                    --genome-dir {genome_dir} \\
                                    --contextLen_L 14 \\
                                    --contextLen_U 0 \\
                                    --context 30 \\
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
    copy_files = {os.path.join(working_directory, "results_with_plexy_and_rpkm_aggregated.tab"): 'input'}
    moveback = {'output': os.path.join(working_directory, "accessibility.tab")}

    calculate_contrafold_command_rendered = template.render(modules=calculate_contrafold_settings.get('modules', None),
                                       command=calculate_contrafold_command,
                                       copy=copy_files,
                                       moveback=moveback,
                                       copydir=copy_dir)
    calculate_contrafold_command = str(calculate_contrafold_command_rendered).format(**{'script': os.path.join(pipeline_directory, calculate_contrafold_script),
                                       'coords': 'input',
                                       'output': 'output',
                                       'genome_dir': settings['general']['genome'],
                                       })
else:
    calculate_contrafold_command = str(calculate_contrafold_command).format(**{'script': os.path.join(pipeline_directory, calculate_contrafold_script),
                                       'coords': os.path.join(working_directory, "results_with_plexy_and_rpkm_aggregated.tab"),
                                       'output': os.path.join(working_directory, "accessibility.tab"),
                                       'genome_dir': settings['general']['genome'],
                                       })
calculate_contrafold_id = jobber.job(calculate_contrafold_command,
                       {'name': 'CalculateCONTRAFold',
                        'uniqueId': True,
                        'options': [('q', calculate_contrafold_settings.get('queue', 'short.q')),
                                    ('l', "h_vmem=%s" % calculate_contrafold_settings.get('mem_req', '2G'))]
                        })

#
# Calulate Flanks composition

calculate_flanks_settings = settings['tasks']['CalculateFlanks']
calculate_flanks_script = 'scripts/Flanks_Composition.py'
calculate_flanks_command = """python {script} \\
                                    --coords {coords} \\
                                    --out {output} \\
                                    --genome-dir {genome_dir} \\
                                    --contextLen 30 \\
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
    copy_files = {os.path.join(working_directory, "results_with_plexy_and_rpkm_aggregated.tab"): 'input'}
    moveback = {'output': os.path.join(working_directory, "flanks.tab")}

    calculate_flanks_command_rendered = template.render(modules=calculate_flanks_settings.get('modules', None),
                                       command=calculate_flanks_command,
                                       copy=copy_files,
                                       moveback=moveback,
                                       copydir=copy_dir)
    calculate_flanks_command = str(calculate_flanks_command_rendered).format(**{'script': os.path.join(pipeline_directory, calculate_flanks_script),
                                       'coords': 'input',
                                       'output': 'output',
                                       'genome_dir': settings['general']['genome'],
                                       })
else:
    calculate_flanks_command = str(calculate_flanks_command).format(**{'script': os.path.join(pipeline_directory, calculate_flanks_script),
                                       'coords': os.path.join(working_directory, "results_with_plexy_and_rpkm_aggregated.tab"),
                                       'output': os.path.join(working_directory, "flanks.tab"),
                                       'genome_dir': settings['general']['genome'],
                                       })
calculate_flanks_id = jobber.job(calculate_flanks_command,
                       {'name': 'CalculateFlanks',
                        'uniqueId': True,
                        'options': [('q', calculate_flanks_settings.get('queue', 'short.q')),
                                    ('l', "h_vmem=%s" % calculate_flanks_settings.get('mem_req', '2G'))]
                        })

jobber.endGroup()


# Calculate Probability

calculate_probability_settings = settings['tasks']['CalculateProbability']
calculate_probability_command = "python %s --input %s --output %s --accessibility %s --flanks %s --model %s" % (
                                             os.path.join(pipeline_directory, "scripts/rg-calculate-probability.py"),
                                             os.path.join(working_directory, 'results_with_plexy_and_rpkm_aggregated.tab'),
                                             os.path.join(working_directory, 'results_with_probability.tab'),
                                             os.path.join(working_directory, "accessibility.tab"),
                                             os.path.join(working_directory, "flanks.tab"),
                                             settings['general']['model'])
calculate_probability_id = jobber.job(calculate_probability_command, {'name': "CalculateProbability",
                                    'options': [('q', calculate_probability_settings.get('queue', 'short.q')),
                                                ('l', "h_vmem=%s" % calculate_probability_settings.get('mem_req', '2G'))],
                                    'dependencies': [calculate_features_group_id,
                                                     aggregate_by_site_id]})

#
# Make statistics and plots for results

make_plots_command = """
                    python {script} \\
                        --results-probability-complex {complex_res} \\
                        --results-raw {raw_results} \\
                        --snoRNAs {snornas} \\
                        --type {type} \\
                        --genome-dir {genome_dir} \\
                        --dir {dir} \\
                        -v
                    """.format(**{'script': os.path.join(pipeline_directory, 'scripts/rg-make-stats-for-results.py'),
                                  'complex_res': os.path.join(working_directory, 'results_with_probability.tab'),
                                  'raw_results': os.path.join(working_directory, 'raw_reads_results.tab'),
                                  'snornas': settings['general']['snoRNAs'],
                                  'type': settings['general']['type'],
                                  'genome_dir': settings['general']['genome'],
                                  'dir': os.path.join(working_directory, "Plots")
                                  })
make_plots_id = jobber.job(make_plots_command, {'name': "MakePlots",
                                    'dependencies': [merge_raw_results_id,
                                                     calculate_probability_id]})


#
# Convert results to BED
convert_to_bed_tuple =  (os.path.join(pipeline_directory, 'scripts/rg-convert-to-bed.py'),
                         os.path.join(working_directory, "results_with_probability.tab"),
                         os.path.join(working_directory, "results_with_probability.bed"))

convert_to_bed_command = "python %s --input %s --output %s" % convert_to_bed_tuple


convert_to_bed_id = jobber.job(convert_to_bed_command, {'name': "ConvertToBed",
                                        'dependencies': [calculate_probability_id]})

#
# Annotate results
annotate_results_settings = settings['tasks']['AnnotateResults']
annotate_results_tuple =  (os.path.join(pipeline_directory, 'scripts/rg-annotate-positions.py'),
                           os.path.join(working_directory, "results_with_probability.bed"),
                           os.path.join(working_directory, "results_with_probability_annotated.tab"),
                           settings['general']['annotations_genes'],
                           settings['general']['annotations_regions'],
                           )

annotate_results_command = "python %s --input %s --output %s --genes %s --regions %s" % annotate_results_tuple


annotate_results_id = jobber.job(annotate_results_command, {'name': "AnnotateResults",
                                    'options': [('q', annotate_results_settings.get('queue', 'short.q')),
                                                ('l', "h_vmem=%s" % annotate_results_settings.get('mem_req', '2G'))],
                                        'dependencies': [convert_to_bed_id]})


############################### END PIPELINE GROUP ###################################
jobber.endGroup()

# Before launching we print the command to stop the pipeline
print "In order to stop the pipeline run a command:"
print "python %s/jobber_server.py -command delete -jobId %i" % (jobber_path, pipeline_id)

#You need to always launch, otherwise jobs wont get executed.
jobber.launch(pipeline_id)
