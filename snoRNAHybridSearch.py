import os
import random
import sys
import shutil
import traceback
from errno import EEXIST

from configobj import ConfigObj
from jinja2 import Template
from argparse import ArgumentParser, RawTextHelpFormatter

parser = ArgumentParser(description=__doc__, formatter_class=RawTextHelpFormatter)
subparsers = parser.add_subparsers(help="Commands", dest='command')
run_parser = subparsers.add_parser("run", help="Run a pipeline")
run_parser.add_argument("-v",
                    "--verbose",
                    dest="verbose",
                    action="store_true",
                    default=False,
                    help="Be loud!")
run_parser.add_argument("--config",
                    dest="config",
                    required=True,
                    help="Config file")
run_parser.add_argument("--name-suffix",
                    dest="name_suffix",
                    default="test_run",
                    help="Suffix to add to pipeline name in order to easily differentiate between different run, defaults to test_run")

clean_parser = subparsers.add_parser("clean", help="Clean after previous run")
clean_parser.add_argument("-v",
                          "--verbose",
                          dest="verbose",
                          action="store_true",
                          default=False,
                          help="Be loud!")
clean_parser.add_argument("-y",
                          "--yes",
                          dest="yes",
                          action="store_true",
                          default=False,
                          help="Force deletion of files.")
clean_parser.add_argument("--make-backup",
                          dest="make_backup",
                          action="store_true",
                          default=False,
                          help="Instead of deleting file with results make its backup")

try:
    options = parser.parse_args()
except Exception, e:
    parser.print_help()


# redefine a functions for writing to stdout and stderr to save some writting
syserr = sys.stderr.write
sysout = sys.stdout.write

# Define directories
working_directory = os.getcwd()
pipeline_directory = os.path.dirname(os.path.abspath(__file__))
output_directory = os.path.join(working_directory, "output")
for_features_directory = os.path.join(working_directory, "output/ForFeatures")
index_directory = os.path.join(working_directory, "index")
plots_directory = os.path.join(working_directory, "Plots")
plexy_directory = os.path.join(working_directory, "Input")


if options.command == 'clean':
    try:
        if options.yes:
            is_sure = "yes"
        else:
            is_sure = raw_input("Do you really want to delete previous run (yes/no)?:  ")
        if is_sure.upper().startswith("Y"):
            dirs_to_delete = [plots_directory,
                              index_directory,
                              plexy_directory,
                              output_directory]
            for directory in dirs_to_delete:
                try:
                    if options.verbose:
                        syserr("Deleting %s\n" % directory)
                    shutil.rmtree(directory)
                except OSError:
                    if options.verbose:
                        syserr(" -> no such a directory: %s\n" % directory)
            files_to_backup = ["results_with_probability_annotated.tab"]
            files_to_delete = ["search_anchors.stats",
                               "accessibility.tab",
                               "snoRNAs.rpkm",
                               "flanks.tab",
                               "snoRNAs.bed",
                               "results_with_probability.bed",
                               "results_with_score_and_rpkm_aggregated.tab",
                               "results_with_probability.tab",
                               "snoRNAs.fa",
                               "raw_reads_results.tab",
                               "search_anchors_for_statistics",
                               "results_with_score.tab",
                               "results_with_score_and_rpkm.tab",
                               "anchors.tab",
                               "for_index.clustered",
                               "for_index.fasta",
                               "mapped_reads_annotated_with_snornas.tab"]
            if options.make_backup:
                for f in files_to_backup:
                    if options.verbose:
                        syserr("Backing up %s\n" % os.path.join(working_directory, f))
                    try:
                        os.rename(os.path.join(working_directory, f), os.path.join(working_directory, f + ".bak"))
                    except OSError, e:
                        if options.verbose:
                            syserr(" -> no such a file: %s\n" % os.path.join(working_directory, f))
            else:
                files_to_delete.extend(files_to_backup)
            for f in files_to_delete:
                if options.verbose:
                    syserr("Removing %s\n" % os.path.join(working_directory, f))
                try:
                    os.remove(os.path.join(working_directory, f))
                except OSError, e:
                    if options.verbose:
                        syserr(" -> no such file: %s\n" % os.path.join(working_directory, f))
            if options.verbose:
                syserr("All output files and directories were cleaned\n")
    except Exception as e:
        syserr(traceback.format_exc())
    finally:
        sys.exit()



def mkdir_p(path_to_dir):
    try:
        os.makedirs(path_to_dir)
    except OSError as e: # Python >2.5
        if e.errno == EEXIST and os.path.isdir(path_to_dir):
            sys.stderr.write("Directory %s already exists. Skipping.\n" % path_to_dir)
        else:
            raise e

settings = ConfigObj(options.config).dict()
mkdir_p(output_directory)
mkdir_p(index_directory)
mkdir_p(plots_directory)
mkdir_p(plexy_directory)
mkdir_p(for_features_directory)

from Jobber import JobClient

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
                                             os.path.join(pipeline_directory, "scripts/rg-generate-input-for-plexy-or-rnasnoop.py"),
                                             settings['general']['snoRNAs'],
                                             plexy_directory,
                                             settings['general']['type'])
generate_plexy_id = jobber.job(generate_plexy_command, {'name': "GenerateInput"})

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
                   os.path.join(working_directory, "snoRNAs.bed"),
                                     settings['general']['type'])

calculate_snorna_expression_command = "python %s --input %s --output %s --library %s --snoRNAs %s --quantile 0.0 --type %s" % calculate_snorna_expression_tuple


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
                                                    ('l', "membycore=%s" % settings['tasks']['ClusterReads'].get('mem_req', '2G'))]
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
                                                            ('l', "membycore=%s" % settings['tasks']['MakeFastaFromClusters'].get('mem_req', '2G'))]
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
                                                            ('l', "membycore=%s" % settings['tasks']['BuildBowtieIndex'].get('mem_req', '2G'))]
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
merge_results_command = "cat {output_dir}/*.scorebed > {cwd}/results_with_score.tab".format(output_dir=output_directory,
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
#                                                                os.path.join(working_directory, "results_with_score.tab"),
#                                                                os.path.join(working_directory, "results_with_score_clustered.tab")
#                                              )
# generate_fasta_id = jobber.job(generate_fasta_command, {'name': "GenerateFasta"})

# Add RPKM to results
#
add_rpkm_command = "python %s --input %s --output %s --rpkm %s --annotated-reads %s --type %s -v" % (
                                             os.path.join(pipeline_directory, "scripts/rg-add-rpkm-to-score.py"),
                                             os.path.join(working_directory, "results_with_score.tab"),
                                             os.path.join(working_directory, "results_with_score_and_rpkm.tab"),
                                             os.path.join(working_directory, "snoRNAs.rpkm"),
                                             os.path.join(working_directory, "mapped_reads_annotated_with_snornas.tab"),
                                             settings['general']['type']
                                             )
add_rpkm_id = jobber.job(add_rpkm_command, {'name': "AddRPKM",
                                            'dependencies': [merge_results_id]})

# Cluster results
#
aggregate_by_site_command = "python %s --input %s --output %s --type %s" % (
                               os.path.join(pipeline_directory, "scripts/rg-aggregate-scored-results.py"),
                                                 os.path.join(working_directory, "results_with_score_and_rpkm.tab"),
                                                 os.path.join(working_directory, "results_with_score_and_rpkm_aggregated.tab"),
                                                 settings['general']['type']
                               )
aggregate_by_site_id = jobber.job(aggregate_by_site_command, {'name': "AggregateBySite",
                                                              'dependencies': [add_rpkm_id]})
#
# Features for model
#

#First step is to split the file with results
split_res_tuple =  (os.path.join(pipeline_directory, "scripts/rg-split-file-into-chunks.py"),
                    os.path.join(working_directory, "results_with_score_and_rpkm_aggregated.tab"),
                    for_features_directory,
                    5000,
                    "part_",
                    ".result"
                    )
split_res_command = "python %s --input %s --dir %s --lines %s --prefix %s --suffix %s -v" % split_res_tuple
split_res_files_id = jobber.job(split_res_command, {'name': "SplitResult",
                                                    'dependencies': [aggregate_by_site_id]})

calculate_features_group_id = jobber.startGroup({'name': "FeaturesGroup",
                                                 'dependencies': [split_res_files_id]})

#We call the script that will generate the jobs that will analyse the split files. We pass the id of the group
#and the folder where the script will find the splitted files.
feat_tuple = (os.path.join(pipeline_directory, "run-features.py"),
                  for_features_directory,
                  calculate_features_group_id,
                  os.path.abspath(options.config),
                  working_directory)
feat_command = "python %s --input-dir %s --group-id %s --config %s --working-dir %s -v" % feat_tuple
jobber.job(feat_command, {'name': "CreateFeaturesJobs"})

jobber.endGroup()


# Calculate Probability

calculate_probability_settings = settings['tasks']['CalculateProbability']
calculate_probability_command = "python %s --input %s --output %s --accessibility %s --flanks %s --model %s" % (
                                             os.path.join(pipeline_directory, "scripts/rg-calculate-probability.py"),
                                             os.path.join(working_directory, 'results_with_score_and_rpkm_aggregated.tab'),
                                             os.path.join(working_directory, 'results_with_probability.tab'),
                                             os.path.join(working_directory, "accessibility.tab"),
                                             os.path.join(working_directory, "flanks.tab"),
                                             settings['general']['model'])
calculate_probability_id = jobber.job(calculate_probability_command, {'name': "CalculateProbability",
                                    'options': [('q', calculate_probability_settings.get('queue', 'short.q')),
                                                ('l', "membycore=%s" % calculate_probability_settings.get('mem_req', '4G'))],
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
                                    'options': [('l', "membycore=4G")],
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
                           settings['general']['annotations_repeats'],
                           )

annotate_results_command = "python %s --input %s --output %s --genes %s --regions %s --repeats %s" % annotate_results_tuple


annotate_results_id = jobber.job(annotate_results_command, {'name': "AnnotateResults",
                                    'options': [('q', annotate_results_settings.get('queue', 'short.q')),
                                                ('l', "membycore=%s" % annotate_results_settings.get('mem_req', '2G'))],
                                        'dependencies': [convert_to_bed_id]})


############################### END PIPELINE GROUP ###################################
jobber.endGroup()

# Before launching we print the command to stop the pipeline
print "In order to stop the pipeline run a command:"
print "jobber_server -command delete -jobId %i, or use web interface" % (pipeline_id)

#You need to always launch, otherwise jobs wont get executed.
jobber.launch(pipeline_id)
