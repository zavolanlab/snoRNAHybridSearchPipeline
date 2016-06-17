#This script generates the jobs of the pipeline that will do the analysis on the splitted files.
import glob
import os
import sys
from configobj import ConfigObj
from jinja2 import Template
from Jobber import JobClient
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

def main(options):
    settings = ConfigObj(options.config).dict()

    working_directory = options.working_dir
    index_directory = os.path.join(working_directory, 'index')
    output_directory = options.input_dir
    pipeline_directory = os.path.dirname(os.path.abspath(__file__))


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

    contrafold_group = jobber.startGroup({'name': 'CONTRAfold'})
    for f in glob.glob(options.input_dir + "/*.result"):
        input_name = os.path.splitext(f)[0]
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
            copy_files = {f: 'input'}
            moveback = {'output': os.path.join(output_directory, input_name + ".accessibility")}

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
                                               'coords': f,
                                               'output': os.path.join(output_directory, input_name + ".accessibility"),
                                               'genome_dir': settings['general']['genome'],
                                               })
        calculate_contrafold_id = jobber.job(calculate_contrafold_command,
                               {'name': 'CalculateCONTRAFold',
                                'uniqueId': True,
                                'options': [('q', calculate_contrafold_settings.get('queue', 'short.q')),
                                            ('l', "membycore=%s" % calculate_contrafold_settings.get('mem_req', '2G'))]
                                })
    jobber.endGroup()


    #
    # Calulate Flanks composition

    flanks_group = jobber.startGroup({'name': 'Flanks'})
    for f in glob.glob(options.input_dir + "/*.result"):
        input_name = os.path.splitext(f)[0]
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
            copy_files = {f: 'input'}
            moveback = {'output': os.path.join(output_directory, input_name + ".flanks")}

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
                                               'coords': f,
                                               'output': os.path.join(output_directory, input_name + ".flanks"),
                                               'genome_dir': settings['general']['genome'],
                                               })
        calculate_flanks_id = jobber.job(calculate_flanks_command,
                               {'name': 'CalculateFlanks',
                                'uniqueId': True,
                                'options': [('q', calculate_flanks_settings.get('queue', 'short.q')),
                                            ('l', "membycore=%s" % calculate_flanks_settings.get('mem_req', '2G'))]
                                })

    jobber.endGroup()

    # We merge the files into our result file after analysis finishes
    merge_contrafold_command = "cat {output_dir}/*.accessibility > {cwd}/accessibility.tab".format(output_dir=output_directory,
                                                                                                       cwd=working_directory)
    merge_contrafold_id = jobber.job(merge_contrafold_command, {'name': "MergeCONTRAfold",
                                     'dependencies': [contrafold_group]})


    # We merge the files into our result file after analysis finishes
    merge_flanks_command = "cat {output_dir}/*.flanks > {cwd}/flanks.tab".format(output_dir=output_directory,
                                                                                 cwd=working_directory)
    merge_flanks_id = jobber.job(merge_flanks_command, {'name': "MergeFlanks",
                                     'dependencies': [flanks_group]})


    jobber.endGroup()
    jobber.launch(options.group_id)
if __name__ == '__main__':
    try:
        options = parser.parse_args()
    except Exception, e:
        parser.print_help()
        sys.exit()
    main(options)
