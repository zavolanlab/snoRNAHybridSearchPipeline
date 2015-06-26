import collections
import subprocess
import tempfile
import datetime
import logging
import random
import drmaa
import yaml
import time
import stat
import copy
import os
from ruffus.cmdline import MESSAGE


class NoConfigurationFileException(Exception):
    def __init__(self, *errmsg):
        Exception.__init__(self, *errmsg)


class DRMAAException(Exception):
    def __init__(self, *errmsg):
        Exception.__init__(self, *errmsg)

def update_dict(orig_dict, new_dict):
    """
    Recursively merge or update dict-like objects.

    Args:
        orig_dict (dict): dict you want to update
        new_dict (dict): dict used for update

    Returns:
        dict (dict): updated orig_dict

    Original code at: http://stackoverflow.com/questions/3232943/update-value-of-a-nested-dictionary-of-varying-depth

    >>> update({'k1': {'k2': 2}}, {'k1': {'k2': {'k3': 3}}, 'k4': 4})
    {'k1': {'k2': {'k3': 3}}, 'k4': 4}
    """
    for key, val in new_dict.iteritems():
        if isinstance(val, collections.Mapping):
            tmp = update_dict(orig_dict.get(key, { }), val)
            orig_dict[key] = tmp
        elif isinstance(val, (list, long)):
            if orig_dict[key]:
                orig_dict[key] = (orig_dict[key] + val)
            else:
                orig_dict[key] = val
        else:
            orig_dict[key] = new_dict[key]
    return orig_dict

def load_configs(default, updated):
    """
    Load configuration files.
    First default configuration is loaded (if found) and after it is updated with local file.

    Args:
        default (string): path to default configuration file
        updated (string): path to local configuration file

    Returns:
        configuration dictionary

    Raises:
        NoConfigurationFileException
    """
    is_default = True
    is_updated = True
    # load default config file
    try:
        with open(default) as handler:
            yaml_file = yaml.load(handler)
    except EnvironmentError, e:
        is_default = False
        print str(e)
    # load config file for specific run
    try:
        with open(updated) as handler:
            yaml_file_new = yaml.load(handler)
    except EnvironmentError, e:
        is_updated = False
        print str(e)
    if not is_default and not is_updated:
        raise NoConfigurationFileException("No configuration file specified: neither default nor local")

    # update yaml_file
    return update_dict(yaml_file, yaml_file_new)



def get_run_options(params, task_name, other_updates = None):
    """
    Return run options - if they do not exist make default options and
    return it. You can specify also inside the function what you would like update
    in case the values are calculated on flight.

    Args:
        params (dict): run parameters read from configuration file
        task_name (string): name of the task
        other_updates (dict): fresh (on flight inside the task) generated parameters to be used

    Returns:
        default_options (dict): configuration dict for this run
    """
    default_options = {'to_cluster': False,
                       'job_queue': None,
                       'job_priority': None,
                       'job_name': None,
                       'job_other_options': None,
                       'touch_only': None,
                       'expand_statement': "",
                       'expand_after_statement': "",
                       'job_array': None}
    try:
        run_options = copy.deepcopy(params['task_options'][task_name]['run_options'])
        updated = update_dict(default_options, run_options)
        if other_updates:
            return update_dict(updated, other_updates)
        else:
            return updated
    except KeyError:
        return default_options


def read_stdout_stderr_from_files( stdout_path, stderr_path, logger = None, cmd_str = "", tries=5):
    '''
    Reads the contents of two specified paths and returns the strings
    Thanks to paranoia approach contributed by Andreas Heger:
        - Retry just in case file system hasn't committed.
        - Logs error if files are missing: No big deal?
        - Cleans up files afterwards

    Args:
        stdout_path (string): path to the stdout
        stderr_path (string): path to the stderr

    Kwargs:
        logger (logger): logging object to log
        cmd_str (string): command used
        tries (int): number of trials to read the files

    Returns:
        stdout (list): stdout read into list
        stderr (list): stderr read into list
    '''
    #
    #   delay up to 10 seconds until files are ready
    #
    for xxx in range(tries):
        if os.path.exists( stdout_path ) and os.path.exists( stderr_path ):
            break
        time.sleep(2)

    try:
        stdout = open( stdout_path, "r" ).readlines()
    except IOError, msg:
        if logger:
            logger.warning( "could not open stdout: %s for \n%s" % (msg, cmd_str))
        stdout = []

    try:
        stderr = open( stderr_path, "r" ).readlines()
    except IOError, msg:
        if logger:
            logger.warning( "could not open stderr: %s for \n%s" % (msg, cmd_str))
        stderr = []

    #
    #   cleanup ignoring errors
    #
    try:
        os.unlink( stdout_path )
        os.unlink( stderr_path )
    except OSError, msg:
        pass

    return stdout, stderr

def expand_statement_for_run( statement, exec_prefix = "", exec_suffix = "", ignore_pipe_errors = False ):
    '''
    add prefix and suffix to statement.

    Args:
        statement (string): command to use

    Kwargs:
        exec_prefix (string): prefix to add
        exec_suffix (string): suffix to add
        ignore_pipe_errors (bool): ignore pipe errors

    Returns:
        whole_statement (string): extended statement
    '''

    if ignore_pipe_errors:
        return "".join( (exec_prefix, statement) )
    else:
        return "".join( (exec_prefix, statement, exec_suffix) )


def setup_drmaa_job(job_queue = None, job_priority = None, job_name = None,
                    job_other_options = "", logger = None, drmaa_session = None):
    """
    setup drmaa job

    Kwargs:
        job_queue (string): queue for the job
        job_priority (string): queue priority
        job_name (string): name of the job
        job_other_options (string): additional options
        logger (logger): logging object to log
        drmaa_session (Session): DRMAA session object

    Returns:
        job_template (Template): DRMAA job template
    """

    job_template = drmaa_session.createJobTemplate()
    job_template.workingDirectory = os.getcwd()
    #job_template.jobEnvironment = { 'BASH_ENV' : '~/.bashrc' } #  if this is on it somehow causes problems
    job_template.args = []

    # optional job parameters
    job_queue      = "-q " + job_queue            if job_queue                   else ""
    job_priority  = ("-p %d " % job_priority)   if not job_priority is None   else ""
    job_name            = "-N " + job_name                  if job_name                         else ""

    job_template.nativeSpecification = "-V {job_queue} {job_priority} {job_name} {job_other_options}" .format(
                                        job_queue      = job_queue,
                                        job_priority  = job_priority,
                                        job_name            = job_name,
                                        job_other_options   = job_other_options)

    # separate stdout and stderr
    job_template.joinFiles=False

    return job_template


def build_job_script(statement, logger = None, exec_prefix=None, exec_suffix=None):
    '''
    Build script to use with as a job

    Args:
        statement (string): command to use

    Kwargs:
        logger (logger): logging object to log
        exec_prefix (string): prefix to add
        exec_suffix (string): suffix to add

    Returns:
        job_script_path (string): path to job script
        stdout_path (string): path to stdout
        stderr_path (string): path to stderr
    '''
    time_stmp_str = "_".join(map(str, datetime.datetime.now().timetuple()[0:6]))
    tmpfile = tempfile.NamedTemporaryFile(mode='w+b', prefix='drmaa_script_' + time_stmp_str + "_", dir = os.getcwd(),  delete = False)

    tmpfile.write("#!/bin/bash\n" )
    tmpfile.write(expand_statement_for_run(statement=statement, exec_prefix=exec_prefix, exec_suffix=exec_suffix) + "\n" )
    tmpfile.close()

    job_script_path = os.path.abspath( tmpfile.name )
    stdout_path = job_script_path + ".stdout"
    stderr_path = job_script_path + ".stderr"

    os.chmod( job_script_path, stat.S_IRWXG | stat.S_IRWXU )

    return (job_script_path, stdout_path, stderr_path)


def run_job_using_drmaa (statement, job_queue = None, job_priority = None, job_name = None,
                         job_other_options = "", logger = None, drmaa_session = None, expand_statement="",
                         expand_after_statement="", job_array=None):

    """
    Runs specified command remotely using drmaa

    Args:
        statement (string): command to use

    Kwargs:
        job_queue (string): queue for the job
        job_priority (string): queue priority
        job_name (string): name of the job
        job_other_options (string): additional options
        logger (logger): logging object to log
        drmaa_session (Session): DRMAA session object
        expand_statement (string): a string used to expand statement (eg. modules to load)
        job_array (list): list of 1-based integers in the format [start, stop, increment]

    Returns:
        stdout (list): stdout read into list
        stderr (list): stderr read into list

    Raises:
        DRMAAException
    """
    # define some variables to be used after
    exec_prefix = """
    set -e

    handler() {
    echo "Job $SGE_TASK_ID receives signal : $1" 1>&2;
        exit 1;
    }

    trap "handler Usr1" SIGUSR1
    trap "handler Stop" SIGSTOP
    trap "handler Usr2" SIGUSR2
    trap "handler Term" SIGTERM
    trap "handler Quit" SIGQUIT
    trap "handler  Int" SIGINT
    trap "handler  Bus" SIGBUS
    trap "handler  Seg" SIGSEGV
    """
    exec_prefix = exec_prefix.replace("\n    ", "\n")
    exec_prefix += "\n" + expand_statement + "\n"

    exec_suffix = '\n' + expand_after_statement

    #
    #   make job template
    #
    job_template = setup_drmaa_job(job_queue, job_priority, job_name, job_other_options, logger, drmaa_session)

    #
    #   make job script
    #
    job_script_path, stdout_path, stderr_path = build_job_script(statement, logger, exec_prefix, exec_suffix)
    job_template.remoteCommand  = job_script_path
    job_template.outputPath     = ":"+ stdout_path
    job_template.errorPath      = ":" + stderr_path


    #
    #   Run job and wait
    #
    if job_array:
            start, end, increment = job_array
            if logger:
                logger.log(MESSAGE, "starting an array job: %i-%i,%i" % (start, end, increment ))
            # sge works with 1-based, closed intervals
            jobids = drmaa_session.runBulkJobs( job_template, start, end, increment )
            if logger:
                logger.log(MESSAGE, "%i array jobs have been submitted as jobid %s" % (len(jobids), jobids[0]) )
            retval = drmaa_session.synchronize(jobids, drmaa.Session.TIMEOUT_WAIT_FOREVER, True)
    else:
        jobid = drmaa_session.runJob(job_template)
        if logger:
            logger.log(MESSAGE, "job has been submitted with jobid %s" % str(jobid ))

        try:
            retval = drmaa_session.wait(jobid, drmaa.Session.TIMEOUT_WAIT_FOREVER)
        except Exception, msg:
            # ignore message 24 in PBS
            # code 24: drmaa: Job finished but resource usage information and/or termination status could not be provided.":
            if not msg.message.startswith("code 24"): raise
            retval = None


    #
    #   Read output
    #
    stdout, stderr = read_stdout_stderr_from_files( stdout_path, stderr_path, logger, statement)

    #
    #   Throw if failed
    #
    if not job_array:
        if retval and retval.exitStatus != 0:
            logger.log(logging.ERROR, "The drmaa command was terminated by signal %i:\n"
                                   "The original command was:\n%s\n"
                                   "The stderr was: \n%s\n" %
                                     (retval.exitStatus, statement, "".join( stderr)) )
            raise DRMAAException( "The drmaa command was terminated by signal %i:\n"
                                   "The original command was:\n%s\n"
                                   "The stderr was: \n%s\n" %
                                     (retval.exitStatus, statement, "".join( stderr)) )


    #
    #   clean up job template
    #
    drmaa_session.deleteJobTemplate(job_template)

    #
    #   Cleanup job script
    #
    try:
        os.unlink( job_script_path )
        # pass
    except OSError:
        if logger:
            logger.warning( "Temporary job script wrapper '%s' missing (and ignored) at clean-up" % job_script_path )

    return stdout, stderr


def run_job_locally (statement, logger = None):
    """
    Runs specified command locally instead of drmaa
    """
    process = subprocess.Popen(  statement,
                                 cwd = os.getcwd(),
                                 shell = True,
                                 stdin = subprocess.PIPE,
                                 stdout = subprocess.PIPE,
                                 stderr = subprocess.PIPE )

    stdout, stderr = process.communicate()

    if process.returncode != 0:
        raise DRMAAException("The locally run command was terminated by signal %i:\n"
                               "The original command was:\n%s\n"
                               "The stderr was: \n%s\n" %
                                 (-process.returncode, statement, "".join( stderr)))

    return stdout.splitlines(True), stderr.splitlines(True)


def run(statement, logger=None, session=None, randseed=None, **kwargs):
    """
    Runs specified command either using drmaa

    Args:
        statement (string): command to use

    Kwargs:
        logger (logger): logging object to log
        kwargs (dict): set of parameters important to run
    """

    if not kwargs["to_cluster"]:
        return run_job_locally(statement, logger)
    else:
        if session is None:
            raise DRMAAException("Provide shared session for the jobs")

        if randseed:
            random.seed(randseed)
        else:
            random.seed()
        time.sleep(random.randint(0, 10))
        return run_job_using_drmaa(statement=statement,
                                   job_queue=kwargs["job_queue"],
                                   job_priority=kwargs["job_priority"],
                                   job_name=kwargs["job_name"],
                                   job_other_options=kwargs["job_other_options"],
                                   job_array=kwargs["job_array"],
                                   logger=logger,
                                   drmaa_session=session,
                                   expand_statement=kwargs["expand_statement"],
                                   expand_after_statement=kwargs["expand_after_statement"])

