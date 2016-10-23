#! /bin/sh
#
# run_test.sh
# Copyright (C) 2016 Rafal Gumienny <r.gumienny@unibas.ch>
#
# Distributed under terms of the GPL license.
#


# end when error
set -e
# raise error when variable is unset
set -u
# raise error when in pipe
set -o pipefail



usage()
{
cat << EOF
Usage: $0 -d <string> [-r] [-c] [-p <string>] [-f <string>] [-m "modA modB modC"]

This script will start the run the calculations for snoRNA chimeras for human.

OPTIONS:
   -h                  Show this message.
   -r                  Run test.
   -c                  Run clean up.
   -d                  Absolute path to the data directory that accompanies this repository.
   -p                  Path to PLEXY (how to call plexy.pl script). Defaults to plexy.pl.
   -c                  Path to CONTRAfold (how to call contrafold). Defaults to contrafold.
   -e                  Executer. Defaults to drmaa. Another option is local.
   -m                  Modules. A list of modules to load (if HPC or environment requires).

Note on modules. I have used following modules on our HPC environment in order to
fulfill dependencies:
    - Python/2.7.5-goolf-1.4.10
    - GCC/4.7.2
    - DRMAA/0.7.6-goolf-1.4.10-Python-2.7.5
    - HTSeq/0.6.1p1-goolf-1.4.10-Python-2.7.5
    - Bowtie2/2.2.6-goolf-1.4.10
    - OpenBLAS/0.2.6-gompi-1.4.10-LAPACK-3.4.2
    - BEDTools/2.25.0-goolf-1.4.10
    - CONTRAfold/2.02-goolf-1.4.10
    - ViennaRNA/2.1.8-goolf-1.4.10
    - SAMtools/1.2-goolf-1.4.10

EOF
}

directory=""
clean=""
run=""
plexy="plexy.pl"
contrafold="contrafold"
executer="drmaa"
modules=""

while getopts "hrcd:p:f:e:m:" opt
do
   case "${opt}" in
      h) usage; exit 1;;
      r) run="run";;
      c) clean="clean";;
      d) directory=$OPTARG;;
      p) plexy=$OPTARG;;
      f) contrafold=$OPTARG;;
      e) executer=$OPTARG;;
      m) modules=$OPTARG;;
   esac
done


if [ "$clean" == "clean" ]
then
    echo "####################### CLEANING ###########################"
    echo
    rm -rf config.ini
    python ../snoRNAHybridSearch.py clean -y -v
    echo
    echo "############################################################"
    exit 0
fi

cwd=$(pwd)

if [ "$run" == "run" ]
then
    if [ "$directory" == "" ] ; then
        echo "####################### ERROR ##############################"
        echo
        echo "Please specify path to the data directory you downloaded!"
        echo
        echo "#############################################################"
        echo
        usage
        exit 1
    fi
    echo "####################### RUNNING ###########################"
    echo
    rm -rf config.ini
    python ../snoRNAHybridSearch.py clean -y
    sed "s~<path_to_data_directory>~${directory}~g" config_template.ini > config.ini
    sed -i "s~<PLEXYpath>~${plexy}~g" config.ini
    sed -i "s~<CWD>~${cwd}~g" config.ini
    sed -i "s~<CONTRAfoldpath>~${contrafold}~g" config.ini
    sed -i "s~<DefaultExecuter>~${executer}~g" config.ini
    python ../snoRNAHybridSearch.py run --config config.ini --name-suffix test_run --modules ${modules} -v
    echo
    echo "############################################################"
fi
