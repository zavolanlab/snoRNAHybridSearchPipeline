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
Usage: $0 -d <string> [-r] [-c] [-p <string>] [-f <string>]

This script will start the run the calculations for snoRNA chimeras for human.

OPTIONS:
   -h                  Show this message.
   -r                  Run test.
   -c                  Run clean up.
   -d                  Absolute path to the data directory that accompanies this repository.
   -p                  Path to PLEXY (how to call plexy.pl script). Defaults to plexy.pl.
   -c                  Path to CONTRAfold (how to call contrafold). Defaults to contrafold.
   -e                  Executer. Defaults to drmaa. Another option is local.
EOF
}

if [ $# -lt 1 ] ; then
    usage
    exit 1
fi

directory=""
clean=""
run=""
plexy="plexy.pl"
contrafold="contrafold"
executer="drmaa"

while getopts "hrcd:p:f:e:" opt
do
   case "${opt}" in
      h) usage;;
      r) run="run";;
      c) clean="clean";;
      d) directory=$OPTARG;;
      p) plexy=$OPTARG;;
      f) contrafold=$OPTARG;;
      f) executer=$OPTARG;;
   esac
done

if [ "$directory" == "" ] ; then
    echo "Please specify path to the data directory you downloaded!"
    usage
    exit 1
fi

cwd=$(pwd)

if [ "$clean" == "clean" ]
then
    echo "Cleaning"
    rm -rf config.ini
    python ../snoRNAHybridSearch.py clean -y
elif [ "$run" == "run" ]
then
    echo "Running test"
    rm -rf config.ini
    python ../snoRNAHybridSearch.py clean -y
    sed "s~<path_to_data_directory>~${directory}~g" config_template.ini > config.ini
    sed -i "s~<PLEXYpath>~${plexy}~g" config.ini
    sed -i "s~<CWD>~${cwd}~g" config.ini
    sed -i "s~<CONTRAfoldpath>~${contrafold}~g" config.ini
    sed -i "s~<DefaultExecuter>~${executer}~g" config.ini
    python ../snoRNAHybridSearch.py run --config config.ini --name-suffix test_run
else
    echo "Please specify if you would like to run or to clean data with -r/-c options."
    exit 1
fi
