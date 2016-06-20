#! /bin/sh
#
# run.sh
# Copyright (C) 2015 Rafal Gumienny <r.gumienny@unibas.ch>
#
# Distributed under terms of the GPL license.
#


# end when error
set -e
# raise error when variable is unset
set -u
# raise error when in pipe
set -o pipefail


for i in 1 2 3 4 5 6
do
	python test_MetaProfile.py &> ${i}.log &
done
