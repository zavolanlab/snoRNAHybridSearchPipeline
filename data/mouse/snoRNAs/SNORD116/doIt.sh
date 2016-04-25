#! /bin/sh
#
# doIt.sh
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


grep SNORD116 ../../mmu_cd_box_snoRNAs_mm10.tab | ruby -ane 'puts ">#{$F[3]}\n#{$F[6]}"' > SNORD116_snOPY.fa
# we know that SNORD116 in mouse is located on chr7 so we check only those snoRNAs
grep snoRNA ../../annotations/genes.gff3 | grep "^7" | rg-gtf-to-bed.py --chrom-out NCBI | rg-extract-sequences.py --format fasta --genome-dir ../../mm10_ensembl/ > chr7_snoRNAs.fa

python rg-make-alignment-matrix.py --a ./SNORD116_snOPY.fa --b ./chr7_snoRNAs.fa > snoRNAs_from_matrix.tab
cut -f1-6 snoRNAs_from_matrix.tab | ruby -ane 'puts [$F[0], $F[1], $F[2], "snoRNA:SNORD116_"+$F[3].split("_")[0], $F[4], $F[5]].join("\t")' > snoRNAs_ENSEMBL.bed
cat ./snoRNAs_ENSEMBL.bed | rg-extract-sequences.py --format fasta --genome-dir ../../mm10_ensembl/ > snoRNAs_ENSEMBL.fa
