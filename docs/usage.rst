Usage
*****

Basic usage
===========

Command to launch the pipeline is as follows:

.. code-block:: bash

    python snoRNAHybridSearch.py run --config congig.ini --name-suffix name_of_the_run

All parameters for the script:

.. argparse::
    :ref: snoRNAHybridSearch.parser
    :prog: snoRNAHybridSearch

Preparing config file and input files
=====================================

Copy config_example.ini from snoRNAHybridSearchPipeline directory to your working directory (directory
where you want to perform calculation, WD):

.. code-block:: bash

    cd Your/Working/Direcory
    cp Path/To/snoRNAHybridSearchPipeline/config_example.ini config.ini

Set all the necessary paths in your config.ini file as indicated in the comments inside the file. The most importand are:
 * **unmapped_reads**: "Absolute/Path/To/unmapped_reads.fa" - an abs path to an input FASTA file with sequences
   that were unmapped in sequencing experiment - see the example file in additional data.
 * **bed_for_index**: "Absolute/Path/To/mapped_reads.bed" - abs path to a BED file with the positions of
   mapped reads in the experiment - see the example file in additional data.
 * **PLEXY_bin**: "Absolute/Path/To/plexy.pl" - path to PLEXY binary (or how you invoke it in the bash)
 * **contrafold_binary**: "contrafold" - path to CONTRAfold binary (or how you invoke in the bash)

Model path:
 * **model**: "Path/To/snoRNAHybridSearch/data/model.bin" - abs path to the model used to calculate
   probability (you can find it in the pipeline directory named model.bin)

snoRNA table:
 * **snoRNAs**: "Absolute/Path/To/snoRNAs_table.tab" - abs path to the table containing all the necessary
   information abut snoRNAs. This table is provided with
   pipeline additional data and for human it is located in the snoRNAHybridSearchData/human/snoRNAs/snoRNAs.tab.
   We have also prepared the table for mouse located in the snoRNAHybridSearchData/mouse/snoRNAs/snoRNAs_table.tab.
   You can also prepare your own snoRNA input - please follow the conventions in the table and pay
   attention to columns described in the README file.

Additional "chromosomes":
 * This files has to be also split into separate FASTA sequences and those sequences has
   to be put into directory with genome. By default, genome directory that can be downloaded
   additionally contains these sequences already prepared.
 * **rRNAs**: "Absolute/Path/To/rRNAs.fa" # rRNA sequences. This is provided with the pipeline
   in data directory, although own can be used. The location for human is
   snoRNAHybridSearchData/human/TargetRNAs/rRNAs_hsa.fa and for mouse snoRNAHybridSearchData/mouse/rRNAs_mmu.fa.
 * **tRNAs**: "Absolute/Path/To/tRNAs.fa" # tRNA sequences. This is provided with
   the pipeline in data directory, although own can be used. The location for human is
   snoRNAHybridSearchData/human/TargetRNAs/tRNAs_hsa.fa and for mouse snoRNAHybridSearchData/mouse/tRNAs_mmu.fa
 * **snRNAs**: "Absolute/Path/To/snRNAs.fa" # snRNA sequences. This is provided with
   the pipeline in data directory, although own can be used. The location for human is
   snoRNAHybridSearchData/human/TargetRNAs/snRNAs_hsa.fa and for mouse snoRNAHybridSearchData/mouse/TargetRNAs/snRNAs_mmu.fa

Annotation files:
 * Annotation files are used to annotate found target positions. They are generated from corresponding ENSEMBL/GENECODE
   gff3 files or downloaded from NCBI. These files can be found in the annotations subdirectory in given species data directory.
 * **annotations_genes**: "Absolute/Path/To/Annotations/genes.gff3". This file is generated from ENSEMBL/GENECODE file and contains information
   abut genes - not transcripts::

    1   pseudogene  gene    11869   14412   .   +   .   gene_id "ENSG00000223972"; gene_name "DDX11L1"; gene_source "ensembl_havana"; gene_biotype "pseudogene";
    1   pseudogene  gene    14363   29806   .   -   .   gene_id "ENSG00000227232"; gene_name "WASH7P"; gene_source "ensembl_havana"; gene_biotype "pseudogene";
    1   lincRNA gene    29554   31109   .   +   .   gene_id "ENSG00000243485"; gene_name "MIR1302-10"; gene_source "ensembl_havana"; gene_biotype "lincRNA";
    1   lincRNA gene    34554   36081   .   -   .   gene_id "ENSG00000237613"; gene_name "FAM138A"; gene_source "ensembl_havana"; gene_biotype "lincRNA";
    1   pseudogene  gene    52473   54936   .   +   .   gene_id "ENSG00000268020"; gene_name "OR4G4P"; gene_source "ensembl_havana"; gene_biotype "pseudogene";
    1   pseudogene  gene    62948   63887   .   +   .   gene_id "ENSG00000240361"; gene_name "OR4G11P"; gene_source "havana"; gene_biotype "pseudogene";
    1   protein_coding  gene    69091   70008   .   +   .   gene_id "ENSG00000186092"; gene_name "OR4F5"; gene_source "ensembl_havana"; gene_biotype "protein_coding";
    1   lincRNA gene    89295   133566  .   -   .   gene_id "ENSG00000238009"; gene_name "RP11-34P13.7"; gene_source "havana"; gene_biotype "lincRNA";
    1   lincRNA gene    89551   91105   .   -   .   gene_id "ENSG00000239945"; gene_name "RP11-34P13.8"; gene_source "havana"; gene_biotype "lincRNA";
    1   pseudogene  gene    131025  134836  .   +   .   gene_id "ENSG00000233750"; gene_name "CICP27"; gene_source "havana"; gene_biotype "pseudogene";

 * **annotations_regions**: "Absolute/Path/To/Annotations/regions.gff3". This file is generated from ENSEMBL/GENECODE file and contains information
   abut the regions in the genes and transcripts like introns, exons, and UTRS::

    1   ensembl_havana  exon    69091   70008   .   +   .   Parent=mRNA_ENST00000335137
    1   ensembl_havana  CDS 69091   70008   .   +   .   Parent=mRNA_ENST00000335137
    1   ensembl exon    134901  135802  .   -   .   Parent=mRNA_ENST00000423372
    1   ensembl intron  135803  137620  .   -   .   Parent=mRNA_ENST00000423372
    1   ensembl exon    137621  139379  .   -   .   Parent=mRNA_ENST00000423372
    1   ensembl three_prime_UTR 134901  135802  .   -   .   Parent=mRNA_ENST00000423372
    1   ensembl three_prime_UTR 137621  138529  .   -   .   Parent=mRNA_ENST00000423372
    1   ensembl CDS 138530  139309  .   -   .   Parent=mRNA_ENST00000423372
    1   ensembl five_prime_UTR  139310  139379  .   -   .   Parent=mRNA_ENST00000423372
    1   ensembl_havana  exon    367640  368634  .   +   .   Parent=mRNA_ENST00000426406

 * **annotations_repeats**: "Absolute/Path/To/Annotations/repeats.gtf". It is a file downloaded from NCBI table browser::

    chr1    hg19_rmsk   exon    16777161    16777470    2147.000000 +   .   gene_id "AluSp"; transcript_id "AluSp";
    chr1    hg19_rmsk   exon    25165801    25166089    2626.000000 -   .   gene_id "AluY"; transcript_id "AluY";
    chr1    hg19_rmsk   exon    33553607    33554646    626.000000  +   .   gene_id "L2b"; transcript_id "L2b";
    chr1    hg19_rmsk   exon    50330064    50332153    12545.000000    +   .   gene_id "L1PA10"; transcript_id "L1PA10";
    chr1    hg19_rmsk   exon    58720068    58720973    8050.000000 -   .   gene_id "L1PA2"; transcript_id "L1PA2";
    chr1    hg19_rmsk   exon    75496181    75498100    10586.000000    +   .   gene_id "L1MB7"; transcript_id "L1MB7";
    chr1    hg19_rmsk   exon    83886031    83886750    980.000000  -   .   gene_id "ERVL-E-int"; transcript_id "ERVL-E-int";
    chr1    hg19_rmsk   exon    100662896   100663391   1422.000000 -   .   gene_id "L2a"; transcript_id "L2a";
    chr1    hg19_rmsk   exon    117440427   117440514   532.000000  +   .   gene_id "L1ME1"; transcript_id "L1ME1";
    chr1    hg19_rmsk   exon    117440495   117441457   4025.000000 +   .   gene_id "L1ME1"; transcript_id "L1ME1_dup1";

Please refere to Annotations/README file for more details on how to generate these files.

Others:
 * **reads_per_file**: number of reads in the split files
 * **anchor_length**: the lenght of the "seed" prepared from snoRNAs which will be searched initially in the unmapped sequences
 * If you would like to run it on cluster follow instructions in the configuration file and ask your admin what parameters you need to set
   up before (like DRMAA path, modules necessary, queues names etc.). All these parameters can be set up in config.ini. To run it locally it
   might take substantial amount of time to perform all calculations.


Example
=======

To test the pipeline go to the tests directory and run:

.. code-block:: bash

    cd Path/To/snoRNAHybridSearch/tests
    bash run_test.sh -h

.. note::

    Usage: ./run_test.sh -d <string> [-r] [-c] [-p <string>] [-f <string>]

    This script will start the run the calculations for snoRNA chimeras for human.

    OPTIONS:
       -h                  Show this message.
       -r                  Run test.
       -c                  Run clean up.
       -d                  Absolute path to the data directory that accompanies this repository.
       -p                  Path to PLEXY (how to call plexy.pl script). Defaults to plexy.pl.
       -c                  Path to CONTRAfold (how to call contrafold). Defaults to contrafold.
       -e                  Executer. Defaults to drmaa. Another option is local.

And if you have installed all the dependancies to default locations (PLEXY, CONTRAfold etc.) run:

.. code-block:: bash

    bash run_test.sh -d /Absolute/Path/To/snoRNAHybridSearchData -r
