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

Preparing input files
=====================

Preparing config file
=====================

Copy config_example.ini from MIRZA-G directory to your working directory (directory
where you want to perform calculation, WD):

.. code-block:: bash

    cd Your/Working/Direcory
    cp Path/To/snoRNAHybridSearch/config_example.ini config.ini

Set all the necessary paths in your config.ini file as indicated in the comments inside the file. The most importand are:
 * **unmapped_reads**: "Path/To/unmapped_reads.fa" - an abs path to an input FASTA file with sequences that were unmapped in sequencing experiment
 * **bed_for_index**: "Path/To/mapped_reads.bed" - abs path to a BED file with the positions of mapped reads in the experiment
 * **PLEXY_bin**: "Path/To/plexy.pl" - path to PLEXY binary (or how you invoke it in the bash)

Model path:
 * **model**: "Path/To/snoRNAHybridSearch/data/model.bin" - abs path to the model used to calculate probability (you can find it in the pipeline/data directory)

Additionally when you would like to calculate with evolutionary conservation you have to make sure that the variable run_only_MIRZA in CalculateMIRZA task is set to “no” instead of “yes” and that you provide proper paths with aligned UTRs and evolutionary tree:
 * **phylogenetic_tree**: "Path/To/MIRZA-G/data/human_tree.nh" - abspath to provided phylogenetic tree
 * **alignment_directory**: "Path/To/MIRZA-G/data/HumanAlignments/" - abspath to provided human alignments directory. If you downloaded package from CLIPz website
      this directory is already in the MIRZA-G directory. If you downloaded from GitHub you have to download it additionally.

If you would like to run it on cluster follow instructions in the configuration file and ask your admin what parameters you need to set
up before (like DRMAA path, modules necessary, queues names etc.). All these parameters can be set up in config.ini.

To run it locally it takes ~70 to 90 seconds for one miRNA without conservation calculation and ~170 seconds with calculation (This
might be substantial amount of time (up to half an hour per miRNA) for worse processors).


Example
=======

To test the pipeline go to the tests directory and run:

.. code-block:: bash

    cd Path/To/snoRNAHybridSearch/tests
    bash rg_run_test.sh help

.. note::

    Usage: rg_run_test.sh clean/run ['CONTRAfold/binary/path']

And if you have installed MIRZA and CONTRAfold to default locations (MIRZA and contrafold) run:

.. code-block:: bash

    bash rg_run_test.sh run
