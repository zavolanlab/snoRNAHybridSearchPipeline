Pipeline flow
*************

CD-box snoRNAs
==============

1. Split the input
------------------

At first split the input unmapped sequences into manageable chunks.

.. argparse::
    :ref: scripts.rg_split_fasta.parser
    :prog: rg_split_fasta

2. Generate various files from snoRNAs
--------------------------------------

i. Make FASTA
^^^^^^^^^^^^^
  .. argparse::
      :ref: scripts.rg_generate_fasta.parser
      :prog: rg_generate_fasta

ii. Generate separate files
^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. argparse::
    :ref: scripts.rg_generate_input_for_plexy_or_rnasnoop.parser
    :prog: rg_generate_input_for_plexy_or_rnasnoop


iii. Make BED
^^^^^^^^^^^^^

.. argparse::
    :ref: scripts.rg_generate_snoRNA_bed.parser
    :prog: rg_generate_snoRNA_bed

3. Annotate with snoRNAs
------------------------

Annotate input BED file used for generation of clusters
with snoRNAs.

.. argparse::
    :ref: scripts.rg_annotate_bed.parser
    :prog: rg_annotate_bed


4. Calculate snoRNA expression
------------------------------

.. argparse::
    :ref: scripts.rg_calculate_snoRNA_RPKM.parser
    :prog: rg_calculate_snoRNA_RPKM

5. Prepare anchors
------------------

.. argparse::
    :ref: scripts.rg_prepare_anchors.parser
    :prog: rg_prepare_anchors


6. Build Bowtie2 index
----------------------

i. Cluster reads
^^^^^^^^^^^^^^^^

.. argparse::
    :ref: scripts.rg_cluster_reads.parser
    :prog: rg_cluster_reads

ii. Make FASTA
^^^^^^^^^^^^^^

Prepare FASTA file from clustered reads

.. argparse::
    :ref: scripts.rg_extract_sequences.parser
    :prog: rg_extract_sequences

iii. Build index
^^^^^^^^^^^^^^^^

The index is build with following command:

.. code-block:: bash

    bowtie2-build input.fa path/to/index/bowtie_index 2> /dev/null


7. Run analysis
---------------

For each part split in first task an analysis is run.

i. Search anchors
^^^^^^^^^^^^^^^^^

.. argparse::
    :ref: scripts.rg_search_anchor_and_make_alignments.parser
    :prog: rg_search_anchor_and_make_alignments

ii. Make statistics
^^^^^^^^^^^^^^^^^^^

This is set of two tasks:
 a. Merging the files from anchor search
 b. Making statistics with following script:

.. argparse::
    :ref: scripts.rg_make_stats_for_search.parser
    :prog: rg_make_stats_for_search


iii. Convert to FASTA
^^^^^^^^^^^^^^^^^^^^^

.. argparse::
    :ref: scripts.rg_convert_tab_to_fasta.parser
    :prog: rg_convert_tab_to_fasta


iv. Map reads
^^^^^^^^^^^^^

Map target parts to the cluster with following command:

.. code-block:: bash

    bowtie2 -x ./index/bowtie_index -f -D100 -L 13 -i C,1 --local -k 10 -U input.anchorfasta -S output.sam

v. Convert result to BED
^^^^^^^^^^^^^^^^^^^^^^^^

Convert result from mapping into BED file with following command:

.. code-block:: bash

    samtools view -S input.sam -b -u | bamToBed -tag AS | grep -P "\t\+" > output


vi. Filter BED
^^^^^^^^^^^^^^

.. argparse::
    :ref: scripts.rg_filter_bed.parser
    :prog: rg_filter_bed

vi. Reasign chromosome
^^^^^^^^^^^^^^^^^^^^^^

.. argparse::
    :ref: scripts.rg_get_true_chromosome_positions.parser
    :prog: rg_get_true_chromosome_positions

vii. Append sequence
^^^^^^^^^^^^^^^^^^^^

The same script as for the FASTA extraction from Bowtie2 index.

viii. Calculate PLEXY
^^^^^^^^^^^^^^^^^^^^^

.. argparse::
    :ref: scripts.rg_check_hybrids_with_plexy.parser
    :prog: rg_check_hybrids_with_plexy

ix. Calculate RNAduplex
^^^^^^^^^^^^^^^^^^^^^^^

.. argparse::
    :ref: scripts.rg_check_hybrids_with_rnaduplex.parser
    :prog: rg_check_hybrids_with_rnaduplex


8. Analyse RNAduplex results
----------------------------

RNAduplex and PLEXY results goes slightly different analysis.

i. Merge results
^^^^^^^^^^^^^^^^

Nothing to add

ii. Cluster results
^^^^^^^^^^^^^^^^^^^

.. argparse::
    :ref: scripts.rg_cluster_results.parser
    :prog: rg_cluster_results

iii. Annotate results
^^^^^^^^^^^^^^^^^^^^^

.. argparse::
    :ref: scripts.rg_annotate_positions.parser
    :prog: rg_annotate_positions


iv. Make statistics
^^^^^^^^^^^^^^^^^^^

.. argparse::
    :ref: scripts.rg_make_plots_for_rnaduplex.parser
    :prog: rg_make_plots_for_rnaduplex


9. Analyse PLEXY
----------------

i. Merge results
^^^^^^^^^^^^^^^^

.. code-block:: bash

    cat output/*.scorebed > results_with_score.tab

ii. Merge raw results
^^^^^^^^^^^^^^^^^^^^^

.. code-block:: bash

    cat output/*.truechrombed > raw_reds_results.tab

iii. Append RPKM
^^^^^^^^^^^^^^^^

.. argparse::
    :ref: scripts.rg_add_rpkm_to_score.parser
    :prog: rg_add_rpkm_to_score


iv. Aggregate results by site
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. argparse::
    :ref: scripts.rg_aggregate_scored_results.parser
    :prog: rg_aggregate_scored_results


v. Calculate features
^^^^^^^^^^^^^^^^^^^^^

For each of the site calculate features: accessibility and flanks composition. The PLEXY
is already calculated.

vi. Calculate probability
^^^^^^^^^^^^^^^^^^^^^^^^^

.. argparse::
    :ref: scripts.rg_calculate_probability.parser
    :prog: rg_calculate_probability

vii. Make plots
^^^^^^^^^^^^^^^

.. argparse::
    :ref: scripts.rg_make_stats_for_results.parser
    :prog: rg_make_stats_for_results

viii. Convert to BED
^^^^^^^^^^^^^^^^^^^^

.. argparse::
    :ref: scripts.rg_convert_to_bed.parser
    :prog: rg_convert_to_bed

ix. Annotate results
^^^^^^^^^^^^^^^^^^^^

.. argparse::
    :ref: scripts.rg_annotate_positions.parser
    :prog: rg_annotate_positions


Miscellaneous
=============

Those scripts are not used (yet) or are used to calculate HACA-box
snoRNAs chimeras. For the sake of documentation they are placed here.

.. argparse::
    :ref: scripts.rg_annotate_results_bed.parser
    :prog: rg_annotate_results_bed

.. argparse::
    :ref: scripts.rg_append_genes_and_names.parser
    :prog: rg_append_genes_and_names

.. argparse::
    :ref: scripts.rg_check_hybrids_with_rnasnoop.parser
    :prog: rg_check_hybrids_with_rnasnoop

.. argparse::
    :ref: scripts.rg_compare_results_to_original.parser
    :prog: rg_compare_results_to_original

.. argparse::
    :ref: scripts.rg_convert_to_asmbed.parser
    :prog: rg_convert_to_asmbed

.. argparse::
    :ref: scripts.rg_convert_to_coords.parser
    :prog: rg_convert_to_coords

.. argparse::
    :ref: scripts.rg_convert_unmapped_to_fasta.parser
    :prog: rg_convert_unmapped_to_fasta

.. argparse::
    :ref: scripts.rg_correlate_expression_with_hybrids.parser
    :prog: rg_correlate_expression_with_hybrids

.. argparse::
    :ref: scripts.rg_filter_reads_for_clustering.parser
    :prog: rg_filter_reads_for_clustering

.. argparse::
    :ref: scripts.rg_generate_haca_stems_for_rnasnoop.parser
    :prog: rg_generate_haca_stems_for_rnasnoop

.. argparse::
    :ref: scripts.rg_get_search_info.parser
    :prog: rg_get_search_info

.. argparse::
    :ref: scripts.rg_get_snoRNA_gff.parser
    :prog: rg_get_snoRNA_gff

.. argparse::
    :ref: scripts.rg_make_cd_snoRNAs_families.parser
    :prog: rg_make_cd_snoRNAs_families

.. argparse::
    :ref: scripts.rg_shuffle_fasta_sequences.parser
    :prog: rg_shuffle_fasta_sequences

.. argparse::
    :ref: scripts.rg_split_file_into_chunks.parser
    :prog: rg_split_file_into_chunks
