.. _tutorial_crispr_screening_PRJNA358686:


===============================================================================
CROP-Seq; 1:1:1 Mixture of DNMT3B, MBD1, and TET2 Knockout Cell Lines (HEK293T)
===============================================================================

Dataset: CROP-Seq; 1:1:1 Mixture of DNMT3B, MBD1, and TET2 Knockout Cell Lines (HEK293T)

Datlinger, P., Rendeiro, A.F., Schmidl, C., Krausgruber, T., Traxler, P., Klughammer, J., Schuster, L.C., Kuchler, A., Alpar, D., and Bock, C. (2017). `Pooled CRISPR screening with single-cell transcriptome readout`_. Nat. Methods *14*, 297–301.

.. _`Pooled CRISPR screening with single-cell transcriptome readout`: https://doi.org/10.1038/nmeth.4177

|


Preparation
===========

Download fastq files from `European Nucleotide Archive`_.

.. _`European Nucleotide Archive`: https://www.ebi.ac.uk/ena/browser/view/PRJNA521522?show=reads

.. code-block:: console

    $ curl -O ftp.sra.ebi.ac.uk/vol1/fastq/SRR516/009/SRR5163029/SRR5163029_1.fastq.gz

    $ curl -O ftp.sra.ebi.ac.uk/vol1/fastq/SRR516/009/SRR5163029/SRR5163029_2.fastq.gz

|


Mapping
=======

This dataset is single-cell RNA-Seq on HEK293T cell lines mixed in equal proportions knocked out for DNMT3B, MBD1, and TET2 (`Fig. 1h`_). The platform of this dataset is Drop-Seq. The details about the original data processing can be found here_. Briefly, raw data were process with Drop-seq Tools v1.12 software. The first 12 bases on read 1 are cell barcodes, followed by 8 bases UMIs. Captured single cell transcripts are on read 2.

.. _`Fig. 1h`: https://www.nature.com/articles/nmeth.4177/figures/1

.. _here: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM2450588

I re-processed these raw reads to get cell-associated barcodes for downstream analysis. Briefly, raw reads were trimmed using `Trim Galore!`_ and mapped to the human reference genome `refdata-gex-GRCh38-2020-A`_ using STAR v2.7.8a. Plasmid `CROPseq-Guide-Puro`_ sequence was not included in the reference.

.. _`Trim Galore!`: https://www.bioinformatics.babraham.ac.uk/projects/trim_galore

.. _`refdata-gex-GRCh38-2020-A`: https://support.10xgenomics.com/single-cell-gene-expression/software/release-notes/build#GRCh38_2020A

.. _`CROPseq-Guide-Puro`: https://www.addgene.org/86708/

.. csv-table:: Mapping statistics
   :widths: 70, 20
   :header-rows: 0

    "Number of Read Pairs","227,621,653"
    "Number of Read Pairs, After Trimming","201,527,916"
    "Reads Mapped to Genome: Unique+Multiple","0.7654"
    "Reads Mapped to Genome: Unique","0.696169"
    "Estimated Number of Cells","1,199"
    "Fraction of Reads in Cells","0.479281"
    "Mean Reads per Cell","45,290"
    "Median Reads per Cell","16,329"
    "Mean UMI per Cell","2,391"
    "Median UMI per Cell","1,118"
    "Mean Genes per Cell","1,042"
    "Median Genes per Cell","716"
    "Total Genes Detected","16,068"

Inspect cell barcodes.

.. code-block:: console

    $ cat cell_barcodes.txt

    AACGGGCATGGG
    AACTGGGCATGG
    AAGACAGCGTGT
    AAGGGCGTACTC
    AATAAATACAAA
    AATAAATACAAC
    AATAAATACAAG
    AATAAATACAAT
    AATCAATCGCAA
    AATCAATCGCAC

Prepare feature barcodes. sgRNA sequences can be found in Supplementary Table 1.

.. code-block:: console

    $ cat feature_barcodes.tsv

    DNMT3B  CAGGATTGGGGGCGAGTCGG
    MBD1    ATAGGTGTCTGAGCGTCCAC
    TET2    CAGGACTCACACGACTATTC

|


Barcode extraction
==================

The transcripts derived from `CROPseq-Guide-Puro`_ and captured by Drop-Seq beads contain sgRNA sequences. There are no secondary libraries built on top of this single-cell RNA-Seq library for sgRNA enrichment. The transcripts derived from `CROPseq-Guide-Puro`_ are captured by the ployA tails. Therefore, the locations of sgRNA on read 2 vary. We need to extract the sgRNA sequences from read 2.

``qc`` mode is used for sgRNA extraction. Use ``-n`` to specify the number of reads to analyze, ``None`` means all the reads. Use ``-t`` to set the number of threads. By default, the diagnostic results and plots are generated in the ``qc`` directory (set by ``--output_directory``), and full length of read 1 and read 2 are searched against reference cell and feature barcodes, respectively. The per base content of both read pairs and the distribution of matched barcode positions are summarized. Use ``-r1_c`` and/or ``-r2_c`` to limit the search range for read 1 and read 2 respectively. Use ``-cb_n`` and/or ``-fb_n`` to set the mismatch tolerance for cell and feature barcode matching (default ``3``).

.. code-block:: console

    $ fba qc \
        -1 SRR5163029_1.fastq.gz \
        -2 SRR5163029_2.fastq.gz \
        -w cell_barcodes.txt \
        -f feature_barcodes.tsv \
        -cb_m 1 \
        -fb_m 1 \
        -cb_n 15 \
        -fb_n 15 \
        -r1_c 0,12 \
        -t $SLURM_CPUS_ON_NODE \
        --num_reads None

The first 12 bases are cell barcodes and the following 8 bases are UMIs. Based on the base content plot, the GC content of cell barcodes and UMIs are quite even.

.. image:: Pyplot_read1_per_base_seq_content.png
   :width: 1200px
   :align: center

|

.. image:: Pyplot_read1_barcodes_starting_ending.png
   :width: 1200px
   :align: center

As for read 2, based on the per base content, after the first 12 bases, the reads are A enriched.

.. image:: Pyplot_read2_per_base_seq_content.png
   :width: 1200px
   :align: center

|

.. image:: Pyplot_read2_barcodes_starting_ending.png
   :width: 1200px
   :align: center

The detailed ``qc`` results are stored in ``feature_barcoding_output.tsv.gz`` file. ``matching_pos`` columns indicate the matched positions on reads. ``matching_description`` columns indicate mismatches in substitutions:insertions:deletions format.

.. code-block:: console

    $ gzip -dc qc/feature_barcoding_output.tsv.gz | head

    read1_seq       cell_barcode    cb_matching_pos cb_matching_description read2_seq       feature_barcode fb_matching_pos fb_matching_description
    TTTAGGATCGTTtgatgtattttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttcttctttcttttttattctttacaacatcctaccataacata no_match        NA      NA      ATTAAAAATATTGTGGCAGGAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAACAAAAAAAAACAAAAAAAAATCAGCTATATAACCACTAATACTTCTA    NA      NA      NA
    GTCGAAACTCTTaacgggatttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttt no_match        NA      NA      TTATAATGGTTACAAATAAAGCAATAGCATCACAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA    NA      NA      NA
    GTTTACGTGTTCatgggcgattttttttttttttttttttttttttttttaaaaaagttaaaagggggcccgtggggggacaaatagaggggcctagagttccaccccccatcccacaaaaaaaaccctcaccgcacagggcctcgcccct GTTTACGTGTTC    0:12    0:0:0   GGAGTACGGAGAATTCTATAAGAGCTTGACCAATGACTGGGAAGATCACTTGGCAGTGAAGCATTTTTCAGTTGAAGGACAGTTGGAATTCAGAGCCCTTCTATTTGTCCCACGACGTGCTCCTTTTGATCTGTTTGAAAAAAAAAAAAAA    no_match        NA      NA
    CCGTCCTAGTTGgtgtatattttttttgtttttttttttttttcaccgggtcagagctgcccctaagtaccacgtcccgtcccacctttatcggacctcggccaccacaaattgcttatccagagtgcccccctccgcccatcccagactc CCGTCCTAGTTG    0:12    0:0:0   AATTAAGTCTCGTAAAGAACGAGAAGCTGAACTTGGACCTAGGGCAACCGACTTCACCAATGTTTACAGCGAGAATCTTGGTGACGACGTGGATGATGAGCGCCTTAAGGTTCTCTTTGGCAAGTTTGGGCCTGCCTTGAGTGTGCGACTT    no_match        NA      NA
    TTTCAGTATTGGggcgaattttttttttttttttttttttttttttttttttttttttttttttggctagtttttttgtggtttttgcttttggttctctcgtttgccctggagctcccaggtccctttcttgtcctaccataggtaaccc TTTCAGTATTGG    0:12    0:0:0   GGACGAAACACCGATAGGTGTCTGAGCGTCCACGTTTTAGAGCTAGAAATAGCAAGTTAAAATAAGGCTAGTCCGTTATCAACTTGAAAAAGTGGCACCGAGTCGGTGCTTTTTTAAGCTTGGCGTAACTAGATCTTGAGACACTGCTTTT MBD1_ATAGGTGTCTGAGCGTCCAC       13:33 0:0:0
    CTAGGTACCACTagacagtttttttttttttttttttttttttttttttttttttttttttctctatgtgtgcttttttttggctttagtctgtgggtccctagttagccccggcgcccccacgcgcagaacgtgtcttaccacaagaacc CTAGGTACCACT    0:12    0:0:0   TTCTTGGGTAGTTTGCAGTTTTTAAAATTATGTTTTAAAATGGACTATCATATGCTTACCGTAACTTGAAAGTATTTCGATTTCTTGGCTTTATATATCTTGTGGAACGGACGAAACACCGATAGGTGTCTGAGCGTCCACGTTTTAGAGC MBD1_ATAGGTGTCTGAGCGTCCAC       121:1410:0:0
    TCTTCCACTACCgtcccgtcttttttttttttttttttttttttttttttttttttttctttatgtcagttttttttgtgctttagtattgggttcccttgtttgcccgagggctcccaggcccagatttgggctaaccaaagggaccccg TCTTCCACTACC    0:12    0:0:0   ACCGATAGGTGTCTGAGCGTCCACGTTTTAGAGCTAGAAATAGCAAGTTAAAATAAGGCTAGTCCGTTATCAACTTGAAAAAGTGGCACCGAGTCGGTGCTTTTTTAAGCTTGGCGTAACTAGATCTTGAGACACTGCTTTTTGCTTGTAC MBD1_ATAGGTGTCTGAGCGTCCAC       4:24  0:0:0
    CTTAATTTGGTGggaagattttttttttttttttttttttttttttttaagtactttaagtaagctttttttaggctttagccgtgggttcccctgttagcccgggaggtccccgggcccaatctgggcctaacagagaggccccgtacaa CTTAATTTGGTG    0:12    0:0:0   CCGTAACTTGAAAGTATTTCGATTTCTTGGCTTTATATATCTTGTGGAAAGGACGAAACACCGCAGGACTCACACGACTCTTCGTTTTAGAGCTAGCAATAGCAAGTTAAAATAAGGCTAGTCCGTTATCAACTTGAAAAAGTGGCACCGT TET2_CAGGACTCACACGACTATTC       63:83 1:0:0
    TCGTACATACGGtggtttttttttttttttttttttttttttttttttttttttttttttttttttgtttttttttttttttgtttttttttttgtgtcctttgttttcactggggctcccaggtccatatccggtgttaccagagaaacc TCGTACATACGG    0:12    0:0:0   ATCATATGCTTACCGTAACTTGAAAGTATTTCGATTTCTTGGCTTTATATATCTTGTGGAAAGGACGAAACACCGCAGGATTGGGGGCGAGTCGGGTTTTAGAGCTAGAAATAGCAAGTTAAAATAAGGCTAGTCCGTTATCAACTTGAAA DNMT3B_CAGGATTGGGGGCGAGTCGG     75:95 0:0:0

11,597 of 227,621,653 read pairs have valid cell and feature barcodes. Note: This is not a sgRNA enrichment library.

|


Matrix generation
=================

Only fragments with correct (passed the criteria) cell and feature barcodes are included. UMI removal is powered by UMI-tools (`Smith, T., et al. 2017. Genome Res. 27, 491–499.`_). Use ``-us`` to set the UMI starting position on read 1 (default ``16``). Use ``-ul`` to set the UMI length (default ``12``). Fragments with UMI length less than this value are discarded. UMI deduplication method is set by ``-ud`` (default ``directional``). Use ``-um`` to set UMI deduplication mismatch threshold (default ``1``).

.. _`Smith, T., et al. 2017. Genome Res. 27, 491–499.`: http://www.genome.org/cgi/doi/10.1101/gr.209601.116

The generated feature count matrix can be easily imported into well-established single cell analysis packages: Seruat_ and Scanpy_.

.. _Seruat: https://satijalab.org/seurat/

.. _Scanpy: https://scanpy.readthedocs.io/en/stable/

.. code-block:: console

    $ fba count \
        -i feature_barcoding_output.tsv.gz \
        -o matrix_featurecount.csv.gz \
        -us 12 \
        -ul 8

Result summary.

11.76 % (1,364 out of 11,597) of read pairs with valid cell and feature barcodes are unique fragments.

.. code-block:: console

    2022-03-04 23:18:27,501 - fba.__main__ - INFO - fba version: 0.0.x
    2022-03-04 23:18:27,501 - fba.__main__ - INFO - Initiating logging ...
    2022-03-04 23:18:27,501 - fba.__main__ - INFO - Python version: 3.10
    2022-03-04 23:18:27,501 - fba.__main__ - INFO - Using count subcommand ...
    2022-03-04 23:18:31,494 - fba.count - INFO - UMI-tools version: 1.1.2
    2022-03-04 23:18:31,495 - fba.count - INFO - UMI start position on read 1 auto-detected, overriding -us
    2022-03-04 23:18:31,495 - fba.count - INFO - UMI length: 8
    2022-03-04 23:18:31,496 - fba.count - INFO - UMI-tools deduplication threshold: 1
    2022-03-04 23:18:31,496 - fba.count - INFO - UMI-tools deduplication method: directional
    2022-03-04 23:18:31,496 - fba.count - INFO - Header line: read1_seq cell_barcode cb_matching_pos cb_matching_description read2_seq feature_barcode fb_matching_pos fb_matching_description
    2022-03-04 23:18:31,581 - fba.count - INFO - Number of lines processed: 11,597
    2022-03-04 23:18:31,581 - fba.count - INFO - Number of cell barcodes detected: 523
    2022-03-04 23:18:31,582 - fba.count - INFO - Number of features detected: 3
    2022-03-04 23:18:31,608 - fba.count - INFO - Total UMIs after deduplication: 1,364
    2022-03-04 23:18:31,609 - fba.count - INFO - Median number of UMIs per cell: 1.0
    2022-03-04 23:18:31,615 - fba.__main__ - INFO - Done.

|


Demultiplexing
==============


Gaussian mixture model
----------------------

The implementation of demultiplexing method ``2`` (set by ``-dm``) is inspired by the method described on `10x Genomics’ website`_. Use ``-p`` to set the probability threshold for demulitplexing (default ``0.9``). Use ``-nc`` to set the number of positive cells for a feature to be included for demultiplexing (default ``200``).

.. _`10x Genomics’ website`: https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/algorithms/crispr

.. code-block:: console

    $ fba demultiplex \
        -i matrix_featurecount.csv.gz \
        -dm 2 \
        -v \
        -nc 0

.. code-block:: console

    2022-03-04 23:19:05,218 - fba.__main__ - INFO - fba version: 0.0.x
    2022-03-04 23:19:05,219 - fba.__main__ - INFO - Initiating logging ...
    2022-03-04 23:19:05,219 - fba.__main__ - INFO - Python version: 3.10
    2022-03-04 23:19:05,219 - fba.__main__ - INFO - Using demultiplex subcommand ...
    2022-03-04 23:19:15,199 - fba.__main__ - INFO - Skipping arguments: "-q/--quantile", "-cm/--clustering_method"
    2022-03-04 23:19:15,200 - fba.demultiplex - INFO - Output directory: demultiplexed
    2022-03-04 23:19:15,201 - fba.demultiplex - INFO - Demultiplexing method: 2
    2022-03-04 23:19:15,201 - fba.demultiplex - INFO - UMI normalization method: clr
    2022-03-04 23:19:15,201 - fba.demultiplex - INFO - Visualization: On
    2022-03-04 23:19:15,201 - fba.demultiplex - INFO - Visualization method: tsne
    2022-03-04 23:19:15,201 - fba.demultiplex - INFO - Loading feature count matrix: matrix_featurecount.csv.gz ...
    2022-03-04 23:19:15,219 - fba.demultiplex - INFO - Number of cells: 523
    2022-03-04 23:19:15,219 - fba.demultiplex - INFO - Number of positive cells for a feature to be included: 0
    2022-03-04 23:19:15,222 - fba.demultiplex - INFO - Number of features: 3 / 3 (after filtering / original in the matrix)
    2022-03-04 23:19:15,222 - fba.demultiplex - INFO - Features: DNMT3B MBD1 TET2
    2022-03-04 23:19:15,222 - fba.demultiplex - INFO - Total UMIs: 1,364 / 1,364
    2022-03-04 23:19:15,223 - fba.demultiplex - INFO - Median number of UMIs per cell: 1.0 / 1.0
    2022-03-04 23:19:15,223 - fba.demultiplex - INFO - Demultiplexing ...
    2022-03-04 23:19:17,319 - fba.demultiplex - INFO - Generating heatmap ...
    2022-03-04 23:19:17,784 - fba.demultiplex - INFO - Embedding ...
    2022-03-04 23:19:32,256 - fba.__main__ - INFO - Done.

Heatmap of the relative abundance of features (sgRNAs) across all cells. Each column represents a single cell.

.. image:: Pyplot_heatmap_cells_demultiplexed_gm.png
   :alt: Heatmap
   :width: 700px
   :align: center

Preview the demultiplexing result: the numbers of singlets and multiplets. Note: This is not a sgRNA enrichment library.

.. code-block:: python

    In [1]: import pandas as pd

    In [2]: m = pd.read_csv("demultiplexed/matrix_cell_identity.csv.gz", index_col=0)

    In [3]: m.loc[:, m.sum(axis=0) == 1].sum(axis=1)
    Out[3]:
    DNMT3B    141
    MBD1      150
    TET2      158
    dtype: int64

    In [4]: sum(m.sum(axis=0) > 1)
    Out[4]: 74
