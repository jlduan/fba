.. _tutorial_pseudo-bulk_SC3_v3_NextGem_DI_CellPlex_Jurkat_Raji_10K_Multiplex:


fba tutorial
============

Dataset: 10k 1:1 Mixture of Raji and Jurkat Cells Multiplexed, 2 CMOs

The detailed description of this dataset can be found `here`_.

.. _`here`: https://www.10xgenomics.com/resources/datasets/10-k-1-1-mixture-of-raji-and-jurkat-cells-multiplexed-2-cm-os-3-1-standard-6-0-0


Preparation
-----------

Fastq files and feature barcodes are prepared as described :ref:`here <tutorial_cellplex_SC3_v3_NextGem_DI_CellPlex_Jurkat_Raji_10K_Multiplex>`.


QC
--

Threshold: one mismatch
+++++++++++++++++++++++

In the ``qc`` subcommand, if ``-1`` (read 1) is omitted, bulk mode is enabled. The purpose of bulk mode is to help design and qc feature barcoding assays before the actual single cell experiments. For instance, you may want to estimate 1) how many reads have valid feature barcodes. This may reflect the specificity of the primers used for library construction and could suggest the number of reads needed for sequencing; 2) the distribution of feature barcodes. This reflects the biological aspect of the design.

Use ``-2`` to specify read 2, and ``-f`` to specify feature barcodes. Search range on reads 2 can be controlled by ``-r2_c``. In this example, only one mismatch is allowed for feature barcode matching (set by ``-fb_m``). Use ``-n`` to specify the number of reads to analyze (``None`` is to analyze all reads provided in the fastq file). By default, the distribution of feature barcodes detected is summarized in ``qc/feature_barcode_frequency.csv``.

.. code-block:: console

    $ fba qc \
        -2 ../SC3_v3_NextGem_DI_CellPlex_Jurkat_Raji_10K_1_multiplexing_capture_S1_combined_R2_001.fastq.gz \
        -f SC3_v3_NextGem_DI_CRISPR_10K_feature_ref.tsv \
        -r2_c 0,15 \
        -fb_m 1 \
        -n None

The content of ``qc/feature_barcode_frequency.csv``.

+------------------------+-----------+------------+
| feature barcode        | num_reads | percentage |
+------------------------+-----------+------------+
| CMO301_ATGAGGAATTCCTGC | 132435325 | 0.62351149 |
+------------------------+-----------+------------+
| CMO302_CATGCCAATAGAGCG | 79628216  | 0.37489323 |
+------------------------+-----------+------------+
| CMO308_CGGATTCCACATCAT | 320078    | 0.00150694 |
+------------------------+-----------+------------+
| CMO309_GTTGATCTATAACAG | 6445      | 3.03E-05   |
+------------------------+-----------+------------+
| CMO303_CCGTCGTCCAAGCAT | 4047      | 1.91E-05   |
+------------------------+-----------+------------+
| CMO304_AACGTTAATCACTCA | 2199      | 1.04E-05   |
+------------------------+-----------+------------+
| CMO307_AAGCTCGTTGGAAGA | 1735      | 8.17E-06   |
+------------------------+-----------+------------+
| CMO312_ACATGGTCAACGCTG | 1508      | 7.10E-06   |
+------------------------+-----------+------------+
| CMO306_AAGATGAGGTCTGTG | 1323      | 6.23E-06   |
+------------------------+-----------+------------+
| CMO310_GCAGGAGGTATCAAT | 532       | 2.50E-06   |
+------------------------+-----------+------------+
| CMO311_GAATCGTGATTCTTC | 502       | 2.36E-06   |
+------------------------+-----------+------------+
| CMO305_CGCGATATGGTCGGA | 472       | 2.22E-06   |
+------------------------+-----------+------------+


Result summary.

98.3% (212,402,382 / 216,070,514) of reads have valid feature barcodes. CMO301_ATGAGGAATTCCTGC and CMO302_CATGCCAATAGAGCG are the most abundant CMOs. They account for all most all of the valid reads. Although the valid read ratio is 1.663171 (132,435,325 / 79,628,216), cells labeled with them separately are mixed at 1: 1 ratio. See :ref:`here <tutorial_cellplex_SC3_v3_NextGem_DI_CellPlex_Jurkat_Raji_10K_Multiplex>` for more details.


.. code-block:: console

    2021-10-02 02:02:31,092 - fba.__main__ - INFO - fba version: 0.0.11
    2021-10-02 02:02:31,092 - fba.__main__ - INFO - Initiating logging ...
    2021-10-02 02:02:31,092 - fba.__main__ - INFO - Python version: 3.7
    2021-10-02 02:02:31,092 - fba.__main__ - INFO - Using qc subcommand ...
    2021-10-02 02:02:31,873 - fba.__main__ - INFO - Bulk mode enabled: only feature barcodes on reads 2 are analyzed
    2021-10-02 02:02:31,873 - fba.__main__ - INFO - Skipping arguments: "-w/--whitelist", "-cb_m/--cb_mismatches", "-r1_c/--read1_coordinate"
    2021-10-02 02:02:31,875 - fba.qc - INFO - Number of reference feature barcodes: 12
    2021-10-02 02:02:31,875 - fba.qc - INFO - Read 2 coordinates to search: [0, 15)
    2021-10-02 02:02:31,875 - fba.qc - INFO - Feature barcode maximum number of mismatches: 1
    2021-10-02 02:02:31,875 - fba.qc - INFO - Read 2 maximum number of N allowed: inf
    2021-10-02 02:02:31,875 - fba.qc - INFO - Number of read pairs to analyze: all
    2021-10-02 02:02:31,875 - fba.qc - INFO - Matching ...
    2021-10-02 02:04:19,871 - fba.qc - INFO - Reads processed: 10,000,000
    2021-10-02 02:06:06,844 - fba.qc - INFO - Reads processed: 20,000,000
    2021-10-02 02:07:53,987 - fba.qc - INFO - Reads processed: 30,000,000
    2021-10-02 02:09:40,854 - fba.qc - INFO - Reads processed: 40,000,000
    2021-10-02 02:11:27,502 - fba.qc - INFO - Reads processed: 50,000,000
    2021-10-02 02:13:14,277 - fba.qc - INFO - Reads processed: 60,000,000
    2021-10-02 02:15:02,641 - fba.qc - INFO - Reads processed: 70,000,000
    2021-10-02 02:16:51,149 - fba.qc - INFO - Reads processed: 80,000,000
    2021-10-02 02:18:40,463 - fba.qc - INFO - Reads processed: 90,000,000
    2021-10-02 02:20:30,099 - fba.qc - INFO - Reads processed: 100,000,000
    2021-10-02 02:22:19,651 - fba.qc - INFO - Reads processed: 110,000,000
    2021-10-02 02:24:09,364 - fba.qc - INFO - Reads processed: 120,000,000
    2021-10-02 02:25:59,016 - fba.qc - INFO - Reads processed: 130,000,000
    2021-10-02 02:27:48,634 - fba.qc - INFO - Reads processed: 140,000,000
    2021-10-02 02:29:38,323 - fba.qc - INFO - Reads processed: 150,000,000
    2021-10-02 02:31:28,018 - fba.qc - INFO - Reads processed: 160,000,000
    2021-10-02 02:33:17,585 - fba.qc - INFO - Reads processed: 170,000,000
    2021-10-02 02:35:07,168 - fba.qc - INFO - Reads processed: 180,000,000
    2021-10-02 02:36:56,770 - fba.qc - INFO - Reads processed: 190,000,000
    2021-10-02 02:38:46,487 - fba.qc - INFO - Reads processed: 200,000,000
    2021-10-02 02:40:36,129 - fba.qc - INFO - Reads processed: 210,000,000
    2021-10-02 02:41:42,628 - fba.qc - INFO - Number of reads processed: 216,070,514
    2021-10-02 02:41:42,628 - fba.qc - INFO - Number of reads w/ valid feature barcodes: 212,402,382
    2021-10-02 02:41:42,629 - fba.__main__ - INFO - Output file: qc/feature_barcode_frequency.csv
    2021-10-02 02:41:42,645 - fba.__main__ - INFO - Done.


|


Threshold: two mismatches
+++++++++++++++++++++++++

Let's relax the threshold to allow 2 mismatches for feature barcode matching (set by ``-fb_m``).

.. code-block:: console

    $ fba qc \
        -2 ../SC3_v3_NextGem_DI_CellPlex_Jurkat_Raji_10K_1_multiplexing_capture_S1_combined_R2_001.fastq.gz \
        -f SC3_v3_NextGem_DI_CRISPR_10K_feature_ref.tsv \
        -r2_c 0,15 \
        -fb_m 2 \
        -n None


The content of ``qc/feature_barcode_frequency.csv``.

+------------------------+-----------+-------------+
| feature barcode        | num_reads | percentage  |
+------------------------+-----------+-------------+
| CMO301_ATGAGGAATTCCTGC | 133957542 | 0.624153341 |
+------------------------+-----------+-------------+
| CMO302_CATGCCAATAGAGCG | 80322629  | 0.374250203 |
+------------------------+-----------+-------------+
| CMO308_CGGATTCCACATCAT | 323662    | 0.00150805  |
+------------------------+-----------+-------------+
| CMO309_GTTGATCTATAACAG | 6498      | 3.03E-05    |
+------------------------+-----------+-------------+
| CMO303_CCGTCGTCCAAGCAT | 4091      | 1.91E-05    |
+------------------------+-----------+-------------+
| CMO304_AACGTTAATCACTCA | 2225      | 1.04E-05    |
+------------------------+-----------+-------------+
| CMO307_AAGCTCGTTGGAAGA | 1751      | 8.16E-06    |
+------------------------+-----------+-------------+
| CMO312_ACATGGTCAACGCTG | 1535      | 7.15E-06    |
+------------------------+-----------+-------------+
| CMO306_AAGATGAGGTCTGTG | 1351      | 6.29E-06    |
+------------------------+-----------+-------------+
| CMO310_GCAGGAGGTATCAAT | 539       | 2.51E-06    |
+------------------------+-----------+-------------+
| CMO311_GAATCGTGATTCTTC | 507       | 2.36E-06    |
+------------------------+-----------+-------------+
| CMO305_CGCGATATGGTCGGA | 477       | 2.22E-06    |
+------------------------+-----------+-------------+


Result summary.

99.33% (214,622,807 / 216,070,514) of reads have valid feature barcodes.

.. code-block:: console

    2021-10-02 02:02:31,268 - fba.__main__ - INFO - fba version: 0.0.11
    2021-10-02 02:02:31,268 - fba.__main__ - INFO - Initiating logging ...
    2021-10-02 02:02:31,268 - fba.__main__ - INFO - Python version: 3.7
    2021-10-02 02:02:31,268 - fba.__main__ - INFO - Using qc subcommand ...
    2021-10-02 02:02:32,021 - fba.__main__ - INFO - Bulk mode enabled: only feature barcodes on reads 2 are analyzed
    2021-10-02 02:02:32,021 - fba.__main__ - INFO - Skipping arguments: "-w/--whitelist", "-cb_m/--cb_mismatches", "-r1_c/--read1_coordinate"
    2021-10-02 02:02:32,025 - fba.qc - INFO - Number of reference feature barcodes: 12
    2021-10-02 02:02:32,025 - fba.qc - INFO - Read 2 coordinates to search: [0, 15)
    2021-10-02 02:02:32,026 - fba.qc - INFO - Feature barcode maximum number of mismatches: 2
    2021-10-02 02:02:32,026 - fba.qc - INFO - Read 2 maximum number of N allowed: inf
    2021-10-02 02:02:32,026 - fba.qc - INFO - Number of read pairs to analyze: all
    2021-10-02 02:02:32,026 - fba.qc - INFO - Matching ...
    2021-10-02 02:13:36,407 - fba.qc - INFO - Reads processed: 10,000,000
    2021-10-02 02:24:40,718 - fba.qc - INFO - Reads processed: 20,000,000
    2021-10-02 02:35:43,572 - fba.qc - INFO - Reads processed: 30,000,000
    2021-10-02 02:46:45,598 - fba.qc - INFO - Reads processed: 40,000,000
    2021-10-02 02:57:47,743 - fba.qc - INFO - Reads processed: 50,000,000
    2021-10-02 03:08:49,904 - fba.qc - INFO - Reads processed: 60,000,000
    2021-10-02 03:19:52,124 - fba.qc - INFO - Reads processed: 70,000,000
    2021-10-02 03:30:54,289 - fba.qc - INFO - Reads processed: 80,000,000
    2021-10-02 03:41:56,459 - fba.qc - INFO - Reads processed: 90,000,000
    2021-10-02 03:53:01,896 - fba.qc - INFO - Reads processed: 100,000,000
    2021-10-02 04:04:07,940 - fba.qc - INFO - Reads processed: 110,000,000
    2021-10-02 04:15:13,882 - fba.qc - INFO - Reads processed: 120,000,000
    2021-10-02 04:26:19,716 - fba.qc - INFO - Reads processed: 130,000,000
    2021-10-02 04:37:25,780 - fba.qc - INFO - Reads processed: 140,000,000
    2021-10-02 04:48:31,630 - fba.qc - INFO - Reads processed: 150,000,000
    2021-10-02 04:59:36,756 - fba.qc - INFO - Reads processed: 160,000,000
    2021-10-02 05:10:42,247 - fba.qc - INFO - Reads processed: 170,000,000
    2021-10-02 05:21:47,635 - fba.qc - INFO - Reads processed: 180,000,000
    2021-10-02 05:32:53,151 - fba.qc - INFO - Reads processed: 190,000,000
    2021-10-02 05:43:58,739 - fba.qc - INFO - Reads processed: 200,000,000
    2021-10-02 05:55:04,397 - fba.qc - INFO - Reads processed: 210,000,000
    2021-10-02 06:01:48,423 - fba.qc - INFO - Number of reads processed: 216,070,514
    2021-10-02 06:01:48,424 - fba.qc - INFO - Number of reads w/ valid feature barcodes: 214,622,807
    2021-10-02 06:01:48,425 - fba.__main__ - INFO - Output file: qc/feature_barcode_frequency.csv
    2021-10-02 06:01:48,442 - fba.__main__ - INFO - Done.

|
