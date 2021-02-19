
# fba tutorial

Dataset: 10k A375 Cells Transduced with (1) Non-Target and (1) Target sgRNA, Dual Indexed

The detailed description of this dataset can be found at [here](https://support.10xgenomics.com/single-cell-gene-expression/datasets/4.0.0/SC3_v3_NextGem_DI_CRISPR_10K).

<br>

## Preparation

Fastq files and feature barcodes are prepared as described at [here](https://github.com/jlduan/fba/blob/master/examples/crispr_screening/SC3_v3_NextGem_DI_CRISPR_10K/tutorial.md).

<br>

## QC

In the qc subcommand, if `-1` (read 1) is omitted, bulk mode is enabled. The purpose of bulk mode is to help design and qc feature barcoding assays before the actual single cell experiments. For instance, you may want to estimate 1) how many reads have valid feature barcodes. This may reflect the specificity of the primers used for library construction and could suggest the number of reads needed for sequencing; 2) the distribution of feature barcodes. This reflects the biological aspect of the design.

Use `-2` to specify read 2, and `-f` to specify feature barcodes. Search range on reads 2 can be controlled by `-r2_coords`. In this example, only one mismatch is allowed for feature barcode matching (set by `-fb_m`). Use `-n` to specify the number of reads to analyze (`None` is to analyze all reads provided in the fastq file). By default, the distribution of feature barcodes detected is summarized in `qc/feature_barcode_frequency.csv`.

```shell
$ fba qc \
    -2 SC3_v3_NextGem_DI_CRISPR_10K_crispr_S1_combined_R2_001.fastq.gz \
    -f SC3_v3_NextGem_DI_CRISPR_10K_feature_ref_edited.tsv \
    -r2_coords 31,51 \
    -fb_m 1 \
    -n None
```

<br>

The content of `qc/feature_barcode_frequency.csv`.

| feature barcode                    | num_reads  | percentage        |
|------------------------------------|------------|-------------------|
| NON\_TARGET-1_AACGTGCTGACGATGCGGGC | 59,310,228 | 0.697988593137905 |
| RAB1A-2_GCCGGCGAACCAGGAAATAG       | 25,662,834 | 0.302011406862095 |

<br>

Result summary.

58.59% (84,973,062 / 145,032,428) of reads have valid feature barcodes.

```shell
2021-02-17 16:12:35,393 - fba.__main__ - INFO - fba version: 0.0.7
2021-02-17 16:12:35,393 - fba.__main__ - INFO - Initiating logging ...
2021-02-17 16:12:35,393 - fba.__main__ - INFO - Python version: 3.7
2021-02-17 16:12:35,393 - fba.__main__ - INFO - Using qc subcommand ...
2021-02-17 16:12:35,394 - fba.__main__ - INFO - Bulk mode enabled: only feature barcodes on reads 2 are analyzed
2021-02-17 16:12:35,394 - fba.__main__ - INFO - Skipping arguments: "-1", "-w", "-cb_m", "-r1_coords"
2021-02-17 16:12:35,395 - fba.qc - INFO - Number of reference feature barcodes: 2
2021-02-17 16:12:35,395 - fba.qc - INFO - Read 2 coordinates to search: [31, 51)
2021-02-17 16:12:35,395 - fba.qc - INFO - Feature barcode maximum number of mismatches: 1
2021-02-17 16:12:35,395 - fba.qc - INFO - Read 2 maximum number of N allowed: inf
2021-02-17 16:12:35,395 - fba.qc - INFO - Number of read pairs to analyze: all
2021-02-17 16:12:35,395 - fba.qc - INFO - Matching ...
2021-02-17 16:15:07,684 - fba.qc - INFO - Reads processed: 10,000,000
2021-02-17 16:17:39,083 - fba.qc - INFO - Reads processed: 20,000,000
2021-02-17 16:20:09,116 - fba.qc - INFO - Reads processed: 30,000,000
2021-02-17 16:22:38,981 - fba.qc - INFO - Reads processed: 40,000,000
2021-02-17 16:25:11,671 - fba.qc - INFO - Reads processed: 50,000,000
2021-02-17 16:27:44,790 - fba.qc - INFO - Reads processed: 60,000,000
2021-02-17 16:30:18,110 - fba.qc - INFO - Reads processed: 70,000,000
2021-02-17 16:32:51,391 - fba.qc - INFO - Reads processed: 80,000,000
2021-02-17 16:35:24,625 - fba.qc - INFO - Reads processed: 90,000,000
2021-02-17 16:37:57,678 - fba.qc - INFO - Reads processed: 100,000,000
2021-02-17 16:40:30,706 - fba.qc - INFO - Reads processed: 110,000,000
2021-02-17 16:43:03,867 - fba.qc - INFO - Reads processed: 120,000,000
2021-02-17 16:45:37,197 - fba.qc - INFO - Reads processed: 130,000,000
2021-02-17 16:48:10,511 - fba.qc - INFO - Reads processed: 140,000,000
2021-02-17 16:49:27,662 - fba.qc - INFO - Number of reads processed: 145,032,428
2021-02-17 16:49:27,663 - fba.qc - INFO - Number of reads w/ valid feature barcodes: 84,973,062
2021-02-17 16:49:27,664 - fba.__main__ - INFO - Output file: qc/feature_barcode_frequency.csv
2021-02-17 16:49:27,689 - fba.__main__ - INFO - Done.
```

<br>

<br>


Let's relax the threshold to allow 2 mismatches for feature barcode matching (set by `-fb_m`).

```shell
$ fba qc \
    -2 SC3_v3_NextGem_DI_CRISPR_10K_crispr_S1_combined_R2_001.fastq.gz \
    -f SC3_v3_NextGem_DI_CRISPR_10K_feature_ref_edited.tsv \
    -r2_coords 31,51 \
    -fb_m 2 \
    -n None
```

<br>

The content of `qc/feature_barcode_frequency.csv`.

| feature barcode                    | num_reads  | percentage         |
|------------------------------------|------------|--------------------|
| NON\_TARGET-1_AACGTGCTGACGATGCGGGC | 66,334,740 | 0.6613115326075217 |
| RAB1A-2_GCCGGCGAACCAGGAAATAG       | 33,973,113 | 0.3386884673924782 |

<br>

Result summary.

69.16% (100,307,853 / 145,032,428) of reads have valid feature barcodes.

```shell
2021-02-17 16:12:00,407 - fba.__main__ - INFO - fba version: 0.0.7
2021-02-17 16:12:00,407 - fba.__main__ - INFO - Initiating logging ...
2021-02-17 16:12:00,408 - fba.__main__ - INFO - Python version: 3.7
2021-02-17 16:12:00,408 - fba.__main__ - INFO - Using qc subcommand ...
2021-02-17 16:12:00,408 - fba.__main__ - INFO - Bulk mode enabled: only feature barcodes on reads 2 are analyzed
2021-02-17 16:12:00,408 - fba.__main__ - INFO - Skipping arguments: "-1", "-w", "-cb_m", "-r1_coords"
2021-02-17 16:12:00,426 - fba.qc - INFO - Number of reference feature barcodes: 2
2021-02-17 16:12:00,426 - fba.qc - INFO - Read 2 coordinates to search: [31, 51)
2021-02-17 16:12:00,426 - fba.qc - INFO - Feature barcode maximum number of mismatches: 2
2021-02-17 16:12:00,426 - fba.qc - INFO - Read 2 maximum number of N allowed: inf
2021-02-17 16:12:00,426 - fba.qc - INFO - Number of read pairs to analyze: all
2021-02-17 16:12:00,426 - fba.qc - INFO - Matching ...
2021-02-17 16:28:02,710 - fba.qc - INFO - Reads processed: 10,000,000
2021-02-17 16:44:07,554 - fba.qc - INFO - Reads processed: 20,000,000
2021-02-17 17:00:13,431 - fba.qc - INFO - Reads processed: 30,000,000
2021-02-17 17:16:17,034 - fba.qc - INFO - Reads processed: 40,000,000
2021-02-17 17:32:21,635 - fba.qc - INFO - Reads processed: 50,000,000
2021-02-17 17:48:26,948 - fba.qc - INFO - Reads processed: 60,000,000
2021-02-17 18:04:31,050 - fba.qc - INFO - Reads processed: 70,000,000
2021-02-17 18:20:34,413 - fba.qc - INFO - Reads processed: 80,000,000
2021-02-17 18:36:38,778 - fba.qc - INFO - Reads processed: 90,000,000
2021-02-17 18:52:44,033 - fba.qc - INFO - Reads processed: 100,000,000
2021-02-17 19:08:49,500 - fba.qc - INFO - Reads processed: 110,000,000
2021-02-17 19:24:56,356 - fba.qc - INFO - Reads processed: 120,000,000
2021-02-17 19:41:02,072 - fba.qc - INFO - Reads processed: 130,000,000
2021-02-17 19:57:09,967 - fba.qc - INFO - Reads processed: 140,000,000
2021-02-17 20:05:15,665 - fba.qc - INFO - Number of reads processed: 145,032,428
2021-02-17 20:05:15,666 - fba.qc - INFO - Number of reads w/ valid feature barcodes: 100,307,853
2021-02-17 20:05:15,667 - fba.__main__ - INFO - Output file: qc/feature_barcode_frequency.csv
2021-02-17 20:05:15,701 - fba.__main__ - INFO - Done.
```
