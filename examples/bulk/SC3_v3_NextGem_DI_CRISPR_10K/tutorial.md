
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

| feature barcode                    | num_reads  | percentage         |
|------------------------------------|------------|--------------------|
| NON\_TARGET-1_AACGTGCTGACGATGCGGGC | 59,310,228 | 0.6979885931379053 |
| RAB1A-2_GCCGGCGAACCAGGAAATAG       | 25,662,834 | 0.3020114068620947 |

<br>

Result summary.

58.59% (84,973,062 / 145,032,428) of reads have valid feature barcodes.

```shell
2020-10-19 21:55:33,663 - fba.__main__ - INFO - fba version: 0.0.6
2020-10-19 21:55:33,663 - fba.__main__ - INFO - Initiating logging ...
2020-10-19 21:55:33,663 - fba.__main__ - INFO - Python version: 3.7
2020-10-19 21:55:33,663 - fba.__main__ - INFO - Using qc subcommand ...
2020-10-19 21:55:33,663 - fba.__main__ - INFO - Bulk mode enabled: only feature barcodes on reads 2 are analyzed
2020-10-19 21:55:33,663 - fba.__main__ - INFO - Skipping arguments: "-1", "-w", "-cb_m", "-r1_coords"
2020-10-19 21:55:33,664 - fba.qc - INFO - Number of reference feature barcodes: 2
2020-10-19 21:55:33,664 - fba.qc - INFO - Read 2 coordinates to search: [31, 51]
2020-10-19 21:55:33,664 - fba.qc - INFO - Feature barcode maximum number of mismatches: 1
2020-10-19 21:55:33,664 - fba.qc - INFO - Read 2 maximum number of N allowed: inf
2020-10-19 21:55:33,665 - fba.qc - INFO - Number of read pairs to analyze: all
2020-10-19 21:55:33,665 - fba.qc - INFO - Matching ...
2020-10-19 21:58:04,681 - fba.qc - INFO - Reads processed: 10,000,000
2020-10-19 22:00:35,884 - fba.qc - INFO - Reads processed: 20,000,000
2020-10-19 22:03:06,736 - fba.qc - INFO - Reads processed: 30,000,000
2020-10-19 22:05:36,955 - fba.qc - INFO - Reads processed: 40,000,000
2020-10-19 22:08:07,207 - fba.qc - INFO - Reads processed: 50,000,000
2020-10-19 22:10:37,627 - fba.qc - INFO - Reads processed: 60,000,000
2020-10-19 22:13:08,185 - fba.qc - INFO - Reads processed: 70,000,000
2020-10-19 22:15:38,854 - fba.qc - INFO - Reads processed: 80,000,000
2020-10-19 22:18:09,499 - fba.qc - INFO - Reads processed: 90,000,000
2020-10-19 22:20:40,123 - fba.qc - INFO - Reads processed: 100,000,000
2020-10-19 22:23:10,740 - fba.qc - INFO - Reads processed: 110,000,000
2020-10-19 22:25:41,470 - fba.qc - INFO - Reads processed: 120,000,000
2020-10-19 22:28:12,351 - fba.qc - INFO - Reads processed: 130,000,000
2020-10-19 22:30:43,402 - fba.qc - INFO - Reads processed: 140,000,000
2020-10-19 22:31:59,298 - fba.qc - INFO - Number of reads processed: 145,032,428
2020-10-19 22:31:59,299 - fba.qc - INFO - Number of reads w/ valid feature barcodes: 84,973,062
2020-10-19 22:31:59,299 - fba.__main__ - INFO - Output file: qc/feature_barcode_frequency.csv
2020-10-19 22:31:59,339 - fba.__main__ - INFO - Done.
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
2020-10-19 22:32:33,657 - fba.__main__ - INFO - fba version: 0.0.6
2020-10-19 22:32:33,658 - fba.__main__ - INFO - Initiating logging ...
2020-10-19 22:32:33,658 - fba.__main__ - INFO - Python version: 3.7
2020-10-19 22:32:33,658 - fba.__main__ - INFO - Using qc subcommand ...
2020-10-19 22:32:33,658 - fba.__main__ - INFO - Bulk mode enabled: only feature barcodes on reads 2 are analyzed
2020-10-19 22:32:33,658 - fba.__main__ - INFO - Skipping arguments: "-1", "-w", "-cb_m", "-r1_coords"
2020-10-19 22:32:33,670 - fba.qc - INFO - Number of reference feature barcodes: 2
2020-10-19 22:32:33,670 - fba.qc - INFO - Read 2 coordinates to search: [31, 51]
2020-10-19 22:32:33,670 - fba.qc - INFO - Feature barcode maximum number of mismatches: 2
2020-10-19 22:32:33,670 - fba.qc - INFO - Read 2 maximum number of N allowed: inf
2020-10-19 22:32:33,670 - fba.qc - INFO - Number of read pairs to analyze: all
2020-10-19 22:32:33,670 - fba.qc - INFO - Matching ...
2020-10-19 22:49:01,502 - fba.qc - INFO - Reads processed: 10,000,000
2020-10-19 23:05:29,402 - fba.qc - INFO - Reads processed: 20,000,000
2020-10-19 23:21:57,536 - fba.qc - INFO - Reads processed: 30,000,000
2020-10-19 23:38:28,123 - fba.qc - INFO - Reads processed: 40,000,000
2020-10-19 23:54:55,946 - fba.qc - INFO - Reads processed: 50,000,000
2020-10-20 00:11:23,748 - fba.qc - INFO - Reads processed: 60,000,000
2020-10-20 00:27:50,410 - fba.qc - INFO - Reads processed: 70,000,000
2020-10-20 00:44:17,468 - fba.qc - INFO - Reads processed: 80,000,000
2020-10-20 01:00:45,392 - fba.qc - INFO - Reads processed: 90,000,000
2020-10-20 01:17:12,952 - fba.qc - INFO - Reads processed: 100,000,000
2020-10-20 01:33:40,369 - fba.qc - INFO - Reads processed: 110,000,000
2020-10-20 01:50:07,896 - fba.qc - INFO - Reads processed: 120,000,000
2020-10-20 02:06:33,881 - fba.qc - INFO - Reads processed: 130,000,000
2020-10-20 02:23:00,000 - fba.qc - INFO - Reads processed: 140,000,000
2020-10-20 02:31:16,393 - fba.qc - INFO - Number of reads processed: 145,032,428
2020-10-20 02:31:16,394 - fba.qc - INFO - Number of reads w/ valid feature barcodes: 100,307,853
2020-10-20 02:31:16,395 - fba.__main__ - INFO - Output file: qc/feature_barcode_frequency.csv
2020-10-20 02:31:16,427 - fba.__main__ - INFO - Done.
```
