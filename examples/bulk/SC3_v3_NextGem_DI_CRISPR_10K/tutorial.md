
# fba tutorial

Dataset: 10k A375 Cells Transduced with (1) Non-Target and (1) Target sgRNA, Dual Indexed

The detailed description of this dataset can be found at [here](https://support.10xgenomics.com/single-cell-gene-expression/datasets/4.0.0/SC3_v3_NextGem_DI_CRISPR_10K).

<br>

## Preparation

Cell and feature barcodes are prepared as described at [here](https://github.com/jlduan/fba/blob/master/examples/crispr_screening/SC3_v3_NextGem_DI_CRISPR_10K/tutorial.md).

<br>

## QC

In the qc subcommand, if `-1` (read 1) is omitted, bulk mode is enabled. The purpose of bulk mode is to help design and qc feature barcoding assays before the actual single cell sequencing. For instance, you may want to estimate 1) how many reads have valid feature barcodes. This may reflect the specificity of the primers used for library construction and suggest the number of reads needed for sequencing ; 2) the distribution of feature barcodes. This reflects the biological aspect of the design.

Use `-2` to specify reads 2, and `-f` to specify feature barcodes. Search range on reads 2 can be controlled by `-r2_coords`. In this example, only one mismatch is allowed for feature barcode matching (set by `-fb_m`). Use `-n` to specify the number of reads to analyze (`None` is to analyze all reads provided in the fastq file). By default, the distribution of feature barcodes detected is summarized in `qc/feature_barcode_frequency.csv`.

```shell
fba qc \
    -2 SC3_v3_NextGem_DI_CRISPR_10K_crispr_S1_combined_R2_001.fastq.gz \
    -f SC3_v3_NextGem_DI_CRISPR_10K_feature_ref_edited.tsv \
    -r2_coords 31,51 \
    -fb_m 1 \
    -n None
```

The content of `qc/feature_barcode_frequency.csv`.

| feature barcode                    | num_reads  | percentage         |
|------------------------------------|------------|--------------------|
| NON\_TARGET-1_AACGTGCTGACGATGCGGGC | 59,310,228 | 0.6979885931379053 |
| RAB1A-2_GCCGGCGAACCAGGAAATAG       | 25,662,834 | 0.3020114068620947 |


Result summary.

58.59% (84,973,062 / 145,032,428) of reads have valid feature barcodes.

```shell
2020-10-05 20:03:29,462 - fba.__main__ - INFO - fba version: 0.0.5dev
2020-10-05 20:03:29,462 - fba.__main__ - INFO - Initiating logging ...
2020-10-05 20:03:29,462 - fba.__main__ - INFO - Python version: 3.7
2020-10-05 20:03:29,462 - fba.__main__ - INFO - Using qc subcommand ...
2020-10-05 20:03:29,462 - fba.__main__ - INFO - Bulk mode enabled: only feature barcodes on reads 2 are analyzed
2020-10-05 20:03:29,462 - fba.__main__ - INFO - Skipping arguments: "-1", "-w", "-cb_m", "-r1_coords"
2020-10-05 20:03:29,463 - fba.qc - INFO - Number of reference feature barcodes: 2
2020-10-05 20:03:29,463 - fba.qc - INFO - Read 2 coordinates to search: [31, 51]
2020-10-05 20:03:29,463 - fba.qc - INFO - Feature barcode maximum number of mismatches: 1
2020-10-05 20:03:29,463 - fba.qc - INFO - Read 2 maximum number of N allowed: inf
2020-10-05 20:03:29,463 - fba.qc - INFO - Number of read pairs to analyze: all
2020-10-05 20:03:29,463 - fba.qc - INFO - Matching ...
2020-10-05 20:07:46,653 - fba.qc - INFO - Reads processed: 10,000,000
2020-10-05 20:12:03,447 - fba.qc - INFO - Reads processed: 20,000,000
2020-10-05 20:16:20,005 - fba.qc - INFO - Reads processed: 30,000,000
2020-10-05 20:20:36,108 - fba.qc - INFO - Reads processed: 40,000,000
2020-10-05 20:24:52,216 - fba.qc - INFO - Reads processed: 50,000,000
2020-10-05 20:29:08,310 - fba.qc - INFO - Reads processed: 60,000,000
2020-10-05 20:33:24,534 - fba.qc - INFO - Reads processed: 70,000,000
2020-10-05 20:37:40,703 - fba.qc - INFO - Reads processed: 80,000,000
2020-10-05 20:41:56,879 - fba.qc - INFO - Reads processed: 90,000,000
2020-10-05 20:46:13,009 - fba.qc - INFO - Reads processed: 100,000,000
2020-10-05 20:50:29,208 - fba.qc - INFO - Reads processed: 110,000,000
2020-10-05 20:54:45,425 - fba.qc - INFO - Reads processed: 120,000,000
2020-10-05 20:59:01,479 - fba.qc - INFO - Reads processed: 130,000,000
2020-10-05 21:03:17,671 - fba.qc - INFO - Reads processed: 140,000,000
2020-10-05 21:05:26,626 - fba.qc - INFO - Number of reads processed: 145,032,428
2020-10-05 21:05:26,626 - fba.qc - INFO - Number of reads w/ valid feature barcodes: 84,973,062
2020-10-05 21:05:26,627 - fba.__main__ - INFO - Output file: qc/feature_barcode_frequency.csv
2020-10-05 21:05:26,652 - fba.__main__ - INFO - Done.
```

<br>

Let's relax the threshold to allow 2 mismatches for feature barcode matching (set by `-fb_m`).

```shell
fba qc \
    -2 SC3_v3_NextGem_DI_CRISPR_10K_crispr_S1_combined_R2_001.fastq.gz \
    -f SC3_v3_NextGem_DI_CRISPR_10K_feature_ref_edited.tsv \
    -r2_coords 31,51 \
    -fb_m 2 \
    -n None
```

The content of `qc/feature_barcode_frequency.csv`.

| feature barcode                    | num_reads  | percentage         |
|------------------------------------|------------|--------------------|
| NON\_TARGET-1_AACGTGCTGACGATGCGGGC | 66,334,740 | 0.6613115326075217 |
| RAB1A-2_GCCGGCGAACCAGGAAATAG       | 33,973,113 | 0.3386884673924782 |





Result summary.

69.16% (100,307,853 / 145,032,428) of reads have valid feature barcodes.

```shell
2020-10-05 19:45:25,903 - fba.__main__ - INFO - fba version: 0.0.5dev
2020-10-05 19:45:25,903 - fba.__main__ - INFO - Initiating logging ...
2020-10-05 19:45:25,903 - fba.__main__ - INFO - Python version: 3.7
2020-10-05 19:45:25,903 - fba.__main__ - INFO - Using qc subcommand ...
2020-10-05 19:45:25,903 - fba.__main__ - INFO - Bulk mode enabled: only feature barcodes on reads 2 are analyzed
2020-10-05 19:45:25,903 - fba.__main__ - INFO - Skipping arguments: "-1", "-w", "-cb_m", "-r1_coords"
2020-10-05 19:45:25,904 - fba.qc - INFO - Number of reference feature barcodes: 2
2020-10-05 19:45:25,904 - fba.qc - INFO - Read 2 coordinates to search: [31, 51]
2020-10-05 19:45:25,904 - fba.qc - INFO - Feature barcode maximum number of mismatches: 2
2020-10-05 19:45:25,904 - fba.qc - INFO - Read 2 maximum number of N allowed: inf
2020-10-05 19:45:25,904 - fba.qc - INFO - Number of read pairs to analyze: all
2020-10-05 19:45:25,904 - fba.qc - INFO - Matching ...
2020-10-05 20:15:33,885 - fba.qc - INFO - Reads processed: 10,000,000
2020-10-05 20:45:31,636 - fba.qc - INFO - Reads processed: 20,000,000
2020-10-05 21:15:22,443 - fba.qc - INFO - Reads processed: 30,000,000
2020-10-05 21:45:07,698 - fba.qc - INFO - Reads processed: 40,000,000
2020-10-05 22:14:43,093 - fba.qc - INFO - Reads processed: 50,000,000
2020-10-05 22:44:24,259 - fba.qc - INFO - Reads processed: 60,000,000
2020-10-05 23:14:10,449 - fba.qc - INFO - Reads processed: 70,000,000
2020-10-05 23:43:56,872 - fba.qc - INFO - Reads processed: 80,000,000
2020-10-06 00:13:42,274 - fba.qc - INFO - Reads processed: 90,000,000
2020-10-06 00:43:17,054 - fba.qc - INFO - Reads processed: 100,000,000
2020-10-06 01:12:52,270 - fba.qc - INFO - Reads processed: 110,000,000
2020-10-06 01:42:27,214 - fba.qc - INFO - Reads processed: 120,000,000
2020-10-06 02:11:59,673 - fba.qc - INFO - Reads processed: 130,000,000
2020-10-06 02:41:34,684 - fba.qc - INFO - Reads processed: 140,000,000
2020-10-06 02:56:27,953 - fba.qc - INFO - Number of reads processed: 145,032,428
2020-10-06 02:56:27,953 - fba.qc - INFO - Number of reads w/ valid feature barcodes: 100,307,853
2020-10-06 02:56:27,954 - fba.__main__ - INFO - Output file: qc/feature_barcode_frequency.csv
2020-10-06 02:56:27,980 - fba.__main__ - INFO - Done.
```
