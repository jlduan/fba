# Changelog

## **0.0.13 (Jan 1 2023)**

Added

-   Multiple `-i/--input` flags support in `count` module
-   Better single-cell ATAC-seq support
    -   Setting `-ul/--umi_length` to 0 to disable UMI deduplication in `count` module
    -   Option `-cb_rc/--cell_barcode_reverse_complement` to convert the input cell barcode sequences to reverse-complement in `regex`, `extract` and `count` modules

## **0.0.12 (Mar 9 2022)**

Added

-   Non-uniform read lengths support in `qc` module

Fixed

-   Readability of the `qc` results
-   Matplotlib `FixedFormatter` warning

## **0.0.11 (Jun 17 2021)**

Changed

-   Rename options `-r1_coords` and `-r2_coords` to `-r1_c` and `-r2_c`, respectively.
-   Improve logging information.
