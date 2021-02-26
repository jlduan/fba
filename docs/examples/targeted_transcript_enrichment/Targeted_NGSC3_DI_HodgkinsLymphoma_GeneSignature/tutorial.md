
# fba tutorial

Dataset: Hodgkin's Lymphoma, Dissociated Tumor: Targeted, Gene Signature Panel

The detailed description of this dataset can be found at [here](https://support.10xgenomics.com/single-cell-gene-expression/datasets/4.0.0/Targeted_NGSC3_DI_HodgkinsLymphoma_GeneSignature).

<br>

## Preparation

Download fastq files.

```shell
$ wget https://cf.10xgenomics.com/samples/cell-exp/4.0.0/Targeted_NGSC3_DI_HodgkinsLymphoma_GeneSignature/Targeted_NGSC3_DI_HodgkinsLymphoma_GeneSignature_fastqs.tar

$ tar xvf Targeted_NGSC3_DI_HodgkinsLymphoma_GeneSignature_fastqs.tar
```

Download cell barcode info.

These are the cell-associated barcodes determined in the parent [Hodgkin's Lymphoma, Dissociated Tumor: Whole Transcriptome Analysis](https://support.10xgenomics.com/single-cell-gene-expression/datasets/4.0.0/Parent_NGSC3_DI_HodgkinsLymphoma) dataset.


```shell
$ wget https://cf.10xgenomics.com/samples/cell-exp/4.0.0/Parent_NGSC3_DI_HodgkinsLymphoma/Parent_NGSC3_DI_HodgkinsLymphoma_filtered_feature_bc_matrix.tar.gz

$ tar zxvf Parent_NGSC3_DI_HodgkinsLymphoma_filtered_feature_bc_matrix.tar.gz
```

Inspect cell barcodes.

```shell
$ gzip -dc filtered_feature_bc_matrix/barcodes.tsv.gz | head

AAACCCACAGGTTCGC-1
AAACCCACATCGATCA-1
AAACCCACATTGAGCT-1
AAACCCAGTATTTCCT-1
AAACGAACACGGTCTG-1
AAACGAAGTGGGTATG-1
AAACGAATCCTTATCA-1
AAACGCTAGTCTTCGA-1
AAACGCTGTTCTTGTT-1
AAAGAACAGATTCGAA-1
```

Prepare feature sequences.

`map` subcommand is designed to deal with secondary libraries built on top of the whole transcriptome assays. The transcripts of interest in the whole transcriptome libraries are enriched through hybridization or PCR amplification, which enables possibly more sensitive detection. The transcripts can be from endogenous genes or ectopic constructs, for example eGFP, as long as they are captured in the parent libraries.

In this tutorial, this 10x Genomics 'Targeted Gene Expression' library is used as an example. The [bait sequences](https://kb.10xgenomics.com/hc/en-us/articles/360045688071-What-are-the-bait-design-criteria-for-10x-pre-designed-and-custom-panels-) are used as references for read 2 mapping. Ideally, sequences of all possible captured transcribed regions should be included.

Read 1 contains cell barcodes and UMIs. Read 2 is expected captured transcribed regions. For feature references, they should only contain transcribed parts (non-overlapping and no introns for endogenous genes).

```shell
$ wget https://cf.10xgenomics.com/samples/cell-exp/4.0.0/Targeted_NGSC3_DI_HodgkinsLymphoma_GeneSignature/Targeted_NGSC3_DI_HodgkinsLymphoma_GeneSignature_target_panel.csv
```

Inspect feature reference info.

```shell
$ head Targeted_NGSC3_DI_HodgkinsLymphoma_GeneSignature_target_panel.csv

#panel_name=Human Gene Signature Panel
#panel_type=predesigned
#reference_genome=GRCh38
#reference_version=2020-A
#target_panel_file_format=1.0
gene_id,bait_seq,bait_id
ENSG00000000003,AGTTGTGGACGCTCGTAAGTTTTCGGCAGTTTCCGGGGAGACTCGGGGACTCCGCGTCTCGCTCTCTGTGTTCCAATCGCCCGGTGCGGTGGTGCAGGGTCTCGGGCTAGTCATGGCGTC,ENSG00000000003|TSPAN6|1
ENSG00000000003,CCCGTCTCGGAGACTGCAGACTAAACCAGTCATTACTTGTTTCAAGAGCGTTCTGCTAATCTACACTTTTATTTTCTGGATCACTGGCGTTATCCTTCTTGCAGTTGGCATTTGGGGCAA,ENSG00000000003|TSPAN6|2
ENSG00000000003,GGTGAGCCTGGAGAATTACTTTTCTCTTTTAAATGAGAAGGCCACCAATGTCCCCTTCGTGCTCATTGCTACTGGTACCGTCATTATTCTTTTGGGCACCTTTGGTTGTTTTGCTACCTG,ENSG00000000003|TSPAN6|3
ENSG00000000003,CCGAGCTTCTGCATGGATGCTAAAACTGTATGCAATGTTTCTGACTCTCGTTTTTTTGGTCGAACTGGTCGCTGCCATCGTAGGATTTGTTTTCAGACATGAGATTAAGAACAGCTTTAA,ENSG00000000003|TSPAN6|4
```

Re-format.

```shell
$ grep -v '#' Targeted_NGSC3_DI_HodgkinsLymphoma_GeneSignature_target_panel.csv | wc -l
53720

$ cut -d',' -f1,2 Targeted_NGSC3_DI_HodgkinsLymphoma_GeneSignature_target_panel.csv | gsed 's/,/\t/g' | grep -v '#' | head -53719 > Targeted_NGSC3_DI_HodgkinsLymphoma_GeneSignature_target_panel.tsv

$ head Targeted_NGSC3_DI_HodgkinsLymphoma_GeneSignature_target_panel.tsv

ENSG00000000003 AGTTGTGGACGCTCGTAAGTTTTCGGCAGTTTCCGGGGAGACTCGGGGACTCCGCGTCTCGCTCTCTGTGTTCCAATCGCCCGGTGCGGTGGTGCAGGGTCTCGGGCTAGTCATGGCGTC
ENSG00000000003 CCCGTCTCGGAGACTGCAGACTAAACCAGTCATTACTTGTTTCAAGAGCGTTCTGCTAATCTACACTTTTATTTTCTGGATCACTGGCGTTATCCTTCTTGCAGTTGGCATTTGGGGCAA
ENSG00000000003 GGTGAGCCTGGAGAATTACTTTTCTCTTTTAAATGAGAAGGCCACCAATGTCCCCTTCGTGCTCATTGCTACTGGTACCGTCATTATTCTTTTGGGCACCTTTGGTTGTTTTGCTACCTG
ENSG00000000003 CCGAGCTTCTGCATGGATGCTAAAACTGTATGCAATGTTTCTGACTCTCGTTTTTTTGGTCGAACTGGTCGCTGCCATCGTAGGATTTGTTTTCAGACATGAGATTAAGAACAGCTTTAA
ENSG00000000003 GAATAATTATGAGAAGGCTTTGAAGCAGTATAACTCTACAGGAGATTATAGAAGCCATGCAGTAGACAAGATCCAAAATACGTTGCATTGTTGTGGTGTCACCGATTATAGAGATTGGAC
ENSG00000000003 AGATACTAATTATTACTCAGAAAAAGGATTTCCTAAGAGTTGCTGTAAACTTGAAGATTGTACTCCACAGAGAGATGCAGACAAAGTAAACAATGAAGGTTGTTTTATAAAGGTGATGAC
ENSG00000000003 CATTATAGAGTCAGAAATGGGAGTCGTTGCAGGAATTTCCTTTGGAGTTGCTTGCTTCCAACTGATTGGAATCTTTCTCGCCTACTGCCTCTCTCGTGCCATAACAAATAACCAGTATGA
ENSG00000000003 GATAGTGTAACCCAATGTATCTGTGGGCCTATTCCTCTCTACCTTTAAGGACATTTAGGGTCCCCCCTGTGAATTAGAAAGTTGCTTGGCTGGAGAACTGACAACACTACTTACTGATAG
ENSG00000000003 ACCAAAAAACTACACCAGTAGGTTGATTCAATCAAGATGTATGTAGACCTAAAACTACACCAATAGGCTGATTCAATCAAGATCCGTGCTCGCAGTGGGCTGATTCAATCAAGATGTATG
ENSG00000000003 TTTGCTATGTTCTAAGTCCACCTTCTATCCCATTCATGTTAGATCGTTGAAACCCTGTATCCCTCTGAAACACTGGAAGAGCTAGTAAATTGTAAATGAAGTAATACTGTGTTCCTCTTG
```

<br>

## Matrix generation

First, all read 1 are searched against reference cell-associated barcodes. Use `-r1_coords` to set the search range, `-cb_m` to set the mismatching threshold. Read 2 with correct cell barcodes (on reads 1) is mapped to the provided sequences (bwa [Li, H. (2013)](https://arxiv.org/abs/1303.3997) or  bowtie2, [Langmead, B., et al. (2012)](http://dx.doi.org/10.1038/nmeth.1923)). The default aligner is bwa. Only alignments passed mapping quality threshold (set by `--mapq`) are kept for downstream feature counting. UMI deduplication is powered by UMI-tools ([Smith, T., et al. 2017. Genome Res. 27, 491–499.](http://www.genome.org/cgi/doi/10.1101/gr.209601.116)). Use `-us` to set the UMI starting position on read 1. Use `-ul` to set the UMI length. Fragments with UMI length less than this value are discarded. Use `-um` to set mismatch threshold. UMI deduplication method is set by `-ud`.

The generated feature count matrix can be easily imported into well-established single cell analysis packages: [Seruat](https://satijalab.org/seurat/) and [Scanpy](https://scanpy.readthedocs.io/en/stable/).

```shell
$ fba map \
    -1 Targeted_NGSC3_DI_HodgkinsLymphoma_GeneSignature_fastqs/Targeted_NGSC3_DI_HodgkinsLymphoma_GeneSignature_S1_L003_R1_001.fastq.gz \
    -2 Targeted_NGSC3_DI_HodgkinsLymphoma_GeneSignature_fastqs/Targeted_NGSC3_DI_HodgkinsLymphoma_GeneSignature_S1_L003_R2_001.fastq.gz \
    -w filtered_feature_bc_matrix/barcodes.tsv.gz \
    -f Targeted_NGSC3_DI_HodgkinsLymphoma_GeneSignature_target_panel.tsv \
    -o matrix_featurecount.csv.gz \
    -r1_coords 0,16 \
    -cb_m 1 \
    -al bwa \
    --mapq 10 \
    -us 16 \
    -ul 12 \
    -um 1 \
    -ud directional \
    --output_directory barcode_mapping
```

Result summary.

7.67% of total read pairs (2,405,998 of 31,372,024) contribute to the final expression matrix after UMI deduplication. Sequenced quite deep.

```shell
2021-02-17 23:33:59,615 - fba.__main__ - INFO - fba version: 0.0.7
2021-02-17 23:33:59,615 - fba.__main__ - INFO - Initiating logging ...
2021-02-17 23:33:59,615 - fba.__main__ - INFO - Python version: 3.7
2021-02-17 23:33:59,615 - fba.__main__ - INFO - Using map subcommand ...
2021-02-17 23:33:59,863 - fba.map - INFO - bwa version: 0.7.17
2021-02-17 23:34:02,116 - fba.map - INFO - samtools version: 1.3
2021-02-17 23:34:02,145 - fba.map - INFO - Number of reference cell barcodes: 3,394
2021-02-17 23:34:02,145 - fba.map - INFO - Read 1 coordinates to search: [0, 16)
2021-02-17 23:34:02,145 - fba.map - INFO - Cell barcode maximum number of mismatches: 1
2021-02-17 23:34:02,145 - fba.map - INFO - Read 1 maximum number of N allowed: 3
2021-02-17 23:34:02,145 - fba.map - INFO - Matching cell barcodes, read 1 ...
2021-02-17 23:47:07,994 - fba.map - INFO - number of read pairs processed: 31,372,024
2021-02-17 23:47:07,995 - fba.map - INFO - Number of read pairs w/ valid cell barcodes: 28,336,049
2021-02-17 23:47:08,024 - fba.map - INFO - Number of reference features: 1,142
2021-02-17 23:47:08,024 - fba.map - INFO - Number of threads: 56
2021-02-17 23:47:08,024 - fba.map - INFO - Aligning read 2 ...
2021-02-17 23:52:34,225 - fba.map - INFO -
[M::bwa_idx_load_from_disk] read 0 ALT contigs
[M::process] read 6222224 sequences (560000160 bp)...
[M::mem_process_seqs] Processed 6222224 reads in 1411.678 CPU sec, 87.647 real sec
[M::process] read 6222224 sequences (560000160 bp)...
[M::mem_process_seqs] Processed 6222224 reads in 484.034 CPU sec, 11.666 real sec
[M::process] read 6222224 sequences (560000160 bp)...
[M::mem_process_seqs] Processed 6222224 reads in 487.450 CPU sec, 11.070 real sec
[M::process] read 6222224 sequences (560000160 bp)...
[M::mem_process_seqs] Processed 6222224 reads in 457.438 CPU sec, 8.857 real sec
[M::process] read 3447153 sequences (310243770 bp)...
[M::mem_process_seqs] Processed 3447153 reads in 273.114 CPU sec, 8.418 real sec
[main] Version: 0.7.17-r1198-dirty
[main] CMD: /home2/s166631/bin/bwa mem -t 56 -C barcode_mapping/feature_ref.fasta barcode_mapping/modified.fq.gz
[main] Real time: 187.399 sec; CPU: 3189.278 sec
2021-02-17 23:52:34,226 - fba.map - INFO - Generating matrix (UMI deduplication) ...
2021-02-17 23:52:34,226 - fba.map - INFO - UMI-tools version: 1.1.1
2021-02-17 23:52:34,226 - fba.map - INFO - Mapping quality threshold: 10
2021-02-17 23:52:34,226 - fba.map - INFO - UMI starting position on read 1: 16
2021-02-17 23:52:34,226 - fba.map - INFO - UMI length: 12
2021-02-17 23:52:34,226 - fba.map - INFO - UMI-tools deduplication threshold: 1
2021-02-17 23:52:34,226 - fba.map - INFO - UMI-tools deduplication method: directional
2021-02-17 23:54:06,700 - fba.map - INFO - Number of cell barcodes detected: 3,379
2021-02-17 23:54:06,700 - fba.map - INFO - Number of features detected: 1,129
2021-02-17 23:54:06,704 - fba.map - INFO - Total UMIs after deduplication: 2,405,998
2021-02-17 23:54:06,713 - fba.map - INFO - Median number of UMIs per cell: 507.0
2021-02-17 23:54:11,085 - fba.__main__ - INFO - Done.
```
