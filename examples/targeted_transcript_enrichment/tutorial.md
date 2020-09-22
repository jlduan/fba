
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

Download cell barcode info. These are the cell-associated barcodes in the parent [Hodgkin's Lymphoma, Dissociated Tumor: Whole Transcriptome Analysis](https://support.10xgenomics.com/single-cell-gene-expression/datasets/4.0.0/Parent_NGSC3_DI_HodgkinsLymphoma) dataset.


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

Prepare feature barcodes. `map` subcommand is designed to deal with secondary libraries built on top of whole transcriptome assays. The transcripts of interest in the whole transcriptome libraries are enriched through hybridization or PCR amplification, which enables possibly more sensitive detection. The transcripts can be from endogenous genes or ectopic constructs, for example eGFP, as long as they are captured in the parent libraries.

In this tutorial, this 10x Genomics 'Targeted Gene Expression' library is used as example. The bait sequences are used as references. Ideally, sequences of all possible captured transcript regions should be included. Reads 1 contain cell barcodes and UMIs. Reads 2 are enriched transcripts. For feature references, they should only contain transcribed parts (non-overlapping and no introns for endogenous genes).

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

Clean up.

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

First, all read 1 are searched against known cell-associated barcodes. Use `-r1_coords` to set the search range, `-cb_m` to set the mismatching threshold. Reads 2 with correct cell barcodes (on reads 1) are mapped to the provided sequences (bowtie2; [Langmead, B., and Salzberg, S.L. 2012. Nat. Methods 9, 357–359.](http://dx.doi.org/10.1038/nmeth.1923)). Only alignments passed mapping quality threshold (set by `--mapq`) are kept for downstream processing. UMI removal is powered by UMI-tools ([Smith, T., et al. 2017. Genome Res. 27, 491–499.](http://www.genome.org/cgi/doi/10.1101/gr.209601.116)). Use `-us` to set the UMI starting position on read 1. Use `-ul` to set the UMI length. Fragments with UMI length less than this value are discarded. UMI deduplication method is set by `-ud`.

The generated feature count matrix can be easily imported into well-established single cell analysis packages: [Seruat](https://satijalab.org/seurat/) and [Scanpy](https://scanpy.readthedocs.io/en/stable/).


```shell
$ fba map \
    -1 Targeted_NGSC3_DI_HodgkinsLymphoma_GeneSignature_fastqs/Targeted_NGSC3_DI_HodgkinsLymphoma_GeneSignature_S1_L003_R1_001.fastq.gz \
    -2 Targeted_NGSC3_DI_HodgkinsLymphoma_GeneSignature_fastqs/Targeted_NGSC3_DI_HodgkinsLymphoma_GeneSignature_S1_L003_R2_001.fastq.gz \
    -w filtered_feature_bc_matrix/barcodes.tsv.gz \
    -f Targeted_NGSC3_DI_HodgkinsLymphoma_GeneSignature_target_panel.tsv \
    -o matrix_featurecount.csv.gz \
    -r1_coords 0,16 \
    -cb_m 2 \
    --mapq 10 \
    -ul 12 \
    -us 16 \
    -um 1 \
    -ud directional \
    --output_directory barcode_mapping \
    -t $SLURM_CPUS_ON_NODE
```

Result summary.

```shell
2020-09-21 20:44:58,009 - fba.map - INFO - bowtie2 version: 2.4.1
2020-09-21 20:44:58,052 - fba.map - INFO - samtools version: 1.3
2020-09-21 20:44:58,112 - fba.map - INFO - Read 1 coordinates to search: [0, 16]
2020-09-21 20:44:58,112 - fba.map - INFO - Cell barcode maximum number of mismatches: 2
2020-09-21 20:44:58,112 - fba.map - INFO - Mapping quality threshold: 10
2020-09-21 20:44:58,113 - fba.map - INFO - Number of threads: 56
2020-09-21 20:44:58,113 - fba.map - INFO - Chunk size: 10,000
2020-09-21 20:44:58,113 - fba.map - INFO - UMI-tools version: 1.0.1
2020-09-21 20:44:58,113 - fba.map - INFO - UMI-tools deduplication method: directional
2020-09-21 20:44:58,113 - fba.map - INFO - UMI-tools deduplication threshold: 1
2020-09-21 20:44:58,113 - fba.map - INFO - UMI length: 12
2020-09-21 20:44:58,113 - fba.map - INFO - UMI starting position on read 1: 16
2020-09-21 20:44:58,113 - fba.map - INFO - Preparing feature barcodes ...
2020-09-21 20:45:00,738 - fba.map - INFO - Matching cell barcodes ...
2020-09-21 21:19:54,318 - fba.map - INFO - Aligning ...
2020-09-21 21:42:33,714 - fba.map - INFO - Generating matrix (UMI deduplication) ...
2020-09-21 21:43:33,317 - fba.map - INFO - Number of cell barcodes detected: 3,379
2020-09-21 21:43:33,319 - fba.map - INFO - Number of features detected: 1,127
2020-09-21 21:43:33,322 - fba.map - INFO - Total UMIs after deduplication: 1,906,780
2020-09-21 21:43:33,332 - fba.map - INFO - Median number of UMIs per cell: 397.0
```
