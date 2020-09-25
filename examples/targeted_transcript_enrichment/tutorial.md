
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

These are the cell-associated barcodes in the parent [Hodgkin's Lymphoma, Dissociated Tumor: Whole Transcriptome Analysis](https://support.10xgenomics.com/single-cell-gene-expression/datasets/4.0.0/Parent_NGSC3_DI_HodgkinsLymphoma) dataset.


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

Prepare feature barcodes.

`map` subcommand is designed to deal with secondary libraries built on top of whole transcriptome assays. The transcripts of interest in the whole transcriptome libraries are enriched through hybridization or PCR amplification, which enables possibly more sensitive detection. The transcripts can be from endogenous genes or ectopic constructs, for example eGFP, as long as they are captured in the parent libraries.

In this tutorial, this 10x Genomics 'Targeted Gene Expression' library is used as an example. The [bait sequences](https://kb.10xgenomics.com/hc/en-us/articles/360045688071-What-are-the-bait-design-criteria-for-10x-pre-designed-and-custom-panels-) are used as references for read 2 mapping. Ideally, sequences of all possible captured transcribed regions should be included. 

Reads 1 contain cell barcodes and UMIs. Reads 2 are expected captured transcribed regions. For feature references, they should only contain transcribed parts (non-overlapping and no introns for endogenous genes).

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

First, all reads 1 are searched against known cell-associated barcodes. Use `-r1_coords` to set the search range, `-cb_m` to set the mismatching threshold. Reads 2 with correct cell barcodes (on reads 1) are mapped to the provided sequences (bowtie2; [Langmead, B., and Salzberg, S.L. 2012. Nat. Methods 9, 357–359.](http://dx.doi.org/10.1038/nmeth.1923)). Only alignments passed mapping quality threshold (set by `--mapq`) are kept for downstream feature counting. UMI removal is powered by UMI-tools ([Smith, T., et al. 2017. Genome Res. 27, 491–499.](http://www.genome.org/cgi/doi/10.1101/gr.209601.116)). Use `-us` to set the UMI starting position on read 1. Use `-ul` to set the UMI length. Fragments with UMI length less than this value are discarded. UMI deduplication method is set by `-ud`.

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
2020-09-23 21:15:53,528 - fba.map - INFO - bowtie2 version: 2.4.1
2020-09-23 21:15:53,559 - fba.map - INFO - samtools version: 1.3
2020-09-23 21:15:53,614 - fba.map - INFO - Read 1 coordinates to search: [0, 16]
2020-09-23 21:15:53,614 - fba.map - INFO - Cell barcode maximum number of mismatches: 1
2020-09-23 21:15:53,614 - fba.map - INFO - Mapping quality threshold: 10
2020-09-23 21:15:53,614 - fba.map - INFO - Number of threads: 48
2020-09-23 21:15:53,614 - fba.map - INFO - Chunk size: 10,000
2020-09-23 21:15:53,614 - fba.map - INFO - UMI-tools version: 1.0.1
2020-09-23 21:15:53,614 - fba.map - INFO - UMI-tools deduplication method: directional
2020-09-23 21:15:53,614 - fba.map - INFO - UMI-tools deduplication threshold: 1
2020-09-23 21:15:53,614 - fba.map - INFO - UMI length: 12
2020-09-23 21:15:53,614 - fba.map - INFO - UMI starting position on read 1: 16
2020-09-23 21:15:53,614 - fba.map - INFO - Preparing feature references ...
2020-09-23 21:15:56,154 - fba.map - INFO - Matching cell barcodes (read 1) ...
2020-09-23 21:52:48,071 - fba.map - INFO - Aligning (read 2) ...
2020-09-23 22:16:05,486 - fba.map - INFO -
28336522 reads; of these:
  28336522 (100.00%) were unpaired; of these:
    3590828 (12.67%) aligned 0 times
    24425882 (86.20%) aligned exactly 1 time
    319812 (1.13%) aligned >1 times
87.33% overall alignment rate
2020-09-23 22:16:05,487 - fba.map - INFO - Generating matrix (UMI deduplication) ...
2020-09-23 22:17:43,808 - fba.map - INFO - Number of cell barcodes detected: 3,377
2020-09-23 22:17:43,808 - fba.map - INFO - Number of features detected: 1,127
2020-09-23 22:17:43,811 - fba.map - INFO - Total UMIs after deduplication: 2,106,755
2020-09-23 22:17:43,821 - fba.map - INFO - Median number of UMIs per cell: 445.0
```
