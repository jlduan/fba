
[![PyPI](https://img.shields.io/pypi/v/fba?logo=pypi&style=flat-square)](https://pypi.org/project/fba/) [![PyPI - License](https://img.shields.io/pypi/l/fba?style=flat-square)](https://github.com/jlduan/fba/blob/master/LICENSE)

# fba

Tools for feature barcoding analysis

<br>

## Installation

```shell
$ pip install fba
```

<br>

## Usage

```
$ fba
usage: fba [-h]  ...

Tools for feature barcoding analyses
Version: 0.0.5dev

optional arguments:
  -h, --help        show this help message and exit

functions:

    extract         extract feature barcodes
    map             map enriched transcripts
    filter          filter extracted feature barcodes
    count           count feature barcodes per cell
    demultiplex     demultiplex cells based on feature barcoding
    qc              quality control of feature barcoding
    kallisto_wrapper
                    deploy kallisto/bustools for feature barcoding
                    quantification
```

<br>

- __extract__: extract cell and feature barcodes from paired fastq files. For single cell RNA-Seq assay, reads 1 usually contain cell barcodes and UMIs, and reads 2 contain feature barcodes.

- __map__: quantify enriched transcripts (through hybridization or PCR amplification) from parent single cell libraries. Reads 1 contain cell barcodes and UMIs, and reads 2 are transcribed regions of enriched/targeted transcripts of interest. Bowtie2 (Langmead, B., et al. 2012) is used for reads 2 alignment. The quantification (UMI deduplication) of enriched/targeted transcripts is powered by UMI-tools (Smith, T., et al. 2017).

- __filter__: filter extracted cell and feature barcoding results (output of `extract` or `qc`).  Additional fragment selection can be performed through `-cb_seq` and/or `-fb_seq`.

- __count__: count UMIs per feature per cell (UMI deduplication), powered by UMI-tools (Smith, T., et al. 2017). Take the output of `extract` or `filter` as input.

- __demultiplex__: demultiplex cells based on the abundance of features (matrix generated by `count`).

- __qc__: generate diagnostic results through fastq files. If `-1` is omitted, bulk mode is enabled and only reads 2 will be analyzed.

- __kallisto_wrapper__: deploy kallisto/bustools for feature barcoding quantification (just a wrapper) (Bray, N.L., et al. 2016).


<br>

## Example workflow

- Cell surface protein labeling
    - [CITE-Seq; 8k cord blood mononuclear cells with 13 antibodies](https://github.com/jlduan/fba/blob/master/examples/cell_surface_protein_labeling/PRJNA393315/tutorial.md)
    - [1k Human PBMCs Stained with a Panel of TotalSeq B Antibodies, Dual Indexed](https://github.com/jlduan/fba/blob/master/examples/cell_surface_protein_labeling/SC3_v3_NextGem_DI_PBMC_CSP_1K/tutorial.md)

- Cell hashing
    - [Peripheral blood mononuclear cells with 8 antibodies](https://github.com/jlduan/fba/blob/master/examples/cell_hashing/PRJNA423077/tutorial.md)

- CRISPR screening
    - [10k A375 Cells Transduced with (1) Non-Target and (1) Target sgRNA, Dual Indexed](https://github.com/jlduan/fba/blob/master/examples/crispr_screening/SC3_v3_NextGem_DI_CRISPR_10K/tutorial.md)

- Targeted transcript enrichment
    - [Hodgkin's Lymphoma, Dissociated Tumor: Targeted, Gene Signature Panel](https://github.com/jlduan/fba/blob/master/examples/targeted_transcript_enrichment/Targeted_NGSC3_DI_HodgkinsLymphoma_GeneSignature/tutorial.md)

- Bulk
    - [10k A375 Cells Transduced with (1) Non-Target and (1) Target sgRNA, Dual Indexed](https://github.com/jlduan/fba/blob/master/examples/bulk/SC3_v3_NextGem_DI_CRISPR_10K/tutorial.md)

<br>
