
[![PyPI](https://img.shields.io/pypi/v/fba?logo=pypi&style=flat-square)](https://pypi.org/project/fba/)
[![License](https://img.shields.io/badge/license-MIT-green?logo=data%3Aimage%2Fsvg%2Bxml%3Bbase64%2CPD94bWwgdmVyc2lvbj0iMS4wIiBlbmNvZGluZz0idXRmLTgiPz4KPCEtLXphei0tPgo8c3ZnIHhtbG5zPSJodHRwOi8vd3d3LnczLm9yZy8yMDAwL3N2ZyIgaGVpZ2h0PSIxNjYiIHdpZHRoPSIzMjEiPgo8ZyBzdHJva2Utd2lkdGg9IjM1IiBzdHJva2U9IiNBMzFGMzQiPgo8cGF0aCBkPSJtMTcuNSwwdjE2Nm01Ny0xNjZ2MTEzbTU3LTExM3YxNjZtNTctMTY2djMzbTU4LDIwdjExMyIvPgo8cGF0aCBkPSJtMTg4LjUsNTN2MTEzIiBzdHJva2U9IiM4QThCOEMiLz4KPHBhdGggZD0ibTIyOSwxNi41aDkyIiBzdHJva2Utd2lkdGg9IjMzIi8%2BCjwvZz4KPC9zdmc%2BCg%3D%3D&style=flat-square)](https://github.com/jlduan/fba/blob/master/LICENSE)
[![GitHub Workflow Status (branch)](https://img.shields.io/github/workflow/status/jlduan/fba/ci/master?logo=github-actions&style=flat-square)](https://github.com/jlduan/fba/actions/workflows/ci.yml)
[![CircleCI](https://img.shields.io/circleci/build/github/jlduan/fba/master?logo=circleci&style=flat-square)](https://app.circleci.com/pipelines/github/jlduan/fba)
[![Read the Docs](https://img.shields.io/readthedocs/fba?logo=read-the-docs&style=flat-square)](https://fba.readthedocs.io/en/latest/index.html)
[![Codecov](https://img.shields.io/codecov/c/gh/jlduan/fba?logo=codecov&style=flat-square)](https://app.codecov.io/gh/jlduan/fba)
[![GitHub commits since latest release (by date)](https://img.shields.io/github/commits-since/jlduan/fba/latest?color=9cf&logo=git&logoColor=red&style=flat-square)](https://github.com/jlduan/fba/commits)
[![Zenodo DOI](https://img.shields.io/badge/DOI-10.5281%2Fzenodo.4642814-informational?logo=data:image/svg+xml;base64,PD94bWwgdmVyc2lvbj0iMS4wIiBlbmNvZGluZz0iVVRGLTgiIHN0YW5kYWxvbmU9Im5vIj8+CjxzdmcgeG1sbnM9Imh0dHA6Ly93d3cudzMub3JnLzIwMDAvc3ZnIiB3aWR0aD0iMTMwIiBoZWlnaHQ9IjEzMCI+CjxjaXJjbGUgc3R5bGU9ImZpbGw6I2ZjYjQyNSIgY3g9IjY1IiBjeT0iNjUiIHI9IjY0Ii8+CjxwYXRoIHN0eWxlPSJmaWxsOiMyMzFmMjAiIGQ9Im0gNDkuODE5MTI3LDg0LjU1OTE0OCAtMTEuODU0MzA0LDAgMCwtNC44MjU2NjUgYyAtMS4yMDM1OTQsMS41MTA4OTQgLTQuMDM1NTE1LDMuMDUxMDUzIC01LjI2NDcxNiwzLjc0MjQ4MyAtMi4xNTExMDEsMS4yMDM1ODUgLTUuMDcyMDY2LDEuOTg3MjI1IC03LjgxMjE2MSwxLjk4NzIyNSAtNC40MzAyNDYsMCAtOC4zNzM5MjUsLTEuMzk5NTM5IC0xMS44MzEwNTcsLTQuNDQ2OTI0IC00LjEyMjk0NjQsLTMuNjM2Mzg5IC02LjA2MDI0NTUsLTkuMTk1NzYgLTYuMDYwMjQ1NSwtMTUuMTg4MTEzIDAsLTYuMDk0NzkxIDIuMTEyNjkxMywtMTAuOTYwMzgxIDYuMzM4MDY0NSwtMTQuNTk2NzYgMy4zNTQ2OTUsLTIuODkzNzQ1IDcuNDU3MDg5LC01LjIwOTc5NSAxMS44MTA1MDUsLTUuMjA5Nzk1IDIuNTM1MjMxLDAgNS42NjE4MDcsMC4yMjczNjMgNy44ODk3MzgsMS4zMDI5MTMgMS4yODA0MTQsMC42MTQ2MDEgMy41NzI2MjgsMi4wNjA3MjEgNC45Mjk4NzIsMy40NjkxNzkgbCAwLC0yNS40MjAxNzcgMTEuODU0MzA0LDAgeiBtIC0xMi4xMTk5LC0xOC42OTI1ODQgYyAwLC0yLjI1MzUzOCAtMC42MTgyNTgsLTQuOTUxNTU1IC0yLjIwNTk3MywtNi41MTM2NjMgLTEuNTg3NzI0LC0xLjU4NzcyNCAtNC40NzQxNTMsLTIuOTk2MTgyIC02LjcyNzY5MSwtMi45OTYxODIgLTIuNTA5NjE1LDAgLTQuODM0NDc2LDEuODI1NTExIC02LjQ0NzgwNywzLjcyMDUzNSAtMS4zMDYwMzEsMS41MzY1MDEgLTEuOTU5MDQxLDMuOTA1MjY5IC0xLjk1OTA0MSw1Ljg3NzExNCAwLDEuOTcxODM1IDAuNzQwODE1LDQuMTY1MDA0IDIuMDQ2ODM2LDUuNzAxNTA1IDEuNTg3NzE0LDEuODk1MDI1IDMuMjk3OTg1LDMuMTkzNzM5IDUuODMzMjE2LDMuMTkzNzM5IDIuMjc5MTQ1LDAgNC45ODk5NjUsLTAuOTU2NjYyIDYuNTUyMDgzLC0yLjUxODc3IDEuNTg3NzE0LC0xLjU2MjEwOCAyLjkwODM3NywtNC4xODUxMzQgMi45MDgzNzcsLTYuNDY0Mjc4IHoiLz4KPHBhdGggc3R5bGU9ImZpbGw6I2ZmZiIgZD0ibSAxMDUuNDI3NjQsMjUuNjE3OTE4IGMgLTEuOTcxODQsMCAtMy42NDkxOSwwLjY5MTQyIC01LjAzMjA0LDIuMDc0MjcxIC0xLjM1NzI0NywxLjM1NzI0NSAtMi4wMzU4NjQsMy4wMjE3NzkgLTIuMDM1ODY0LDQuOTkzNjMzIDAsMS45NzE4MzUgMC42Nzg2MTcsMy42NDkxOTMgMi4wMzU4NjQsNS4wMzIwMzQgMS4zODI4NSwxLjM4Mjg2MSAzLjA2MDIsMi4wNzQyODEgNS4wMzIwNCwyLjA3NDI4MSAxLjk5NzQ0LDAgMy42NzQ3OSwtMC42Nzg2MjcgNS4wMzIwMywtMi4wMzU4NjEgMS4zODI4NSwtMS4zODI4NjEgMi4wNzQyOCwtMy4wNzMwMTIgMi4wNzQyOCwtNS4wNzA0NTQgMCwtMS45NzE4NTQgLTAuNjkxNDMsLTMuNjM2Mzg4IC0yLjA3NDI4LC00Ljk5MzYzMyAtMS4zODI4NSwtMS4zODI4NTEgLTMuMDYwMiwtMi4wNzQyNzEgLTUuMDMyMDMsLTIuMDc0MjcxIHogTSA3NC4yMTkzODMsNDUuNTA3OTIxIGMgLTcuMzIzOTkyLDAgLTEyLjk3MDYyNSwyLjI4MzAwOSAtMTYuOTM5OTIxLDYuODQ4OTQ5IC0zLjI3Nzg3NiwzLjc4MjQzOCAtNC45MTY4MDMsOC4xMTgyNTIgLTQuOTE2ODAzLDEzLjAwODQwNiAwLDUuNDMwNDgxIDEuNjI2MTI0LDEwLjAwOTgzNCA0Ljg3ODM4MywxMy43MzgyMzYgMy45NDM2ODksNC41Mzg5MTggOS40NzUwOTMsNi44MDg2MjIgMTYuNTk0MjEsNi44MDg2MjIgNy4wOTM1MTIsMCAxMi42MTIxMjIsLTIuMjY5NzA0IDE2LjU1NTgwMSwtNi44MDg2MjIgMy4yNTIyNTksLTMuNzI4NDAyIDQuODc4MzkzLC04LjE5OTMgNC44NzgzOTMsLTEzLjQxMzY0OCAwLC01LjE2MDMyMyAtMS42Mzg5MzgsLTkuNjA0NjAyIC00LjkxNjgwMywtMTMuMzMyOTk0IC00LjAyMDUwOSwtNC41NjU5NCAtOS4zOTgyNjMsLTYuODQ4OTQ5IC0xNi4xMzMyNiwtNi44NDg5NDkgeiBtIDI0LjkwODYwMywxLjM4NjY4NiAwLDM3LjYzNDY3NiAxMi41OTkzMDQsMCAwLC0zNy42MzQ2NzYgLTEyLjU5OTMwNCwwIHogTSA3My44MzUyNTIsNTYuOTc1OTgxIGMgMi4zMDQ3NTIsMCA0LjI2Mzc5MywwLjg1MjMzNyA1Ljg3NzEyNCwyLjU1NDQyNiAxLjYzODkyOCwxLjY3NTA3NiAyLjQ1ODQwMiwzLjcyNzg4MSAyLjQ1ODQwMiw2LjE1OTQ1NyAwLDIuNDU4NTc4IC0wLjgwNjY3MSw0LjUzODAyMiAtMi40MTk5OTIsNi4yNDAxMTEgLTEuNjEzMzMxLDEuNjc1MDg2IC0zLjU4NTE3NSwyLjUxNDA5OSAtNS45MTU1MzQsMi41MTQwOTkgLTIuNjEyMDUxLDAgLTQuNzM3NTQ2LC0xLjAyNzM2NiAtNi4zNzY0NzQsLTMuMDgwNjgyIC0xLjMzMTYzNywtMS42NDgwNTMgLTEuOTk3NDUxLC0zLjUzOTE1NCAtMS45OTc0NTEsLTUuNjczNTI4IDAsLTIuMTA3MzYyIDAuNjY1ODE0LC0zLjk4NTEzOCAxLjk5NzQ1MSwtNS42MzMyMDEgMS42Mzg5MjgsLTIuMDUzMzE2IDMuNzY0NDIzLC0zLjA4MDY4MiA2LjM3NjQ3NCwtMy4wODA2ODIgeiIvPgo8L3N2Zz4K&style=flat-square)](https://doi.org/10.5281/zenodo.4642814)

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

optional arguments:
  -h, --help        show this help message and exit

functions:

    extract         extract cell and feature barcodes
    map             map enriched transcripts
    filter          filter extracted barcodes
    count           count feature barcodes per cell
    demultiplex     demultiplex cells based on feature abundance
    qc              quality control of feature barcoding assay
    kallisto_wrapper
                    deploy kallisto/bustools for feature barcoding
                    quantification
```

<br>

- __extract__: extract cell and feature barcodes from paired fastq files. For single cell assays, read 1 usually contains cell partitioning and UMI information, and read 2 contains feature information.
- __map__: quantify enriched transcripts (through hybridization or PCR amplification) from parent single cell libraries. Read 1 contains cell partitioning and UMI information, and read 2 contains transcribed regions of enriched/targeted transcripts of interest.  BWA (Li, H. 2013) or Bowtie2 (Langmead, B., et al. 2012) is used for read 2 alignment. The quantification (UMI deduplication) of enriched/targeted transcripts is powered by UMI-tools (Smith, T., et al. 2017).
- __filter__: filter extracted cell and feature barcodes (output of `extract` or `qc`). Additional fragment filter/selection can be applied through `-cb_seq` and/or `-fb_seq`.
- __count__: count UMIs per feature per cell (UMI deduplication), powered by UMI-tools (Smith, T., et al. 2017). Take the output of `extract` or `filter` as input.
- __demultiplex__: demultiplex cells based on the abundance of features (matrix generated by `count` as input).
- __qc__: generate diagnostic information. If `-1` is omitted, bulk mode is enabled and only read 2 will be analyzed.
- __kallisto_wrapper__: deploy kallisto/bustools for feature barcoding quantification (just a wrapper) (Bray, N.L., et al. 2016).

<br>

## Workflow example

- Cell surface protein labeling
    - [CITE-Seq; 8k cord blood mononuclear cells with 13 antibodies](https://fba.readthedocs.io/en/latest/tutorials/cell_surface_protein_labeling/PRJNA393315/tutorial.html#tutorial-cell-surface-protein-labeling-prjna393315)
    - [1k Human PBMCs Stained with a Panel of TotalSeq B Antibodies, Dual Indexed](https://fba.readthedocs.io/en/latest/tutorials/cell_surface_protein_labeling/SC3_v3_NextGem_DI_PBMC_CSP_1K/tutorial.html#tutorial-cell-surface-protein-labeling-sc3-v3-nextgem-di-pbmc-csp-1k)
- Cell hashing
    - [Peripheral blood mononuclear cells with 8 antibodies](https://fba.readthedocs.io/en/latest/tutorials/cell_hashing/PRJNA423077/tutorial.html#tutorial-cell-hashing-prjna423077)
- CRISPR screening
    - [10k A375 Cells Transduced with (1) Non-Target and (1) Target sgRNA, Dual Indexed](https://fba.readthedocs.io/en/latest/tutorials/crispr_screening/SC3_v3_NextGem_DI_CRISPR_10K/tutorial.html#tutorial-crispr-screening-sc3-v3-nextgem-di-crispr-10k)
- Targeted transcript enrichment
    - [Hodgkin's Lymphoma, Dissociated Tumor: Targeted, Gene Signature Panel](https://fba.readthedocs.io/en/latest/tutorials/targeted_transcript_enrichment/Targeted_NGSC3_DI_HodgkinsLymphoma_GeneSignature/tutorial.html#tutorial-targeted-transcript-enrichment-targeted-ngsc3-di-hodgkinslymphoma-genesignature)
- Bulk
    - [10k A375 Cells Transduced with (1) Non-Target and (1) Target sgRNA, Dual Indexed](https://fba.readthedocs.io/en/latest/tutorials/bulk/SC3_v3_NextGem_DI_CRISPR_10K/tutorial.html#tutorial-bulk-sc3-v3-nextgem-di-crispr-10k)

<br>
