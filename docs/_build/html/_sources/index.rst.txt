:og:description:
   fba is a flexible and streamlined toolbox for quality control,
   quantification, demultiplexing of various single-cell feature
   barcoding assays. It can be applied to customized feature barcoding
   specifications, including different CRISPR constructs or targeted
   enriched transcripts.

:og:image:
   https://oup.silverchair-cdn.com/oup/backfile/Content_public/Journal/bioinformatics/37/22/10.1093_bioinformatics_btab375/2/btab375f1.jpeg?Expires=1674767963&Signature=fs5S-bzHNkgN5tWD2FuCC3dMTaSCMqxLj9BcETS0QfJ~Dd2kJJ5h9j7auTFwiTiqJDaBBtWwBFAvlKKEkmkEOzpD-2tcsgnenYqaCKK1orz6nI2OT0pG50npsL8YpjgSuNcRnwQb4rEEW67T2qBOxymLbdMB-Qc5R2Z1YSYjoJKz2klCQYboohwP6eZg4g2zhPgyJDK7FDvQqfNwYDtJWSGf2~L6RltKrmxs6qeiT23r8lDOxdun1KC6dJqOS8MIhQoZ2err7SE-JXqdAmVH89AmQZyALDsstb2ZqSEyw5KyPWidgIc3PPeCaSLNWydW7LIDztZZzyLRwQH59UWZTg__&Key-Pair-Id=APKAIE5G5CRDK6RD3PGA

..
   fba documentation master file, created by
   sphinx-quickstart on Tue Apr  6 18:44:42 2021.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

#####################################
 Welcome to ``fba``'s documentation!
#####################################

``fba``: a flexible and streamlined package for single-cell feature
barcoding assays.

.. image:: https://img.shields.io/pypi/v/fba.svg?logo=pypi
   :target: https://pypi.python.org/pypi/fba

.. image:: https://img.shields.io/conda/vn/bioconda/fba?logo=anaconda
   :alt: Conda (channel only)
   :target: https://bioconda.github.io/recipes/fba/README.html

.. image:: https://img.shields.io/github/v/release/jlduan/fba?color=orange&logo=github
   :alt: GitHub release (latest by date)
   :target: https://github.com/jlduan/fba/releases

.. image:: https://img.shields.io/badge/license-MIT-yellow.svg
   :target: https://github.com/jlduan/fba/blob/master/LICENSE

.. image:: https://github.com/jlduan/fba/actions/workflows/ci.yml/badge.svg?branch=master
   :target: https://github.com/jlduan/fba/actions/workflows/ci.yml

.. image:: https://circleci.com/gh/jlduan/fba/tree/master.svg?style=svg
   :target: https://circleci.com/gh/jlduan/fba/tree/master

.. image:: https://readthedocs.org/projects/fba/badge/?version=latest
   :target: https://fba.readthedocs.io/en/latest/

.. image:: https://codecov.io/gh/jlduan/fba/branch/master/graph/badge.svg?token=H3189R59G0
   :target: https://codecov.io/gh/jlduan/fba

.. image:: https://api.codeclimate.com/v1/badges/d52a1e18eb229ca39da2/maintainability
   :alt: Maintainability
   :target: https://codeclimate.com/github/jlduan/fba/maintainability

.. image:: https://static.pepy.tech/personalized-badge/fba?period=total&units=international_system&left_color=grey&right_color=yellow&left_text=downloads
   :target: https://pypi.org/project/fba

.. image:: https://img.shields.io/conda/dn/bioconda/fba?logo=anaconda
   :alt: Conda
   :target: https://bioconda.github.io/recipes/fba/README.html

.. image:: https://zenodo.org/badge/292596737.svg
   :target: https://zenodo.org/badge/latestdoi/292596737

.. toctree::
   :maxdepth: 3
   :hidden:

   installation
   user_guide
   tutorials/index

.. toctree::
   :caption: Development
   :hidden:

   changelog
   GitHub Repository <https://github.com/jlduan/fba>

..

   **工欲善其事，必先利其器。—— 论语·卫灵公**

************
 Quickstart
************

Installation
============

.. code:: console

   $ pip install fba

Usage
=====

.. code:: console

   $ fba

   usage: fba [-h]  ...

   Tools for feature barcoding analysis

   optional arguments:
   -h, --help         show this help message and exit

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

**********
 Citation
**********

Jialei Duan, Gary C Hon, **FBA: feature barcoding analysis for single
cell RNA-Seq**, *Bioinformatics*, Volume 37, Issue 22, 15 November 2021,
Pages 4266–4268. DOI: `10.1093/bioinformatics/btab375`_. PMID:
33999185_.

.. _10.1093/bioinformatics/btab375: https://doi.org/10.1093/bioinformatics/btab375

.. _33999185: https://pubmed.ncbi.nlm.nih.gov/33999185

.. raw:: html

   <span class="__dimensions_badge_embed__" data-doi="10.1093/bioinformatics/btab375" data-hide-zero-citations="true" data-style="small_rectangle"></span><script async src="https://badge.dimensions.ai/badge.js" charset="utf-8"></script>
