.. fba documentation master file, created by
   sphinx-quickstart on Tue Apr  6 18:44:42 2021.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to fba's documentation!
===============================


``fba``: Tools for feature barcoding analysis


.. image:: https://img.shields.io/pypi/v/fba.svg
   :target: https://pypi.python.org/pypi/fba

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

.. image:: https://zenodo.org/badge/DOI/10.5281/zenodo.4642815.svg
   :target: https://doi.org/10.5281/zenodo.4642815


Quickstart
----------


Installation
^^^^^^^^^^^^

.. code-block::

   $ pip install fba

Usage
^^^^^

.. code-block::

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


Contents
--------

.. :caption: Contents

.. toctree::
   :maxdepth: 3

   installation
   user_guide

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
