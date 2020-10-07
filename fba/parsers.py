# parser.py

import sys
import argparse
import numpy as np
from itertools import chain
from fba import __version__


def coords(s):
    try:
        s = ''.join([i.strip() for i in chain.from_iterable(s)]).split(',')
        return [int(i) for i in s]

    except TypeError:
        raise argparse.ArgumentTypeError('Coordinates must be start,end')


def parse_args(args=sys.argv[1:]):

    # create top-level parser
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=(
            'Tools for feature barcoding analyses\n'
            + "Version: "
            + __version__
        )
    )

    # create sub-level parser
    subparsers = parser.add_subparsers(
        title='functions',
        dest='command',
        metavar=''
    )

    add_extract_subparser(subparsers)
    add_map_subparser(subparsers)
    add_filter_subparser(subparsers)
    add_count_subparser(subparsers)
    add_demultiplex_subparser(subparsers)
    add_qc_subparser(subparsers)
    add_kallisto_subparser(subparsers)

    if len(sys.argv) > 1:
        if (sys.argv[1] == '--version' or sys.argv[1] == '-v'):
            print("fba version: %s" % __version__)
            exit()
    else:
        args = parser.parse_args(['-h'])
        exit()

    return parser.parse_args(args)


# extract
def add_extract_subparser(subparsers):
    parser = subparsers.add_parser(
        'extract',
        help='extract feature barcodes'
    )

    parser.add_argument(
        '-1',
        '--read1',
        dest='read1',
        required=True,
        type=str,
        help='specify fastq file for read 1'
    )

    parser.add_argument(
        '-2',
        '--read2',
        dest='read2',
        required=True,
        type=str,
        help='specify fastq file for read 2'
    )

    parser.add_argument(
        '-w', '--whitelist',
        dest='whitelist',
        required=True,
        type=str,
        help='specify a whitelist of accepted cell barcodes'
    )

    parser.add_argument(
        '-f', '--feature_ref',
        dest='feature_ref',
        required=True,
        type=str,
        help='specify a reference of feature barcodes'
    )

    parser.add_argument(
        '-o', '--output',
        dest='output',
        required=True,
        help='specify an output file'
    )

    parser.add_argument(
        '-r1_coords', '--read1_coords',
        dest='read1_coords',
        required=False,
        default=(0, 16),
        type=coords,
        help='specify expected coordinates \'start,end\' of read 1 to search. For example, \'0,16\': starts at 0, stops at 15. Nucleotide bases outside the range will be masked as lower case in output. The default is \'0,16\''
    )

    parser.add_argument(
        '-r2_coords', '--read2_coords',
        dest='read2_coords',
        required=True,
        default=None,
        type=coords, help='see \'-r1_coords\''
    )
    parser.add_argument(
        '-cb_m',
        '--cb_mismatches',
        dest='cell_barcode_mismatches',
        required=False,
        type=int,
        default=1,
        help='specify cell barcode mismatching threshold. The default is 1'
    )

    parser.add_argument(
        '-fb_m',
        '--fb_mismatches',
        dest='feature_barcode_mismatches',
        required=False,
        type=int,
        default=1,
        help='specify feature barcode mismatching threshold. The default is 1'
    )

    parser.add_argument(
        '-cb_n',
        '--cb_num_n_threshold',
        dest='cb_num_n_threshold',
        required=False,
        type=int,
        default=3,
        help='specify maximum number of ambiguous nucleotides allowed for read 1. The default is 3'
    )

    parser.add_argument(
        '-fb_n',
        '--fb_num_n_threshold',
        dest='fb_num_n_threshold',
        required=False,
        type=int,
        default=3,
        help='specify maximum number of ambiguous nucleotides allowed for read 2. The default is 3'
    )

    parser.add_argument(
        '-e',
        '--exhaustive',
        dest='exhaustive',
        required=False,
        action='store_true',
        help='specify whether to search all the barcodes meeting the mismatching criteria and select the best one. If not specified, the first one will be selected'
    )


# map
def add_map_subparser(subparsers):
    parser = subparsers.add_parser(
        'map',
        help='map enriched transcripts',
        description='Quantify enriched transcripts (hybridization or PCR amplification) from single cell libraries.'
    )

    parser.add_argument(
        '-1',
        '--read1',
        dest='read1',
        required=True,
        type=str,
        help='specify fastq file for read 1'
    )

    parser.add_argument(
        '-2',
        '--read2',
        dest='read2',
        required=True,
        type=str,
        help='specify fastq file for read 2'
    )

    parser.add_argument(
        '-w', '--whitelist',
        dest='whitelist',
        required=True,
        type=str,
        help='specify a whitelist of accepted cell barcodes'
    )

    parser.add_argument(
        '-f', '--feature_ref',
        dest='feature_ref',
        required=True,
        type=str,
        help='specify a reference of feature barcodes'
    )

    parser.add_argument(
        '-o', '--output',
        dest='output',
        required=True,
        help='specify an output file'
    )

    parser.add_argument(
        '-r1_coords', '--read1_coords',
        dest='read1_coords',
        required=False,
        default=(0, 16),
        type=coords,
        help='specify expected coordinates \'start,end\' of read 1 to search for cell barcodes. The default is \'0,16\''
    )

    parser.add_argument(
        '-cb_m',
        '--cb_mismatches',
        dest='cell_barcode_mismatches',
        required=False,
        type=int,
        default=1,
        help='specify cell barcode mismatching threshold. The default is 1'
    )

    parser.add_argument(
        '-cb_n',
        '--cb_num_n_threshold',
        dest='cb_num_n_threshold',
        required=False,
        type=int,
        default=3,
        help='specify maximum number of ambiguous nucleotides allowed for read 1. The default is 3'
    )

    parser.add_argument(
        '--mapq', '--mapping_quality',
        dest='mapq',
        required=False,
        default=10,
        type=int,
        help='specify minimal mapping quality required for feature mapping. The default is 10'
    )

    parser.add_argument(
        '-us',
        '--umi_start',
        dest='umi_pos_start',
        required=False,
        type=int,
        default=16,
        help='specify expected UMI starting postion on read 1. The default is 16'
    )

    parser.add_argument(
        '-ul',
        '--umi_length',
        dest='umi_length',
        required=False,
        type=int,
        default=12,
        help='specify the length of UMIs on read 1. Reads with UMI length less than this value will be discarded. The default is 12'
    )

    parser.add_argument(
        '-um',
        '--umi_mismatches',
        dest='umi_mismatches',
        required=False,
        type=int,
        default=1,
        help='specify the maximun edit distance allowed for UMIs on read 1 for deduplication. The default is 1'
    )

    parser.add_argument(
        '-ud',
        '--umi_deduplication_method',
        dest='umi_deduplication_method',
        required=False,
        type=str,
        choices=['unique', 'percentile',
                 'cluster', 'adjacency', 'directional'],
        default='directional',
        help='specify UMI deduplication method (powered by UMI-tools. Smith, T., et al. 2017). The default is \'directional\''
    )

    parser.add_argument(
        '--output_directory',
        dest='output_directory',
        required=False,
        type=str,
        default='barcode_mapping',
        help='specify a temp directory. The default is \'./barcode_mapping\''
    )

    parser.add_argument(
        '-t', '--threads',
        dest='threads',
        required=False,
        type=int,
        default=None,
        help='specify number of threads to launch. The default is to use all available'
    )

    parser.add_argument(
        '--num_n_ref',
        dest='num_n_ref',
        required=False,
        type=int,
        default=0,
        help='specify the number of Ns to separate sequences belonging to the same feature. The default is 0'
    )


# filter
def add_filter_subparser(subparsers):
    parser = subparsers.add_parser(
        'filter',
        help='filter extracted feature barcodes')

    parser.add_argument(
        '-i',
        '--input',
        dest='input',
        required=True,
        type=str,
        help='specify an input file. The output of extract or qc subcommand'
    )

    parser.add_argument(
        '-o',
        '--output',
        dest='output',
        required=True,
        type=str,
        help='specify an output file'
    )

    parser.add_argument(
        '-cb_s',
        '--cb_start',
        dest='cell_barcode_pos_start',
        required=False,
        type=int,
        default=0,
        help='specify expected cell barcode starting postion on read 1. The default is 0'
    )

    parser.add_argument(
        '-cb_m',
        '--cb_mismatches',
        dest='cell_barcode_mismatches',
        required=False,
        type=int,
        default=1,
        help='specify cell barcode mismatching threshold. The default is 1'
    )

    parser.add_argument(
        '-cb_ls',
        '--cb_left_shift',
        dest='cell_barcode_left_shift',
        required=False,
        type=int,
        default=1,
        help='specify the maximum left shift allowed for cell barcode. The default is 1'
    )

    parser.add_argument(
        '-cb_rs',
        '--cb_right_shift',
        dest='cell_barcode_right_shift',
        required=False,
        type=int,
        default=1,
        help='specify the maximum right shift allowed for cell barcode. The default is 1'
    )

    parser.add_argument(
        '-cb_seq',
        '--cb_extra_seq',
        dest='cell_barcode_extra_seq',
        required=False,
        type=str,
        default=None,
        help='specify an extra constant sequence to filter on read 1. The default is None'
    )

    parser.add_argument(
        '-cb_seq_m',
        '--cb_extra_seq_mismatches',
        dest='cell_barcode_extra_seq_mismatches',
        required=False,
        type=int,
        default=None,
        help='specify the maximun edit distance allowed for the extra constant sequence on read 1 for filtering. The default is off'
    )

    parser.add_argument(
        '-fb_s',
        '--fb_start',
        dest='feature_barcode_pos_start',
        required=False,
        type=int,
        default=10,
        help='specify expected feature barcode starting postion on read 2. The default is 10'
    )

    parser.add_argument(
        '-fb_m',
        '--fb_mismatches',
        dest='feature_barcode_mismatches',
        required=False,
        type=int,
        default=1,
        help='specify feature barcode mismatching threshold. The default is 1'
    )

    parser.add_argument(
        '-fb_ls',
        '--fb_left_shift',
        dest='feature_barcode_left_shift',
        required=False,
        type=int,
        default=1,
        help='specify the maximum left shift allowed for feature barcode. The default is 1'
    )

    parser.add_argument(
        '-fb_rs',
        '--fb_right_shift',
        dest='feature_barcode_right_shift',
        required=False,
        type=int,
        default=1,
        help='specify the maximum right shift allowed for feature barcode. The default is 1'
    )

    parser.add_argument(
        '-fb_seq',
        '--fb_extra_seq',
        dest='feature_barcode_extra_seq',
        required=False,
        type=str,
        default=None,
        help='specify an extra constant sequence to filter on read 2. The default is None'
    )

    parser.add_argument(
        '-fb_seq_m',
        '--fb_extra_seq_mismatches',
        dest='feature_barcode_extra_seq_mismatches',
        required=False,
        type=int,
        default=None,
        help='specify the maximun edit distance allowed for the extra constant sequence on read 2. The default is off'
    )


# count
def add_count_subparser(subparsers):
    parser = subparsers.add_parser(
        'count',
        help='count feature barcodes per cell'
    )

    parser.add_argument(
        '-i',
        '--input',
        dest='input',
        required=True,
        type=str,
        help='specify an input file'
    )

    parser.add_argument(
        '-o',
        '--output',
        dest='output',
        required=True,
        type=str,
        help='specify an output file'
    )

    parser.add_argument(
        '-us',
        '--umi_start',
        dest='umi_pos_start',
        required=False,
        type=int,
        default=16,
        help='specify expected UMI starting postion on read 1. The default is 16'
    )

    parser.add_argument(
        '-ul',
        '--umi_length',
        dest='umi_length',
        required=False,
        type=int,
        default=12,
        help='specify the length of UMIs on read 1. Reads with UMI length less than this value will be discarded. The default is 12'
    )
    parser.add_argument(
        '-um',
        '--umi_mismatches',
        dest='umi_mismatches',
        required=False,
        type=int,
        default=1,
        help='specify the maximun edit distance allowed for UMIs on read 1 for deduplication. The default is 1'
    )

    parser.add_argument(
        '-ud',
        '--umi_deduplication_method',
        dest='umi_deduplication_method',
        required=False,
        type=str,
        choices=['unique', 'percentile',
                 'cluster', 'adjacency', 'directional'],
        default='directional',
        help='specify UMI deduplication method (powered by UMI-tools. Smith, T., et al. 2017). The default is \'directional\''
    )


# demultiplex
def add_demultiplex_subparser(subparsers):
    parser = subparsers.add_parser(
        'demultiplex',
        help='demultiplex cells based on feature barcoding'
    )

    parser.add_argument(
        '-i',
        '--input',
        dest='input',
        required=True,
        type=str,
        help='specify an input file (feature count matrix). The output of count subcommand'
    )

    parser.add_argument(
        '--output_directory',
        dest='output_directory',
        required=False,
        type=str,
        default='demultiplexed',
        help='specify a output directory. The default is \'./demultiplexed\''
    )

    parser.add_argument(
        '-m',
        dest='method',
        required=False,
        type=int,
        default=1,
        choices=[1],
        help='specify demultiplexing method. \'1\': Stoeckius et al. 2018. The default is 1'
    )

    parser.add_argument(
        '-q',
        dest='quantile',
        required=False,
        type=float,
        default=0.9999,
        help='specify quantile for the probability mass function (0.9 to 1 recommended). The default is 0.9999'
    )


# qc
def add_qc_subparser(subparsers):
    parser = subparsers.add_parser(
        'qc',
        help='quality control of feature barcoding'
    )

    parser.add_argument(
        '-1',
        '--read1',
        dest='read1',
        required=False,
        type=str,
        help='specify fastq file for read 1'
    )

    parser.add_argument(
        '-2',
        '--read2',
        dest='read2',
        required=True,
        type=str,
        help='specify fastq file for read 2. If only read 2 file is provided, bulk mode is enabled (skipping arguments \'-1\', \' -w\', \'-cb_m\', \'-r1_coords\', must provide \'-r2_coords\' and \'-fb_m\'). In bulk mode, reads 2 will be searched against reference feature barcodes and read count for each feature barcode will be summarized'
    )

    parser.add_argument(
        '-w', '--whitelist',
        dest='whitelist',
        required=False,
        type=str,
        help='specify a whitelist of accepted cell barcodes'
    )

    parser.add_argument(
        '-f', '--feature_ref',
        dest='feature_ref',
        required=True,
        type=str,
        help='specify a reference of feature barcodes'
    )

    parser.add_argument(
        '-r1_coords', '--read1_coords',
        dest='read1_coords',
        required=False,
        default=None,
        type=coords,
        help='specify coordinates \'start,end\' of read 1 to search (doesn\'t need to be the exact expected barcode range). The default is to use all the nucleotide bases. Nucleotide bases outside the range will be masked as lower case in output'
    )

    parser.add_argument(
        '-r2_coords', '--read2_coords',
        dest='read2_coords',
        required=False,
        default=None,
        type=coords, help='see \'-r1_coords\''
    )

    parser.add_argument(
        '-cb_m',
        '--cb_mismatches',
        dest='cell_barcode_mismatches',
        required=False,
        type=int,
        default=3,
        help='specify cell barcode mismatching threshold. The default is 3'
    )

    parser.add_argument(
        '-fb_m',
        '--fb_mismatches',
        dest='feature_barcode_mismatches',
        required=False,
        type=int,
        default=3,
        help='specify feature barcode mismatching threshold. The default is 3'
    )

    parser.add_argument(
        '-cb_n',
        '--cb_num_n_threshold',
        dest='cb_num_n_threshold',
        required=False,
        type=int,
        default=np.Inf,
        help='specify maximum number of ambiguous nucleotides allowed for read 1. The default is no limit'
    )

    parser.add_argument(
        '-fb_n',
        '--fb_num_n_threshold',
        dest='fb_num_n_threshold',
        required=False,
        type=int,
        default=np.Inf,
        help='specify maximum number of ambiguous nucleotides allowed for read 2. The default is no limit'
    )

    parser.add_argument(
        '-t', '--threads',
        dest='threads',
        required=False,
        type=int,
        default=None,
        help='specify number of threads for barcode extraction. The default is to use all available'
    )

    parser.add_argument(
        '-n', '--num_reads',
        dest='num_reads',
        required=False,
        # type=int,
        default=100_000,
        help='specify number of reads for analysis. Set to \'None\' will analyze all the reads. The default is 100,000'
    )

    parser.add_argument(
        '--chunk_size',
        dest='chunk_size',
        required=False,
        type=int,
        default=50_000,
        help='specify the chunk size for multiprocessing. The default is 50,000'
    )

    parser.add_argument(
        '--output_directory',
        dest='output_directory',
        required=False,
        type=str,
        default='qc',
        help='specify a output directory. The default is \'./qc\''
    )


# kallisto
def add_kallisto_subparser(subparsers):
    parser = subparsers.add_parser(
        'kallisto_wrapper',
        help='deploy kallisto/bustools for feature barcoding quantification',
        description='kallisto'
    )

    parser.add_argument(
        '-1',
        '--read1',
        dest='read1',
        required=True,
        type=str,
        help='specify fastq file for read 1'
    )

    parser.add_argument(
        '-2',
        '--read2',
        dest='read2',
        required=True,
        type=str,
        help='specify fastq file for read 2'
    )

    parser.add_argument(
        '-w', '--whitelist',
        dest='whitelist',
        required=True,
        type=str,
        help='specify a whitelist of accepted cell barcodes'
    )

    parser.add_argument(
        '-f', '--feature_ref',
        dest='feature_ref',
        required=True,
        type=str,
        help='specify a reference of feature barcodes'
    )

    parser.add_argument(
        '--technology',
        dest='technology',
        required=False,
        type=str,
        default='10xv3',
        help='specify feature barcoding technology. The default is 10xv3'
    )

    parser.add_argument(
        '-o', '--output',
        dest='output',
        required=True,
        help='specify an output file'
    )

    parser.add_argument(
        '-t', '--threads',
        dest='threads',
        required=False,
        type=int,
        default=1,
        help='specify number of kallisto/bustools threads to launch. The default is 1'
    )

    parser.add_argument(
        '--output_directory',
        dest='output_directory',
        required=False,
        type=str,
        default='barcode_mapping',
        help='specify a temp directory. The default is \'./kallisto\''
    )
