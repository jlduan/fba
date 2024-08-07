# parsers.py

import argparse
import sys

from fba import __version__


def coords(s):
    try:
        # from itertools import chain
        # s = ''.join([i.strip() for i in chain.from_iterable(s)]).split(',')
        # return [int(i) for i in s]
        return [int(i) for i in s.split(",")]

    except TypeError:
        raise argparse.ArgumentTypeError("Coordinate format must be start,end")


def parse_args(args=sys.argv[1:]):
    # create top-level parser
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=(
            "Tools for single-cell feature barcoding analysis\n"
            + "Version: "
            + __version__
        ),
    )

    # create sub-level parser
    subparsers = parser.add_subparsers(title="functions", dest="command", metavar="")

    add_extract_subparser(subparsers)
    add_map_subparser(subparsers)
    add_filter_subparser(subparsers)
    add_count_subparser(subparsers)
    add_demultiplex_subparser(subparsers)
    add_qc_subparser(subparsers)
    add_kallisto_subparser(subparsers)

    if len(sys.argv) > 1:
        if sys.argv[1] == "--version" or sys.argv[1] == "-v":
            print("fba version: %s" % __version__)
            exit()
    else:
        args = parser.parse_args(["-h"])
        exit()

    return parser.parse_args(args)


# extract
def add_extract_subparser(subparsers):
    parser = subparsers.add_parser(
        "extract",
        help="extract cell and feature barcodes",
        description=(
            "Extract cell and feature barcodes from paired fastq files. "
            "For single cell assays, "
            "read 1 usually contains cell partitioning and UMI information, "
            "and read 2 contains feature information."
        ),
    )

    parser.add_argument(
        "-1",
        "--read1",
        dest="read1",
        required=True,
        type=str,
        help="specify fastq file for read 1",
    )

    parser.add_argument(
        "-2",
        "--read2",
        dest="read2",
        required=True,
        type=str,
        help="specify fastq file for read 2",
    )

    parser.add_argument(
        "-w",
        "--whitelist",
        dest="whitelist",
        required=True,
        type=str,
        help="specify a whitelist of accepted cell barcodes",
    )

    parser.add_argument(
        "-f",
        "--feature_ref",
        dest="feature_ref",
        required=True,
        type=str,
        help="specify a reference of feature barcodes",
    )

    parser.add_argument(
        "-o",
        "--output",
        dest="output",
        required=True,
        help="specify an output file",
    )

    parser.add_argument(
        "-r1_c",
        "--read1_coordinate",
        dest="read1_coordinate",
        required=False,
        default=(0, 16),
        type=coords,
        help=(
            "specify coordinate 'start,end' of read 1 to search. "
            "For example, '0,16': starts at 0, stops at 15. "
            "Nucleotide bases outside the range will be masked "
            "as lowercase in the output. Default (0,16)"
        ),
    )

    parser.add_argument(
        "-r2_c",
        "--read2_coordinate",
        dest="read2_coordinate",
        required=True,
        default=None,
        type=coords,
        help="see '-r1_c/--read1_coordinate'",
    )
    parser.add_argument(
        "-cb_m",
        "--cb_mismatches",
        dest="cell_barcode_mismatches",
        required=False,
        type=int,
        default=1,
        help="specify cell barcode mismatching threshold. Default (1)",
    )

    parser.add_argument(
        "-fb_m",
        "--fb_mismatches",
        dest="feature_barcode_mismatches",
        required=False,
        type=int,
        default=1,
        help="specify feature barcode mismatching threshold. Default (1)",
    )

    parser.add_argument(
        "-cb_n",
        "--cb_num_n_threshold",
        dest="cb_num_n_threshold",
        required=False,
        type=int,
        default=3,
        help=(
            "specify maximum number of ambiguous nucleotides "
            "allowed for read 1. Default (3)"
        ),
    )

    parser.add_argument(
        "-fb_n",
        "--fb_num_n_threshold",
        dest="fb_num_n_threshold",
        required=False,
        type=int,
        default=3,
        help=(
            "specify maximum number of ambiguous nucleotides "
            "allowed for read 2. Default (3)"
        ),
    )

    parser.add_argument(
        "-cb_rc",
        "--cell_barcode_reverse_complement",
        dest="cb_reverse_complement",
        required=False,
        action="store_true",
        help=(
            "specify to convert cell barcode sequences "
            "into their reverse-complement counterparts for processing."
        ),
    )


# map
def add_map_subparser(subparsers):
    parser = subparsers.add_parser(
        "map",
        help="map enriched transcripts",
        description=(
            "Quantify enriched transcripts "
            "(through hybridization or PCR amplification) "
            "from parent single cell libraries. "
            "Read 1 contains cell partitioning and UMI information, "
            "and read 2 contains transcribed regions of "
            "enriched/targeted transcripts of interest. "
            "BWA (Li, H. 2013) or Bowtie2 (Langmead, B., et al. 2012) is used "
            "for read 2 alignment. The quantification (UMI deduplication) "
            "of enriched/targeted transcripts is powered by "
            "UMI-tools (Smith, T., et al. 2017)."
        ),
    )

    parser.add_argument(
        "-1",
        "--read1",
        dest="read1",
        required=True,
        type=str,
        help="specify fastq file for read 1",
    )

    parser.add_argument(
        "-2",
        "--read2",
        dest="read2",
        required=True,
        type=str,
        help="specify fastq file for read 2",
    )

    parser.add_argument(
        "-w",
        "--whitelist",
        dest="whitelist",
        required=True,
        type=str,
        help="specify a whitelist of accepted cell barcodes",
    )

    parser.add_argument(
        "-f",
        "--feature_ref",
        dest="feature_ref",
        required=True,
        type=str,
        help="specify a reference of feature barcodes",
    )

    parser.add_argument(
        "-o",
        "--output",
        dest="output",
        required=True,
        help="specify an output file",
    )

    parser.add_argument(
        "-r1_c",
        "--read1_coordinate",
        dest="read1_coordinate",
        required=False,
        default=(0, 16),
        type=coords,
        help=(
            "specify coordinate 'start,end' of read 1 to search. "
            "For example, '0,16': starts at 0, stops at 15. "
            "Nucleotide bases outside the range will be masked "
            "as lowercase in the output. Default (0,16)"
        ),
    )

    parser.add_argument(
        "-cb_m",
        "--cb_mismatches",
        dest="cell_barcode_mismatches",
        required=False,
        type=int,
        default=1,
        help="specify cell barcode mismatching threshold. Default (1)",
    )

    parser.add_argument(
        "-cb_n",
        "--cb_num_n_threshold",
        dest="cb_num_n_threshold",
        required=False,
        type=int,
        default=3,
        help=(
            "specify maximum number of ambiguous nucleotides "
            "allowed for read 1. Default (3)"
        ),
    )

    parser.add_argument(
        "-al",
        "--aligner",
        dest="aligner",
        required=False,
        type=str,
        choices=["bwa", "bowtie2"],
        default="bwa",
        help="specify aligner for read 2. Default (bwa)",
    )

    parser.add_argument(
        "--mapq",
        "--mapping_quality",
        dest="mapq",
        required=False,
        default=10,
        type=int,
        help=(
            "specify minimal mapping quality required for feature mapping. "
            "Default (10)"
        ),
    )

    parser.add_argument(
        "-us",
        "--umi_start",
        dest="umi_pos_start",
        required=False,
        type=int,
        default=16,
        help="specify expected UMI starting postion on read 1. Default (16)",
    )

    parser.add_argument(
        "-ul",
        "--umi_length",
        dest="umi_length",
        required=False,
        type=int,
        default=12,
        help=(
            "specify the length of UMIs on read 1. "
            "Reads with UMI length less than this value will be discarded. "
            "Default (12)"
        ),
    )

    parser.add_argument(
        "-um",
        "--umi_mismatches",
        dest="umi_mismatches",
        required=False,
        type=int,
        default=1,
        help=(
            "specify the maximun edit distance allowed for UMIs "
            "on read 1 for deduplication. Default (1)"
        ),
    )

    parser.add_argument(
        "-ud",
        "--umi_deduplication_method",
        dest="umi_deduplication_method",
        required=False,
        type=str,
        choices=[
            "unique",
            "percentile",
            "cluster",
            "adjacency",
            "directional",
        ],
        default="directional",
        help=(
            "specify UMI deduplication method "
            "(powered by UMI-tools. Smith, T., et al. 2017). "
            "Default (directional)"
        ),
    )

    parser.add_argument(
        "--output_directory",
        dest="output_directory",
        required=False,
        type=str,
        default="barcode_mapping",
        help="specify a temp directory. Default (./barcode_mapping)",
    )

    parser.add_argument(
        "-t",
        "--threads",
        dest="threads",
        required=False,
        type=int,
        default=None,
        help=(
            "specify number of threads to launch. "
            "The default is to use all available"
        ),
    )

    parser.add_argument(
        "--num_n_ref",
        dest="num_n_ref",
        required=False,
        type=int,
        default=0,
        help=(
            "specify the number of Ns to separate sequences "
            "belonging to the same feature. Default (0)"
        ),
    )


# filter
def add_filter_subparser(subparsers):
    parser = subparsers.add_parser(
        "filter",
        help="filter extracted barcodes",
        description=(
            "Filter extracted cell and feature barcodes "
            "(output of `extract` or `qc`). "
            "Additional fragment filter/selection "
            "can be applied through `-cb_seq` and/or `-fb_seq`."
        ),
    )

    parser.add_argument(
        "-i",
        "--input",
        dest="input",
        required=True,
        type=str,
        help="specify an input file. The output of `extract` or `qc`",
    )

    parser.add_argument(
        "-o",
        "--output",
        dest="output",
        required=True,
        type=str,
        help="specify an output file",
    )

    parser.add_argument(
        "-cb_s",
        "--cb_start",
        dest="cell_barcode_pos_start",
        required=False,
        type=int,
        default=0,
        help=(
            "specify expected cell barcode starting postion on read 1. " "Default (0)"
        ),
    )

    parser.add_argument(
        "-cb_m",
        "--cb_mismatches",
        dest="cell_barcode_mismatches",
        required=False,
        type=int,
        default=1,
        help="specify cell barcode mismatching threshold. Default (1)",
    )

    parser.add_argument(
        "-cb_ls",
        "--cb_left_shift",
        dest="cell_barcode_left_shift",
        required=False,
        type=int,
        default=1,
        help=(
            "specify the maximum left shift allowed for cell barcode. " "Default (1)"
        ),
    )

    parser.add_argument(
        "-cb_rs",
        "--cb_right_shift",
        dest="cell_barcode_right_shift",
        required=False,
        type=int,
        default=1,
        help=(
            "specify the maximum right shift allowed for cell barcode. " "Default (1)"
        ),
    )

    parser.add_argument(
        "-cb_seq",
        "--cb_extra_seq",
        dest="cell_barcode_extra_seq",
        required=False,
        type=str,
        default=None,
        help=(
            "specify an extra constant sequence to filter on read 1. " "Default (None)"
        ),
    )

    parser.add_argument(
        "-cb_seq_m",
        "--cb_extra_seq_mismatches",
        dest="cell_barcode_extra_seq_mismatches",
        required=False,
        type=int,
        default=None,
        help=(
            "specify the maximun edit distance allowed "
            "for the extra constant sequence on read 1 for filtering. "
            "Default (None)"
        ),
    )

    parser.add_argument(
        "-fb_s",
        "--fb_start",
        dest="feature_barcode_pos_start",
        required=False,
        type=int,
        default=10,
        help=(
            "specify expected feature barcode starting postion on read 2. "
            "Default (10)"
        ),
    )

    parser.add_argument(
        "-fb_m",
        "--fb_mismatches",
        dest="feature_barcode_mismatches",
        required=False,
        type=int,
        default=1,
        help=("specify feature barcode mismatching threshold. Default (1)"),
    )

    parser.add_argument(
        "-fb_ls",
        "--fb_left_shift",
        dest="feature_barcode_left_shift",
        required=False,
        type=int,
        default=1,
        help=(
            "specify the maximum left shift allowed for feature barcode. " "Default (1)"
        ),
    )

    parser.add_argument(
        "-fb_rs",
        "--fb_right_shift",
        dest="feature_barcode_right_shift",
        required=False,
        type=int,
        default=1,
        help=(
            "specify the maximum right shift allowed for feature barcode. "
            "Default (1)"
        ),
    )

    parser.add_argument(
        "-fb_seq",
        "--fb_extra_seq",
        dest="feature_barcode_extra_seq",
        required=False,
        type=str,
        default=None,
        help=(
            "specify an extra constant sequence to filter on read 2. " "Default (None)"
        ),
    )

    parser.add_argument(
        "-fb_seq_m",
        "--fb_extra_seq_mismatches",
        dest="feature_barcode_extra_seq_mismatches",
        required=False,
        type=int,
        default=None,
        help=(
            "specify the maximun edit distance allowed "
            "for the extra constant sequence on read 2. Default (None)"
        ),
    )


# count
def add_count_subparser(subparsers):
    parser = subparsers.add_parser(
        "count",
        help="count feature barcodes per cell",
        description=(
            "Count UMIs per feature per cell (UMI deduplication), "
            "powered by UMI-tools (Smith, T., et al. 2017). "
            "Take the output of `extract` or `filter` as input."
        ),
    )

    parser.add_argument(
        "-i",
        "--input",
        dest="input",
        action="append",
        required=True,
        type=str,
        help="specify input files. Multiple '-i' flags can be used.",
    )

    parser.add_argument(
        "-o",
        "--output",
        dest="output",
        required=True,
        type=str,
        help="specify an output file",
    )

    parser.add_argument(
        "-us",
        "--umi_start",
        dest="umi_pos_start",
        required=False,
        type=int,
        default=16,
        help=(
            "specify expected UMI starting postion on read 1. "
            "Coordinate is 0-based, half-open. Default (16)"
        ),
    )

    parser.add_argument(
        "-ul",
        "--umi_length",
        dest="umi_length",
        required=False,
        type=int,
        default=12,
        help=(
            "specify the length of UMIs on read 1. "
            "Reads with UMI length shorter than this value will be discarded. "
            "Coordinate is 0-based, half-open. "
            "For example, "
            "'-us 16 -ul 12' means UMI starts at 16 ends at 27. "
            "Default (12)"
        ),
    )

    parser.add_argument(
        "-um",
        "--umi_mismatches",
        dest="umi_mismatches",
        required=False,
        type=int,
        default=1,
        help=(
            "specify the maximun edit distance allowed "
            "for UMIs on read 1 for deduplication. Default (1)"
        ),
    )

    parser.add_argument(
        "-ud",
        "--umi_deduplication_method",
        dest="umi_deduplication_method",
        required=False,
        type=str,
        choices=[
            "unique",
            "percentile",
            "cluster",
            "adjacency",
            "directional",
        ],
        default="directional",
        help=(
            "specify UMI deduplication method "
            "(powered by UMI-tools. Smith, T., et al. 2017). "
            "Default (directional)"
        ),
    )

    parser.add_argument(
        "-cb_rc",
        "--cell_barcode_reverse_complement",
        dest="cb_reverse_complement",
        required=False,
        action="store_true",
        help=(
            "specify to convert cell barcode sequences "
            "into their reverse-complement counterparts in the output."
        ),
    )


# demultiplex
def add_demultiplex_subparser(subparsers):
    parser = subparsers.add_parser(
        "demultiplex",
        help="demultiplex cells based on feature abundance",
        description=(
            "Demultiplex cells based on the abundance of features "
            "(matrix generated by `count` as input)."
        ),
    )

    parser.add_argument(
        "-i",
        "--input",
        dest="input",
        required=True,
        type=str,
        help=("specify an input file (feature count matrix). " "The output of `count`"),
    )

    parser.add_argument(
        "--output_directory",
        dest="output_directory",
        required=False,
        type=str,
        default="demultiplexed",
        help="specify a output directory. Default (./demultiplexed)",
    )

    parser.add_argument(
        "-dm",
        dest="demultiplexing_method",
        required=False,
        type=int,
        default=1,
        choices=[1],
        help=(
            "specify demultiplexing method. '1': Stoeckius et al. 2018. " "Default (1)"
        ),
    )

    parser.add_argument(
        "-nm",
        dest="normalization_method",
        required=False,
        type=str,
        default="clr",
        choices=["clr", "log"],
        help=(
            "specify normalization method. "
            "'clr': centred log-ratio transformation. "
            "'log': log10 transformation. "
            "Default (clr)"
        ),
    )

    parser.add_argument(
        "-q",
        "--quantile",
        dest="quantile",
        required=False,
        type=float,
        default=0.9999,
        help=(
            "specify quantile cutoff for the probability mass function "
            "for demultiplexing method 1. "
            "Default (0.9999)"
        ),
    )

    parser.add_argument(
        "-cm",
        "--clustering_method",
        dest="clustering_method",
        required=False,
        type=str,
        default="kmedoids",
        choices=["kmedoids", "hdbscan"],
        help="specify initial clustering method for demultiplexing method 1. "
        "Default (kmedoids)",
    )

    parser.add_argument(
        "-v",
        "--visualization",
        dest="visualization",
        required=False,
        action="store_true",
        help="specify to visualize demultiplexing result",
    )

    parser.add_argument(
        "-vm",
        "--visualization_method",
        dest="visualization_method",
        required=False,
        type=str,
        choices=["tsne", "umap"],
        default="tsne",
        help=(
            "specify embedding method for "
            "visualization (works if '-v' is given). Default (tsne)"
        ),
    )

    parser.add_argument(
        "-nc",
        "--num_cells",
        dest="num_cells",
        required=False,
        type=int,
        default=200,
        help=(
            "specify minimal number of positive cells "
            "required for a feature to be included for demultiplexing. "
            "Default (200)"
        ),
    )


# qc
def add_qc_subparser(subparsers):
    parser = subparsers.add_parser(
        "qc",
        help="quality control of feature barcoding assay",
        description=(
            "Generate diagnostic information. "
            "If `-1` is omitted, "
            "bulk mode is enabled and only read 2 will be analyzed."
        ),
    )

    parser.add_argument(
        "-1",
        "--read1",
        dest="read1",
        required=False,
        type=str,
        help="specify fastq file for read 1",
    )

    parser.add_argument(
        "-2",
        "--read2",
        dest="read2",
        required=True,
        type=str,
        help=(
            "specify fastq file for read 2. "
            "If only read 2 file is provided, "
            "bulk mode is enabled "
            "(skipping arguments '-w/--whitelist', "
            "'-cb_m/--cb_mismatches', "
            "'-r1_c/--read1_coordinate', "
            "must provide '-r2_c/--read2_coordinate' "
            "and '-fb_m/--fb_mismatches'). "
            "In bulk mode, read 2 will be searched against "
            "reference feature barcodes and "
            "read count for each feature barcode will be summarized"
        ),
    )

    parser.add_argument(
        "-w",
        "--whitelist",
        dest="whitelist",
        required=False,
        type=str,
        help="specify a whitelist of accepted cell barcodes",
    )

    parser.add_argument(
        "-f",
        "--feature_ref",
        dest="feature_ref",
        required=True,
        type=str,
        help="specify a reference of feature barcodes",
    )

    parser.add_argument(
        "-r1_c",
        "--read1_coordinate",
        dest="read1_coordinate",
        required=False,
        default=None,
        type=coords,
        help=(
            "specify coordinate 'start,end' of read 1 to search "
            "(doesn't need to be the exact expected barcode coordinate). "
            "Coordinate is 0-based, half-open. "
            "For example, '0,16': starts at 0, stops at 15. "
            "The default is to use all the nucleotide bases. "
            "Nucleotide bases outside the range will be masked "
            "as lowercase in the output"
        ),
    )

    parser.add_argument(
        "-r2_c",
        "--read2_coordinate",
        dest="read2_coordinate",
        required=False,
        default=None,
        type=coords,
        help="see '-r1_c/--read1_coordinate'",
    )

    parser.add_argument(
        "-cb_m",
        "--cb_mismatches",
        dest="cell_barcode_mismatches",
        required=False,
        type=int,
        default=3,
        help="specify cell barcode mismatching threshold. Default (3)",
    )

    parser.add_argument(
        "-fb_m",
        "--fb_mismatches",
        dest="feature_barcode_mismatches",
        required=False,
        type=int,
        default=3,
        help="specify feature barcode mismatching threshold. Default (3)",
    )

    parser.add_argument(
        "-cb_n",
        "--cb_num_n_threshold",
        dest="cb_num_n_threshold",
        required=False,
        type=int,
        default=float("inf"),
        help=(
            "specify maximum number of ambiguous nucleotides allowed "
            "for read 1. The default is no limit"
        ),
    )

    parser.add_argument(
        "-fb_n",
        "--fb_num_n_threshold",
        dest="fb_num_n_threshold",
        required=False,
        type=int,
        default=float("inf"),
        help=(
            "specify maximum number of ambiguous nucleotides allowed "
            "for read 2. The default is no limit"
        ),
    )

    parser.add_argument(
        "-cb_rc",
        "--cell_barcode_reverse_complement",
        dest="cb_reverse_complement",
        required=False,
        action="store_true",
        help=(
            "specify to convert cell barcode sequences "
            "into their reverse-complement counterparts for processing."
        ),
    )

    parser.add_argument(
        "-t",
        "--threads",
        dest="threads",
        required=False,
        type=int,
        default=None,
        help=(
            "specify number of threads for barcode extraction. "
            "Default is to use all available"
        ),
    )

    parser.add_argument(
        "-n",
        "--num_reads",
        dest="num_reads",
        required=False,
        # type=int,
        default=100_000,
        help=(
            "specify number of reads for analysis. "
            "Set to (None) will analyze all the reads. Default (100,000)"
        ),
    )

    parser.add_argument(
        "--chunk_size",
        dest="chunk_size",
        required=False,
        type=int,
        default=50_000,
        help="specify the chunk size for multiprocessing. Default (50,000)",
    )

    parser.add_argument(
        "--output_directory",
        dest="output_directory",
        required=False,
        type=str,
        default="qc",
        help="specify a output directory. Default (./qc)",
    )


# kallisto
def add_kallisto_subparser(subparsers):
    parser = subparsers.add_parser(
        "kallisto_wrapper",
        help="deploy kallisto/bustools for feature barcoding quantification",
        description=(
            "Deploy kallisto/bustools for feature barcoding quantification "
            "(just a wrapper) (Bray, N.L., et al. 2016)."
        ),
    )

    parser.add_argument(
        "-1",
        "--read1",
        dest="read1",
        required=True,
        type=str,
        help="specify fastq file for read 1",
    )

    parser.add_argument(
        "-2",
        "--read2",
        dest="read2",
        required=True,
        type=str,
        help="specify fastq file for read 2",
    )

    parser.add_argument(
        "-w",
        "--whitelist",
        dest="whitelist",
        required=True,
        type=str,
        help="specify a whitelist of accepted cell barcodes",
    )

    parser.add_argument(
        "-f",
        "--feature_ref",
        dest="feature_ref",
        required=True,
        type=str,
        help="specify a reference of feature barcodes",
    )

    parser.add_argument(
        "--technology",
        dest="technology",
        required=False,
        type=str,
        default="10xv3",
        help="specify feature barcoding technology. The default is 10xv3",
    )

    parser.add_argument(
        "-o",
        "--output",
        dest="output",
        required=True,
        help="specify an output file",
    )

    parser.add_argument(
        "-t",
        "--threads",
        dest="threads",
        required=False,
        type=int,
        default=1,
        help=("specify number of kallisto/bustools threads to launch. " "Default (1)"),
    )

    parser.add_argument(
        "--output_directory",
        dest="output_directory",
        required=False,
        type=str,
        default="barcode_mapping",
        help="specify a temp directory. Default (./kallisto)",
    )
