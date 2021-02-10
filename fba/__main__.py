# __main__.py

import sys
import pandas as pd
from pathlib import Path
from fba import __version__
from fba.levenshtein import extract_feature_barcoding_fastss
from fba.map import map_feature_barcoding
from fba.filter import filter_matching
from fba.count import generate_matrix
from fba.qc import (
    summarize_sequence_content,
    summarize_barcode_positions,
    analyze_bulk
)
from fba.regex import extract_feature_barcoding_regex
from fba.parsers import parse_args
from fba.demultiplex import demultiplex_feature_barcoding
from fba.kallisto import run_kallisto
from fba.utils import open_by_suffix, get_logger


def main():
    args = parse_args()

    logger = get_logger(logger_name=__name__)
    banner = """


    █████▒▄▄▄▄    ▄▄▄
    ▓██   ▒▓█████▄ ▒████▄
    ▒████ ░▒██▒ ▄██▒██  ▀█▄
    ░▓█▒  ░▒██░█▀  ░██▄▄▄▄██
    ░▒█░   ░▓█  ▀█▓ ▓█   ▓██▒
    ▒ ░   ░▒▓███▀▒ ▒▒   ▓▒█░
    ░     ▒░▒   ░   ▒   ▒▒ ░
    ░ ░    ░    ░   ░   ▒
            ░            ░  ░
                ░
    """
    logger.info(banner)
    # print(banner)

    logger.info(f'fba version: {__version__}')
    logger.info('Initiating logging ...')
    logger.info(
        f'Python version: {sys.version_info.major}.{sys.version_info.minor}')

    if not sys.version_info.major == 3 and sys.version_info.minor >= 6:
        logger.critical('Please use Python >= 3.6')
        sys.exit(1)

    if (args.command == 'extract'):
        logger.info('Using extract subcommand ...')

        with open_by_suffix(file_name=args.output, mode='w') as f:

            f.write('\t'.join(
                [
                    'read1_seq',
                    'cell_barcode',
                    'cb_num_mismatches',
                    'read2_seq',
                    'feature_barcode',
                    'fb_num_mismatches'
                ]
            ) + '\n')

            for out in extract_feature_barcoding_fastss(
                    read1_file=args.read1,
                    read2_file=args.read2,
                    cb_file=args.whitelist,
                    fb_file=args.feature_ref,
                    cb_num_mismatches=args.cell_barcode_mismatches,
                    fb_num_mismatches=args.feature_barcode_mismatches,
                    read1_coords=args.read1_coords,
                    read2_coords=args.read2_coords,
                    output_file=args.output,
                    cb_num_n_threshold=args.cb_num_n_threshold,
                    fb_num_n_threshold=args.fb_num_n_threshold
            ):
                f.write(out + '\n')

        logger.info('Done.')

    elif (args.command == 'map'):
        logger.info('Using map subcommand ...')

        matrix_featurecount = map_feature_barcoding(
            read1_file=args.read1,
            read2_file=args.read2,
            cb_file=args.whitelist,
            fb_file=args.feature_ref,
            read1_coords=args.read1_coords,
            num_mismatches=args.cell_barcode_mismatches,
            num_n_threshold=args.cb_num_n_threshold,
            num_n_ref=args.num_n_ref,
            umi_pos_start=args.umi_pos_start,
            umi_length=args.umi_length,
            umi_deduplication_method=args.umi_deduplication_method,
            umi_deduplication_threshold=args.umi_mismatches,
            mapq=args.mapq,
            output_directory=args.output_directory,
            num_threads=args.threads,
            aligner=args.aligner
        )

        matrix_featurecount.to_csv(path_or_buf=args.output,
                                   compression='infer')
        logger.info('Done.')

    elif (args.command == 'filter'):
        logger.info('Using filter subcommand ...')

        _ = filter_matching(
            matching_file=args.input,
            filtered_matching_file=args.output,
            cb_pos_start=args.cell_barcode_pos_start,
            cb_num_mismatches=args.cell_barcode_mismatches,
            cb_left_shift=args.cell_barcode_left_shift,
            cb_right_shift=args.cell_barcode_right_shift,
            cb_extra_seq=args.cell_barcode_extra_seq,
            cb_extra_seq_num_mismatches=args.cell_barcode_extra_seq_mismatches,
            fb_pos_start=args.feature_barcode_pos_start,
            fb_num_mismatches=args.feature_barcode_mismatches,
            fb_left_shift=args.feature_barcode_left_shift,
            fb_right_shift=args.feature_barcode_right_shift,
            fb_extra_seq=args.cell_barcode_extra_seq,
            fb_extra_seq_num_mismatches=args.feature_barcode_extra_seq_mismatches)  # noqa
        logger.info(f'Filtered feature barcoding result: {_}')
        logger.info('Done.')

    elif (args.command == 'count'):
        logger.info('Using count subcommand ...')

        matrix_featurecount = generate_matrix(
            matching_file=args.input,
            umi_pos_start=args.umi_pos_start,
            umi_length=args.umi_length,
            umi_deduplication_method=args.umi_deduplication_method,
            umi_deduplication_threshold=args.umi_mismatches
        )

        matrix_featurecount.to_csv(
            path_or_buf=args.output,
            compression='infer'
        )
        logger.info('Done.')

    elif (args.command == 'demultiplex'):
        logger.info('Using demultiplex subcommand ...')

        _ = demultiplex_feature_barcoding(
            matrix_featurecount_file=args.input,
            output_directory=args.output_directory,
            q=args.quantile,
            initial_clustering_methold=args.clustering_method,
            visualization=args.visualization,
            embeding_method=args.visualization_method,
            seed=42
        )
        logger.info('Done.')

    elif (args.command == 'qc'):
        logger.info('Using qc subcommand ...')

        if not isinstance(args.num_reads, int):
            if args.num_reads.isdigit():
                num_reads = int(args.num_reads)
            elif args.num_reads.upper() == 'NONE':
                num_reads = None
            else:
                sys.exit(1)
        else:
            num_reads = args.num_reads

        if args.read1:
            _ = summarize_sequence_content(
                read1_file=args.read1,
                read2_file=args.read2,
                num_reads=num_reads,
                output_directory=args.output_directory
            )

            OUTPUT_FILE = 'feature_barcoding_output.tsv.gz'
            OUTPUT_FILE = str(Path(args.output_directory) / OUTPUT_FILE)
            with open_by_suffix(file_name=OUTPUT_FILE, mode='w') as f:

                f.write('\t'.join(
                    [
                        'read1_seq',
                        'cell_barcode',
                        'cb_matching_pos',
                        'cb_matching_description',
                        'read2_seq',
                        'feature_barcode',
                        'fb_matching_pos',
                        'fb_matching_description'
                    ]
                ) + '\n')

                for out in extract_feature_barcoding_regex(
                        read1_file=args.read1,
                        read2_file=args.read2,
                        cb_file=args.whitelist,
                        fb_file=args.feature_ref,
                        cb_num_mismatches=args.cell_barcode_mismatches,
                        fb_num_mismatches=args.feature_barcode_mismatches,
                        cb_num_n_threshold=args.cb_num_n_threshold,
                        fb_num_n_threshold=args.fb_num_n_threshold,
                        read1_coords=args.read1_coords,
                        read2_coords=args.read2_coords,
                        num_threads=args.threads,
                        chunk_size=args.chunk_size,
                        num_reads=num_reads):

                    f.write(out + '\n')

            _ = summarize_barcode_positions(
                matching_file=OUTPUT_FILE,
                output_directory=args.output_directory)

        else:
            logger.info('Bulk mode enabled: '
                        'only feature barcodes on reads 2 are analyzed')
            if not args.read2_coords:
                logger.critical('Please specify "-r2_coords" in bulk mode')
                sys.exit(1)

            logger.info(
                'Skipping arguments: "-1", "-w", "-cb_m", "-r1_coords"'
            )

            fb_frequency = analyze_bulk(
                read_file=args.read2,
                read_coords=args.read2_coords,
                fb_file=args.feature_ref,
                num_mismatches=args.feature_barcode_mismatches,
                num_n_threshold=args.fb_num_n_threshold,
                num_reads=num_reads
            )

            Path(args.output_directory).mkdir(exist_ok=True)
            OUTPUT_FILE = 'feature_barcode_frequency.csv'
            OUTPUT_FILE = str(Path(args.output_directory) / OUTPUT_FILE)
            logger.info(f'Output file: {OUTPUT_FILE}')

            fb_frequency = pd.DataFrame.from_dict(
                data=fb_frequency,
                orient='index',
                columns=['num_reads']).sort_values(
                by='num_reads',
                ascending=False
            )
            fb_frequency['percentage'] = fb_frequency['num_reads'] / sum(
                fb_frequency['num_reads'])
            fb_frequency.to_csv(path_or_buf=OUTPUT_FILE)
        logger.info('Done.')

    elif (args.command == 'kallisto_wrapper'):

        logger.info('Using kallisto_wrapper subcommand ...')
        matrix_featurecount = run_kallisto(
            read1_file=args.read1,
            read2_file=args.read2,
            cb_file=args.whitelist,
            fb_file=args.feature_ref,
            technology=args.technology,  # '10xv3',
            output_directory=args.output_directory,  # 'kallisto',
            num_threads=args.threads)

        matrix_featurecount.to_csv(
            path_or_buf=args.output,
            compression='infer'
        )


if __name__ == "__main__":
    main()
