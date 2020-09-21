# __main__.py

import sys
import numpy as np
from pathlib import Path
from fba.utils import open_by_suffix, get_logger
from fba.parsers import parse_args
from fba.extract import extract_feature_barcoding
from fba.polyleven import extract_feature_barcoding_polyleven
from fba.map import map_feature_barcoding
from fba.filter import filter_matching
from fba.count import generate_matrix
from fba.demultiplex import demultiplex_feature_barcoding
from fba.qc import summarize_sequence_content, summarize_barcode_positions
from fba.kallisto import run_kallisto


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

    logger.info('Initiating logging ...')
    logger.info(
        f'Python version: {sys.version_info.major}.{sys.version_info.minor}')

    if not sys.version_info.major == 3 and sys.version_info.minor >= 6:
        logger.critical('Please use Python >= 3.6')
        sys.exit(1)

    if (args.command == 'extract'):
        logger.info('Using extract subcommand ...')

        if args.method.lower() == 'regex':
            with open_by_suffix(file_name=args.output, mode='w') as f:

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

                for out in extract_feature_barcoding(
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
                        chunk_size=args.chunk_size):

                    f.write(out + '\n')

        elif args.method == 'polyleven':
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

                for out in extract_feature_barcoding_polyleven(
                        read1_file=args.read1,
                        read2_file=args.read2,
                        cb_file=args.whitelist,
                        fb_file=args.feature_ref,
                        cb_num_mismatches=args.cell_barcode_mismatches,
                        fb_num_mismatches=args.feature_barcode_mismatches,
                        read1_coords=args.read1_coords,
                        read2_coords=args.read2_coords,
                        cb_num_n_threshold=args.cb_num_n_threshold,
                        fb_num_n_threshold=args.fb_num_n_threshold,
                        num_threads=args.threads,
                        chunk_size=args.chunk_size):

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
            num_n_ref=args.num_n_ref,
            num_mismatches=args.cell_barcode_mismatches,
            umi_pos_start=args.umi_pos_start,
            umi_length=args.umi_length,
            umi_deduplication_method=args.umi_deduplication_method,
            umi_deduplication_threshold=args.umi_mismatches,
            mapq=args.mapq,
            output_directory=args.output_directory,
            num_threads=args.threads,
            chunk_size=args.chunk_size
        )

        matrix_featurecount.to_csv(
            path_or_buf=args.output,
            compression='infer'
        )
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
            umi_length=args.umi_length,
            umi_pos_start=args.umi_pos_start,
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
            seed=42
        )
        logger.info('Done.')

    elif (args.command == 'qc'):
        logger.info('Using qc subcommand ...')

        _ = summarize_sequence_content(read1_file=args.read1,
                                       read2_file=args.read2,
                                       num_reads=args.num_reads,
                                       output_directory=args.output_directory)

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

            for out in extract_feature_barcoding(
                    read1_file=args.read1,
                    read2_file=args.read2,
                    cb_file=args.whitelist,
                    fb_file=args.feature_ref,
                    cb_num_mismatches=3,  # args.cell_barcode_mismatches,
                    fb_num_mismatches=3,  # args.feature_barcode_mismatches,
                    cb_num_n_threshold=np.Inf,
                    fb_num_n_threshold=np.Inf,
                    read1_coords=args.read1_coords,
                    read2_coords=args.read2_coords,
                    num_threads=args.threads,
                    chunk_size=args.chunk_size,
                    num_reads=args.num_reads):

                f.write(out + '\n')

        _ = summarize_barcode_positions(matching_file=OUTPUT_FILE,
                                        output_directory=args.output_directory)
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
