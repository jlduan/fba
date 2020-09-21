# count.py

import numpy as np
import pandas as pd
from umi_tools import UMIClusterer
from umi_tools import __version__ as umi_tools_version
from collections import Counter
from .utils import open_by_suffix, get_logger


logger = get_logger(logger_name=__name__)


def generate_matrix(matching_file,
                    umi_length=12,
                    umi_pos_start=None,
                    umi_deduplication_method='directional',
                    umi_deduplication_threshold=1):
    """Generates a matrix based on matching results.

    Parameters
    ----------
    matching_file : string
        The path of the matching result file.
    umi_length : int, optional
        The length of UMI on read 1 after cell barcode. The default is 12.
    umi_pos_start : int, optional
        The starting coordiate of UMI on read 1. If the input matching result
        is from the regex method of extract or map subcommand, the staring
        coordinate will be auto determined.
    umi_pos_start : string, optional
        The UMI dedupliation method used in UMI-tools
        (Smith, T., et al. (2017). Genome Res. 27, 491–499.).
        See https://cgatoxford.wordpress.com/2015/08/14/unique-molecular-identifiers-the-problem-the-solution-and-the-proof
    umi_deduplication_threshold : int, optional
        The mismatch tolerance for UMI deduplication.

    Returns
    -------
    DataFrame
        A pandas DataFrame of feature count table. The columns are cells and
        the rows are features.
    """  # noqa

    logger.info(f'UMI-tools version: {umi_tools_version}')
    logger.info('UMI-tools deduplication method: '
                + f'{umi_deduplication_method}')
    logger.info('UMI-tools deduplication threshold: '
                + f'{umi_deduplication_threshold}')
    logger.info(f'UMI length: {umi_length}')

    matrix_featurecount = {}
    line_counter = int()

    with open_by_suffix(file_name=matching_file) as f:
        header_line = next(f)

        logger.info('Header line: {}'.format(
            header_line.rstrip().replace('\t', ' ')))

        if len(header_line.split('\t')) == 6:
            if umi_pos_start:
                logger.info(
                    f'UMI starting position on read 1: {umi_pos_start}'
                )
            else:
                raise ValueError('Need to specify UMI starting position: -us')
        else:
            logger.info('UMI start position on read 1 auto-detected')

        for line in f:
            i = line.rstrip().split('\t')
            line_counter += 1

            read_seq = i[0]
            cell_barcode = i[1]

            if len(header_line.split('\t')) == 6:
                feature_barcode = i[4]
            else:
                feature_barcode = i[5]
                umi_pos_start = [int(ii) for ii in i[2].split(':')][1]

            umi_pos_end = umi_pos_start + umi_length

            if len(read_seq) >= umi_pos_end:
                umi_seq = read_seq[
                    umi_pos_start:umi_pos_end].upper().encode()

            if cell_barcode not in matrix_featurecount:
                matrix_featurecount[cell_barcode] = {}

            if feature_barcode not in matrix_featurecount[cell_barcode]:
                matrix_featurecount[cell_barcode][feature_barcode] = []

            matrix_featurecount[cell_barcode][
                feature_barcode].append(umi_seq)

    logger.info(f'Number of lines processed: {line_counter:,}')

    cell_barcodes = sorted(matrix_featurecount.keys())
    feature_barcodes = sorted(
        set([ii
             for i in matrix_featurecount
             for ii in matrix_featurecount[i]])
    )
    logger.info(f'Number of cell barcodes detected: {len(cell_barcodes):,}')
    logger.info(f'Number of features detected: {len(feature_barcodes):,}')

    clusterer = UMIClusterer(cluster_method=umi_deduplication_method)
    for i in matrix_featurecount:
        for ii in feature_barcodes:

            umis = matrix_featurecount[i].setdefault(ii, 0)
            if umis:
                matrix_featurecount[i][ii] = len(
                    clusterer(Counter(umis),
                              threshold=umi_deduplication_threshold)
                )

    matrix_featurecount = {i: [matrix_featurecount[i][ii]
                               for ii in feature_barcodes]
                           for i in cell_barcodes}
    matrix_featurecount = pd.DataFrame.from_dict(matrix_featurecount,
                                                 orient='columns')
    matrix_featurecount.index = feature_barcodes

    logger.info('Total UMIs after deduplication: '
                + f'{matrix_featurecount.values.sum():,}')
    logger.info('Median number of UMIs per cell: '
                + f'{np.median(matrix_featurecount.sum(axis=0)):,}')

    return matrix_featurecount


if __name__ == '__main__':
    pass
