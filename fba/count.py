# count.py

from collections import Counter, defaultdict

import numpy as np
import pandas as pd
from umi_tools import UMIClusterer
from umi_tools import __version__ as umi_tools_version

from fba.levenshtein import rev_compl
from fba.utils import get_logger, open_by_suffix

logger = get_logger(logger_name=__name__)


def check_header(matching_file):
    """Checks if the header is correct."""

    headers = list()
    for i in matching_file:
        with open_by_suffix(file_name=i) as f:
            header = next(f)
            headers.append(header)

    headers = [i.strip().split("\t") for i in headers]
    if len(set([len(i) for i in headers])) != 1:
        raise ValueError("The header is not consistent.")
    else:
        return headers[0]


def generate_matrix(
    matching_file,
    umi_pos_start=16,
    umi_length=12,
    umi_deduplication_method="directional",
    umi_deduplication_threshold=1,
    barcode_reverse_complement=False,
):
    """Generates a matrix based on matching results.

    Parameters
    ----------
    matching_file : list
        The paths and names of matching result files.
    umi_pos_start : int, optional
        The starting coordinate of UMI on read 1. If the input matching result
        is from the regex method of extract subcommand, the staring
        coordinate will be auto determined.
    umi_length : int, optional
        The length of UMI on read 1 after cell barcode. The default is 12.
        Set to 0 if no UMIs are present.
    umi_deduplication_method : str, optional
        The UMI deduplication method used in UMI-tools
        (Smith, T., et al. (2017). Genome Res. 27, 491–499.).
        See https://cgatoxford.wordpress.com/2015/08/14/unique-molecular-identifiers-the-problem-the-solution-and-the-proof
    umi_deduplication_threshold : int, optional
        The mismatch tolerance for UMI deduplication.
    barcode_reverse_complement : bool, optional
        Whether to reverse-complement the barcode sequence.
        The default is False.

    Returns
    -------
    DataFrame
        A pandas DataFrame of feature count. The columns are cells and
        the rows are features.
    """  # noqa

    logger.info(f"UMI-tools version: {umi_tools_version}")

    matrix_featurecount = defaultdict(dict)
    line_count = int()

    header = check_header(matching_file=matching_file)

    if umi_length:
        if len(header) == 6:
            if umi_pos_start:
                logger.info(f"UMI starting position on read 1: {umi_pos_start}")
            else:
                logger.critical(
                    "Need to specify UMI starting position on read 1: "
                    '"-us/--umi_start"'
                )
                raise ValueError("Need to specify UMI starting position")
        else:
            logger.info(
                "UMI start position on read 1 auto-detected, "
                'overriding "-us/--umi_start"'
            )

        logger.info(f"UMI length: {umi_length}")

        logger.info(
            "UMI-tools deduplication threshold: "
            f"{umi_deduplication_threshold}"
        )
        logger.info(
            "UMI-tools deduplication method: " f"{umi_deduplication_method}"
        )

    else:
        logger.info(
            "UMI length set to 0, ignoring UMI information. "
            'Skipping arguments: "-us/--umi_start".'
        )

    logger.info("Header: {}".format(" ".join(header)))

    for j in matching_file:
        with open_by_suffix(file_name=j) as f:
            header = next(f)

            for line in f:
                i = line.rstrip().split("\t")
                line_count += 1

                read_seq = i[0]
                cell_barcode = i[1]

                if len(header.split("\t")) == 6:
                    feature_barcode = i[4]
                else:
                    feature_barcode = i[5]
                    umi_pos_start = [int(ii) for ii in i[2].split(":")][1]

                if umi_length:
                    umi_pos_end = umi_pos_start + umi_length
                    if len(read_seq) >= umi_pos_end:
                        umi_seq = (
                            read_seq[umi_pos_start:umi_pos_end].upper().encode()
                        )

                        if (
                            feature_barcode
                            not in matrix_featurecount[cell_barcode]
                        ):
                            matrix_featurecount[cell_barcode][
                                feature_barcode
                            ] = list()

                        matrix_featurecount[cell_barcode][
                            feature_barcode
                        ].append(umi_seq)
                else:
                    if feature_barcode not in matrix_featurecount[cell_barcode]:
                        matrix_featurecount[cell_barcode][
                            feature_barcode
                        ] = int()
                    matrix_featurecount[cell_barcode][feature_barcode] += 1

    logger.info(f"Number of read pairs processed: {line_count:,}")

    cell_barcodes = sorted(matrix_featurecount.keys())
    feature_barcodes = sorted(
        set([ii for i in matrix_featurecount for ii in matrix_featurecount[i]])
    )

    logger.info(f"Number of cell barcodes detected: {len(cell_barcodes):,}")
    logger.info(f"Number of features detected: {len(feature_barcodes):,}")

    if umi_length:
        logger.info("UMI deduplicating ...")

        clusterer = UMIClusterer(cluster_method=umi_deduplication_method)
        for i in matrix_featurecount:
            for ii in feature_barcodes:

                umis = matrix_featurecount[i].setdefault(ii, 0)
                if umis:
                    matrix_featurecount[i][ii] = len(
                        clusterer(
                            Counter(umis), threshold=umi_deduplication_threshold
                        )
                    )

    else:
        logger.info("Counting ...")

        for i in matrix_featurecount:
            for ii in feature_barcodes:

                feature_count = matrix_featurecount[i].setdefault(ii, 0)
                matrix_featurecount[i][ii] = feature_count

    matrix_featurecount = {
        i: [matrix_featurecount[i][ii] for ii in feature_barcodes]
        for i in cell_barcodes
    }
    matrix_featurecount = pd.DataFrame.from_dict(
        matrix_featurecount, orient="columns"
    )
    matrix_featurecount.index = feature_barcodes

    if barcode_reverse_complement:
        matrix_featurecount.columns = [
            rev_compl(i) for i in matrix_featurecount.columns
        ]

    if umi_length:
        logger.info(
            "Total UMIs after deduplication: "
            f"{matrix_featurecount.values.sum():,}"
        )
        logger.info(
            "Median number of UMIs per cell: "
            f"{np.median(matrix_featurecount.sum(axis=0)):,}"
        )
    else:
        logger.info("Total reads: " f"{matrix_featurecount.values.sum():,}")
        logger.info(
            "Median number of reads per cell: "
            f"{np.median(matrix_featurecount.sum(axis=0)):,}"
        )

    return matrix_featurecount


if __name__ == "__main__":
    pass
