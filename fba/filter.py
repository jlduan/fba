# filter.py

import regex
from fba.utils import open_by_suffix, get_logger


logger = get_logger(logger_name=__name__)


def compile_regex_ref_barcodes_single(sequence, num_mismatches=1):
    """Compiles a regular expression pattern, returning a pattern object.

    Parameters
    ----------
    sequence : str
        A DNA string.
    num_mismatches: int, optional
        The fuzziness level, the maximum errors allowed.

    Returns
    -------
    Pattern
        A compiled regex object.
    """

    ref_sequence = regex.compile(
        f'({sequence.upper()}){{e<={num_mismatches}}}',
        regex.BESTMATCH
    )

    return ref_sequence


def is_matched(x,
               barcode_pos_start,
               mismatching_threshold=1,
               left_shift=0,
               right_shift=0,
               sequence_regex=None):
    """Tests whether it is a valid match.

    Parameters
    ----------
    x : list
        A single matching description.
    barcode_pos_start : int
        Expected barcode starting position on read.
    mismatching_threshold : int, optional
        Maximun barcode mismatches allowed.
    left_shift : int, optional
        The deviation allowed for the barcode starting position.
    right_shift : int, optional
        The deviation allowed for the barcode end position. The expected ending
        position is determined based on the expected barcode starting postion
        and its length.
    sequence_regex: regex pattern, optional
        A compiled regex pattern of additional constant sequecne to search on
        read.

    Returns
    -------
    bool
        Whether this matching passes filter or not.

    Notes
    -----
    substitutions
    insertions
    deletions

    CATTTCAAGTGAGTA
    ACATTTCAGGT AGTCAgatgcatcggca
    1:16
    2:0:1

    ACACAGTCATAATGCC
    AGTAC CAGT  TAATGCCgtattttcgacc
    3:16
    0:0:3

    CCAATTTAGCGC CTAC
    CCAATTTAGAGCTCT Ctaatctaaagcg
    0:16
    1:1:1
    """

    assert len(x) == 4
    read_seq, barcode, matching_pos, matching_description = x

    if barcode not in {'no_match', 'n_skipping'}:
        matching_pos = [int(i) for i in matching_pos.split(':')]
        matching_description = [int(i)
                                for i in matching_description.split(':')]

        num_mismatches = len(barcode.split('_')[-1]) - \
            (matching_pos[1] - matching_pos[0]) + \
            sum(matching_description)

        if num_mismatches <= mismatching_threshold:

            if (barcode_pos_start
                    - left_shift) <= matching_pos[0] <= (barcode_pos_start
                                                         + left_shift):

                barcode_pos_end = (barcode_pos_start
                                   + len(barcode.split('_')[-1]))

                if (barcode_pos_end
                    - right_shift) <= matching_pos[1] <= (barcode_pos_end
                                                          + right_shift):

                    if sequence_regex:
                        if sequence_regex.search(read_seq):
                            return True
                    else:
                        return True


def filter_matching(matching_file,
                    filtered_matching_file,
                    cb_pos_start=0,
                    cb_num_mismatches=1,
                    cb_left_shift=1,
                    cb_right_shift=1,
                    cb_extra_seq=None,
                    cb_extra_seq_num_mismatches=None,
                    fb_pos_start=10,
                    fb_num_mismatches=1,
                    fb_left_shift=1,
                    fb_right_shift=1,
                    fb_extra_seq=None,
                    fb_extra_seq_num_mismatches=None):
    """Filters raw cell and feature barcode matching result."""

    with open_by_suffix(file_name=matching_file) as f:
        header_line = next(f).rstrip().split('\t')
        logger.info('Header line: {}'.format(' '.join(header_line)))
        logger.info(
            f'Cell barcode maximum number of mismatches: {cb_num_mismatches}')
        logger.info(('Feature barcode maximum number of mismatches: '
                     + f'{fb_num_mismatches}'))

        with open_by_suffix(file_name=filtered_matching_file, mode='w') as fo:
            fo.write('\t'.join(header_line) + '\n')

            read_counter = [int(), int()]

            if len(header_line) == 6:
                logger.info(
                    'Skipping arguments: '
                    + '"cb_pos_start", "cb_left_shift", "cb_right_shift"')

                logger.info(
                    'Skipping arguments: '
                    + '"fb_pos_start", "fb_left_shift", "fb_right_shift"')

            for line in f:
                read_counter[1] += 1

                i = line.rstrip().split('\t')
                if cb_extra_seq and cb_extra_seq_num_mismatches:
                    cell_barcode_sequence_regex = \
                        compile_regex_ref_barcodes_single(
                            cb_extra_seq,
                            num_mismatches=cb_extra_seq_num_mismatches
                        )
                else:
                    cell_barcode_sequence_regex = None

                if len(header_line) == 8 or len(header_line) == 12:
                    cell_barcode_matching = i[:4]
                    cell_barcode_passed = is_matched(
                        x=cell_barcode_matching,
                        barcode_pos_start=cb_pos_start,
                        mismatching_threshold=cb_num_mismatches,
                        left_shift=cb_left_shift,
                        right_shift=cb_right_shift,
                        sequence_regex=cell_barcode_sequence_regex
                    )
                elif len(header_line) == 6:
                    if int(i[2]) <= cb_num_mismatches:
                        cell_barcode_passed = True
                    else:
                        cell_barcode_passed = False

                if cell_barcode_passed:

                    if fb_extra_seq and fb_extra_seq_num_mismatches:
                        feature_barcode_sequence_regex = \
                            compile_regex_ref_barcodes_single(
                                fb_extra_seq,
                                num_mismatches=fb_extra_seq_num_mismatches
                            )
                    else:
                        feature_barcode_sequence_regex = None

                    if len(header_line) == 8 or len(header_line) == 12:
                        feature_barcode_matching = i[4:8]
                        feature_barcode_passed = is_matched(
                            x=feature_barcode_matching,
                            barcode_pos_start=fb_pos_start,
                            mismatching_threshold=fb_num_mismatches,
                            left_shift=fb_left_shift,
                            right_shift=fb_right_shift,
                            sequence_regex=feature_barcode_sequence_regex
                        )
                    elif len(header_line) == 6:
                        if int(i[5]) <= fb_num_mismatches:
                            feature_barcode_passed = True
                        else:
                            feature_barcode_passed = False

                    if feature_barcode_passed:
                        fo.write(line)
                        read_counter[0] += 1

    logger.info(f'Number of lines processed: {read_counter[1]:,}')
    logger.info(f'Number of lines passed filters: {read_counter[0]:,}')

    return filtered_matching_file


if __name__ == '__main__':
    pass
