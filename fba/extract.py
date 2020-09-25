# extract.py

import regex
from Bio.SeqIO.QualityIO import FastqGeneralIterator
from itertools import islice, repeat
from multiprocessing import Pool
from .utils import open_by_suffix, get_logger


logger = get_logger(logger_name=__name__)


def compile_regex_ref_barcodes_exact(barcodes):
    """Compiles a regular expression pattern, returning a pattern object.

    A list of barcodes are compiled with regex alternation.

    Parameters
    ----------
    barcodes : list
        A list of barcodes.

    Returns
    -------
    Pattern
        A compiled regex patttern of a list of barcodes.
    """

    ref_barcodes = '|'.join([f'({i})' for i in barcodes])
    ref_barcodes = regex.compile(ref_barcodes, regex.BESTMATCH)

    return ref_barcodes


def compile_regex_ref_barcodes_fuzzy(barcodes, num_mismatches=1):
    """Compiles a list of regular expression patterns.

    Parameters
    ----------
    barcodes : list
        A list of barcodes.
    num_mismatches : int, optional
        The maximum number of errors allowd.

    Returns
    -------
    list
        A list of compiled regex pattterns of barcodes.
    """

    ref_barcodes = [regex.compile(f'({i}){{0<e<={num_mismatches}}}',
                                  regex.BESTMATCH) for i in barcodes]

    return ref_barcodes


def compare_against_ref_barcodes(obs_sequence,
                                 compiled_pattern_exact,
                                 compiled_pattern_fuzzy=None):
    """Searches one string for known barcode.

    Parameters
    ----------
    obs_sequence : str
        A string of sequence.
    compiled_pattern_exact :  Pattern
        A compiled pattern object.
    compiled_pattern_fuzzey : list, optional
        A list of compiled pattern objects.

    Returns
    -------
    list
        A list of pattern searching result with four elements.
    """

    matched = compiled_pattern_exact.search(obs_sequence)

    def format_matching(x, y, z=None):
        """Formats matching result."""

        if z:
            return [x,
                    regex.search(
                        r'\((.*?)\)', z.pattern).group(1),
                    ':'.join([str(i) for i in y.span()]),
                    ':'.join([str(i) for i in y.fuzzy_counts])]

        else:
            return [x,
                    y.string[y.start():y.end()],
                    ':'.join([str(i) for i in y.span()]),
                    ':'.join([str(i) for i in y.fuzzy_counts])]

    if matched:
        matched_info = format_matching(
            x=obs_sequence,
            y=matched
        )

    else:
        if compiled_pattern_fuzzy:
            for i in compiled_pattern_fuzzy:
                matched = i.search(obs_sequence)

                if matched:
                    matched_info = format_matching(
                        x=obs_sequence,
                        y=matched,
                        z=i
                    )
                    break

                else:
                    matched_info = [obs_sequence, 'no_match', 'NA', 'NA']
        else:
            matched_info = [obs_sequence, 'no_match', 'NA', 'NA']

    return matched_info


def match_barcodes(obs_sequence,
                   barcodes_compiled_exact,
                   barcodes_compiled_fuzzy=None,
                   num_n_threshold=0):
    """Searches one string for known barcode."""

    matching_result = []
    n_count = obs_sequence.upper().count('N')

    if n_count > num_n_threshold:
        matching_result.extend((obs_sequence, 'n_skipping', 'NA', 'NA'))

    else:
        matched_info = compare_against_ref_barcodes(
            obs_sequence=obs_sequence,
            compiled_pattern_exact=barcodes_compiled_exact,
            compiled_pattern_fuzzy=barcodes_compiled_fuzzy
        )

        matching_result.extend(matched_info)

    return matching_result


def match_barcodes_paired(read_seqs,
                          cb_compiled_exact,
                          fb_compiled_exact,
                          cb_compiled_fuzzy=None,
                          fb_compiled_fuzzy=None,
                          cb_num_n_threshold=0,
                          fb_num_n_threshold=0):
    """Searches a pair of strings for known cell and feature barcodes."""

    assert len(read_seqs) == 4
    read1_seq, read2_seq, read1_seq2, read2_seq2 = read_seqs

    matching_out = match_barcodes(
        obs_sequence=read1_seq,
        barcodes_compiled_exact=cb_compiled_exact,
        barcodes_compiled_fuzzy=cb_compiled_fuzzy,
        num_n_threshold=cb_num_n_threshold
    )

    if matching_out[1] in ['no_match', 'n_skipping']:
        matching_out.append(read2_seq)
        matching_out.extend(('NA', 'NA', 'NA'))

    else:
        fb_matching = match_barcodes(
            obs_sequence=read2_seq,
            barcodes_compiled_exact=fb_compiled_exact,
            barcodes_compiled_fuzzy=fb_compiled_fuzzy,
            num_n_threshold=fb_num_n_threshold
        )

        matching_out.extend(fb_matching)

    matching_out.extend((read1_seq2, read2_seq2))

    return matching_out


def extract_feature_barcoding(read1_file,
                              read2_file,
                              cb_file,
                              fb_file,
                              cb_num_mismatches,
                              fb_num_mismatches,
                              cb_num_n_threshold=2,
                              fb_num_n_threshold=2,
                              read1_coords=None,
                              read2_coords=None,
                              num_threads=1,
                              chunk_size=1000,
                              num_reads=None):
    """Extracts feature barcodes."""

    logger.info(f'regex version: {regex.__version__}')

    logger.info(
        f'Cell barcode maximum number of mismatches: {cb_num_mismatches}')
    logger.info(
        f'Feature barcode maximum number of mismatches: {fb_num_mismatches}')
    logger.info(
        f'Read 1 maximum number of N allowed: {cb_num_n_threshold}')
    logger.info(
        f'Read 2 maximum number of N allowed: {fb_num_n_threshold}')
    logger.info(f'Number of threads: {num_threads}')

    if num_threads > 1:
        logger.info(f'Chunk size: {chunk_size:,}')
    if num_reads:
        logger.info(f'Number of reads to analyze: {num_reads:,}')
        if chunk_size > num_reads:
            chunk_size = num_reads

    with open_by_suffix(file_name=cb_file) as f:
        cell_barcodes = [i.split('-')[0].rstrip() for i in f]
    logger.info('Loading cell barcodes ...')
    logger.info(f'Number of cell barcodes: {len(cell_barcodes):,}')

    with open_by_suffix(file_name=fb_file) as f:
        feature_barcodes = {
            i.rstrip().split('\t')[1]: i.split('\t')[0] for i in f
        }

    logger.info('Loading feature barcodes ...')
    logger.info(f'Number of feature barcodes: {len(feature_barcodes):,}')

    cell_barcodes_compiled_exact = compile_regex_ref_barcodes_exact(
        barcodes=cell_barcodes
    )
    feature_barcodes_compiled_exact = compile_regex_ref_barcodes_exact(
        barcodes=feature_barcodes.keys()
    )

    if cb_num_mismatches:
        cell_barcodes_compiled_fuzzy = compile_regex_ref_barcodes_fuzzy(
            barcodes=cell_barcodes,
            num_mismatches=cb_num_mismatches
        )
    else:
        cell_barcodes_compiled_fuzzy = None

    if fb_num_mismatches:
        feature_barcodes_compiled_fuzzy = compile_regex_ref_barcodes_fuzzy(
            barcodes=feature_barcodes,
            num_mismatches=fb_num_mismatches
        )
    else:
        feature_barcodes_compiled_fuzzy = None

    read1_iter = FastqGeneralIterator(open_by_suffix(file_name=read1_file))
    read2_iter = FastqGeneralIterator(open_by_suffix(file_name=read2_file))

    def get_sequence(read1_iter, read2_iter,
                     read1_coords=read1_coords, read2_coords=read2_coords):
        """Gets sequences."""

        for (_, read1_seq, _), (_, read2_seq, _) in zip(read1_iter,
                                                        read2_iter):

            if read1_coords:
                r1_start, r1_end = read1_coords
                r1 = read1_seq[r1_start: min(r1_end, len(read1_seq))]
            else:
                r1 = read1_seq

            if read2_coords:
                r2_start, r2_end = read2_coords
                r2 = read2_seq[r2_start: min(r2_end, len(read2_seq))]
            else:
                r2 = read2_seq

            yield r1, r2, read1_seq, read2_seq

    def restore_orig_seq(x,
                         read1_coords=read1_coords,
                         read2_coords=read2_coords):
        """Formats matching output, restores original seqs and coordinates."""

        """
        ['TGATCTTAGAACACGT', 'TGATCTTAGAACACGT', '0:16', '2:0:0', 'GGGGGGGGGGGGGGGGAGGGGGCCGGAAAAGAACCCCGAGAGGCCAGCGCCAAACAAAAAAGAACAAAAAAGAGGAAAAAAAAAAAAAAA', 'no_match', 'NA', 'NA', 'TGATCTTAGAACACGTCAGGGTCCTGAA', 'GGGGGGGGGGGGGGGGAGGGGGCCGGAAAAGAACCCCGAGAGGCCAGCGCCAAACAAAAAAGAACAAAAAAGAGGAAAAAAAAAAAAAAA']
        ['TCTCAGCGTATAGTCC', 'TCTCAGCGTATAGTCC', '0:16', '2:0:0', 'AGCGGGCGCATGTTCCCGCTCAACTATACGAACGGCTTTAAGGCCGGTCCTAGCAACCTGAAGGCTTAGGACTATACGCTGAGACTGTCT', 'TGTTCCCGCTCAACT', '10:25', '0:0:0', 'TCTCAGCGTATAGTCCTAAGCCTTCAGG', 'AGCGGGCGCATGTTCCCGCTCAACTATACGAACGGCTTTAAGGCCGGTCCTAGCAACCTGAAGGCTTAGGACTATACGCTGAGACTGTCT']
        ['CGATCGGGTGTGCGCT', 'no_match', 'NA', 'NA', 'CGATCGGCAGTGCGCTCACCTATTAGCGGCTAAGGCGATCTTGAGAGAGCGCACACCCGATCGCTGTCTCTTATACACATCTGACGCTGC', 'NA', 'NA', 'NA', 'CGATCGGGTGTGCGCTCTCTCAAGATCG', 'CGATCGGCAGTGCGCTCACCTATTAGCGGCTAAGGCGATCTTGAGAGAGCGCACACCCGATCGCTGTCTCTTATACACATCTGACGCTGC']
        ['CAACAGTGTAACTAAG', 'CAACAGTGTAACTAAG', '0:16', '2:0:0', 'GGGCAATGTAGCTGCGCTTTCCATTCGAGGCCGGGATTTAAGGCCGGTCCTAGCAANNCGGCTACCCTCTTAGTTACACTGTNGCTGTCT', 'n_skipping', 'NA', 'NA', 'CAACAGTGTAACTAAGAGGGTAGCCGTA', 'GGGCAATGTAGCTGCGCTTTCCATTCGAGGCCGGGATTTAAGGCCGGTCCTAGCAANNCGGCTACCCTCTTAGTTACACTGTNGCTGTCT']
        """ # noqa

        # read1
        if read1_coords:
            r1_start, r1_end = read1_coords
            if ':' in x[2]:
                x[2] = ':'.join(
                    [str(int(i) + r1_start) for i in x[2].split(':')]
                )
            x[0] = x[-2][:r1_start].lower() + x[0] + x[-2][r1_end:].lower()

        # read2
        if read2_coords:
            r2_start, r2_end = read2_coords
            if ':' in x[6]:
                x[6] = ':'.join(
                    [str(int(i) + r2_start) for i in x[6].split(':')]
                )
            x[4] = x[-1][:r2_start].lower() + x[4] + x[-1][r2_end:].lower()

        if x[5] not in {'no_match', 'n_skipping', 'NA'}:
            x[5] = feature_barcodes[x[5]] + '_' + x[5]
        return '\t'.join(x[:-2])

    logger.info('Matching ...')

    reads_ = islice(
        get_sequence(read1_iter, read2_iter),
        0,
        num_reads
    )

    if num_threads == 1:
        for r1, r2, read1_seq, read2_seq in reads_:
            out = match_barcodes_paired(
                read_seqs=(r1, r2, read1_seq, read2_seq),
                cb_compiled_exact=cell_barcodes_compiled_exact,
                fb_compiled_exact=feature_barcodes_compiled_exact,
                cb_compiled_fuzzy=cell_barcodes_compiled_fuzzy,
                fb_compiled_fuzzy=feature_barcodes_compiled_fuzzy,
                cb_num_n_threshold=cb_num_n_threshold,
                fb_num_n_threshold=fb_num_n_threshold
            )

            out = restore_orig_seq(x=out,
                                   read1_coords=read1_coords,
                                   read2_coords=read2_coords)
            yield out

    else:
        items = list(
            islice(
                reads_,
                chunk_size
            )
        )

        with Pool(processes=num_threads) as p:
            while items:
                outs = p.starmap(
                    match_barcodes_paired,
                    zip(items,
                        repeat(cell_barcodes_compiled_exact),
                        repeat(feature_barcodes_compiled_exact),
                        repeat(cell_barcodes_compiled_fuzzy),
                        repeat(feature_barcodes_compiled_fuzzy),
                        repeat(cb_num_n_threshold),
                        repeat(fb_num_n_threshold),
                        )
                )

                outs = [
                    restore_orig_seq(
                        x=i,
                        read1_coords=read1_coords,
                        read2_coords=read2_coords
                    ) for i in outs
                ]
                yield '\n'.join(outs)

                items = list(islice(reads_,
                                    chunk_size))


if __name__ == '__main__':
    pass