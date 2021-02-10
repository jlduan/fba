# levenshtein.py

import sys
import dnaio
from itertools import combinations
from polyleven import levenshtein
from fba.utils import open_by_suffix, get_logger


logger = get_logger(logger_name=__name__)


# In-memory implementation of TinyFastSS
# https://github.com/fujimotos/TinyFastSS


def indexkeys(word, max_dist):
    """Return the set of index keys ("variants") of a word.
    >>> indexkeys('aiu', 1)
    {'aiu', 'iu', 'au', 'ai'}
    """

    res = set()
    wordlen = len(word)
    limit = min(max_dist, wordlen) + 1

    for dist in range(limit):
        variants = combinations(word, wordlen - dist)

        for variant in variants:
            res.add(''.join(variant))

    return res


def bytes2set(b, delimiter=b'\x00', encoding='utf-8'):
    """Deserialize bytes into a set of unicode strings.
    >>> int2byte(b'a\x00b\x00c')
    {u'a', u'b', u'c'}
    """

    if not b:
        return set()

    list_ = b.split(delimiter)

    return set(bword.decode(encoding) for bword in list_)


def set2bytes(s, delimiter=b'\x00', encoding='utf-8'):
    """Serialize a set of unicode strings into bytes.
    >>> set2byte({u'a', u'b', u'c')
    b'a\x00b\x00c'
    """

    list_ = []
    for uword in sorted(s):
        bword = uword.encode(encoding)
        list_.append(bword)

    return delimiter.join(list_)


def create_index(barcodes, num_mismatches=1):
    """Creates FastSS index for barcode simirality search.

    Parameters
    ----------
    barcodes : list
        A list of barcodes to compare against.
    num_mismatches : int, optional
        Maximum levenshtein distance allowd.

    Returns
    -------
    dict
        A dictionary of keys to a list of barcodes. Keys are variants of
        original barcodes.
   """

    d = dict()
    for bc in barcodes:

        for key in indexkeys(bc, num_mismatches):
            bkey = key.encode()
            barcode_set = {bc}

            if bkey in d:
                barcode_set |= bytes2set(d[bkey])

            d[bkey] = set2bytes(barcode_set)

    d = {
        sys.intern(i.decode()): sys.intern(d[i].decode())
        for i in d
    }

    return d


def query_index(seq, barcode_index, num_mismatches=1):
    """Performs a fuzzy barcode search, exhaustive.

    Parameters
    ----------
    seq : str
        A DNA string.
    barcode_index : dict
        A FastSS index of barcodes.
    num_mismatches : int, optional
        Maximum levenshtein distance allowd.

    Returns
    -------
    dict
        A dictionary of fuzzy searching result. Keys are levenshtein distance.
        Values are list of matched barcodes.
    """

    res = {d: [] for d in range(num_mismatches + 1)}

    cands = {barcode_index.get(key) for key in indexkeys(seq, num_mismatches)}
    # cands.discard(None)
    # cands = {i for i in cands if i}

    if seq in cands:
        res[0].append(seq)

    else:
        for cand in cands:
            if cand:
                dist = levenshtein(seq, cand, num_mismatches)

                if dist <= num_mismatches:
                    res[dist].append(cand)

    return res


def select_query(x, read_seq, read_qual):
    """Selects query.

    If multiple matches are found, select the one with the smallest levenshtein
    distance. If multiple matches having the same smallest levenshtein
    distance, select the one with lowest sequecing quality sum at the
    mismatched nucleotides. If there's still a tie, select the first one.

    Parameters
    ----------
    x : dict
        A dictionary of fuzzy searching result. Keys are levenshtein distance.
        Values are list of matched barcodes.
    read_seq : str
        A DNA string (not whole read).
    read_qual : str
        A sequencing quality string (not whole read).

    Returns
    -------
    str
        Matched reference barcode.
    int
        Levenshtein distance.
    """

    # for i in sorted(x.keys()):
    for i in x:
        if x[i]:
            if len(x[i]) == 1:
                return x[i][0], i
            else:
                s = [sum([ord(read_qual[idx]) - 33
                          for idx, val in enumerate(read_seq)
                          if ii[idx] != val]) for ii in x[i]]

                s = [idx for idx, val in enumerate(s) if val == min(s)][0]

                return x[i][s], i


def format_one_query(q, read_seq, read_coords, barcode_dict=None):
    """Formats output.

    Parameters
    ----------
    q : tuple
        A dictionary of fuzzy searching result. Keys are levenshtein distance.
        Values are list of matched barcodes.
    read_seq : str
        A DNA string (full length read).
    read_coords : tuple or list
        The positions of read used for comparison.
    barcode_dict : dict, optional
        Names for the matched barcodes. Keys are barcode sequences. Values are
        alternative names.

    Returns
    -------
    str
        Sequencing read (full length).
    str
        Matched barcode or barcode name.
    str
        Levenshtein distance.
    """

    read_seq_m, _ = q
    x, y = read_coords

    read_seq = (read_seq[:x].lower()
                + read_seq[x:y]
                + read_seq[y:].lower())

    if barcode_dict:
        barcode = barcode_dict[q[0]]
    else:
        barcode = q[0]

    return read_seq, barcode, str(q[1])


def match_barcodes_paired_fastss(read_seqs,
                                 cb_index,
                                 fb_index,
                                 feature_barcodes,
                                 read1_coords,
                                 read2_coords,
                                 cb_num_mismatches=1,
                                 fb_num_mismatches=1,
                                 cb_num_n_threshold=3,
                                 fb_num_n_threshold=3):
    """Searches one read pair for known cell and feature barcodes.

    Parameters
    ----------
    read_seqs : tuple
        A pair of DNA strings including sequencing qualities (4 elements).
    cb_index : dict
        FastSS index for cell barcodes.
    fb_index : dict
        FastSS index for feature barcodes.
    feature_barcodes : dict
        A dictionary of feature barcodes. Keys are feature barcode sequences.
        Values are their names.
    read1_coords : tuple
        The positions of read 1 to compare against cell barcodes.
    read2_coords : tuple
        The positions of read 2 to compare against feature barcodes.
    cb_num_mismatches : int, optional
        Maximum levenshtein distance allowd for cell barcode matching.
    fb_num_mismatches : int, optional
        Maximum levenshtein distance allowd for feature barcode matching.
    cb_num_n_threshold : int, optional
        Maximum Ns allowd for read 1. Read 1 with more Ns than this
        threshold will be skipped.
    fb_num_n_threshold : int, optional
        Maximum Ns allowd for read 2.

    Returns
    -------
    list
        A list of matching result.
    """

    read1_seq, read1_qual, read2_seq, read2_qual = read_seqs

    if read1_seq.count('N') <= cb_num_n_threshold:
        x1, y1 = read1_coords

        cb_queries = query_index(read1_seq[x1: y1],
                                 barcode_index=cb_index,
                                 num_mismatches=cb_num_mismatches)

        cb_matched = select_query(cb_queries,
                                  read1_seq[x1: y1],
                                  read1_qual[x1: y1])

        if cb_matched and read2_seq.count('N') <= fb_num_n_threshold:
            x2, y2 = read2_coords

            fb_queries = query_index(read2_seq[x2: y2],
                                     barcode_index=fb_index,
                                     num_mismatches=fb_num_mismatches)

            fb_matched = select_query(fb_queries,
                                      read2_seq[x2: y2],
                                      read2_qual[x2: y2])
            if fb_matched:
                out = format_one_query(cb_matched, read1_seq, read1_coords) \
                    + format_one_query(fb_matched, read2_seq,
                                       read2_coords, feature_barcodes)

                return out


def extract_feature_barcoding_fastss(read1_file,
                                     read2_file,
                                     cb_file,
                                     fb_file,
                                     read1_coords,
                                     read2_coords,
                                     cb_num_mismatches,
                                     fb_num_mismatches,
                                     output_file,
                                     cb_num_n_threshold=3,
                                     fb_num_n_threshold=3):
    """Extracts feature barcodes."""

    with open_by_suffix(file_name=cb_file) as f:
        cell_barcodes = [i.split('-')[0].rstrip() for i in f]

    with open_by_suffix(file_name=fb_file) as f:
        feature_barcodes = {
            i.rstrip().split('\t')[-1]: i.rstrip().replace('\t', '_')
            for i in f
        }

    logger.info(f'Number of reference cell barcodes: {len(cell_barcodes):,}')
    logger.info(
        f'Number of reference feature barcodes: {len(feature_barcodes):,}'
    )

    logger.info('Read 1 coordinates to search: [' +
                ', '.join([str(i) for i in read1_coords]) + ')')
    logger.info('Read 2 coordinates to search: [' +
                ', '.join([str(i) for i in read2_coords]) + ')')

    logger.info(
        f'Cell barcode maximum number of mismatches: {cb_num_mismatches}')
    logger.info(
        f'Feature barcode maximum number of mismatches: {fb_num_mismatches}')
    logger.info(
        f'Read 1 maximum number of N allowed: {cb_num_n_threshold}')
    logger.info(
        f'Read 2 maximum number of N allowed: {fb_num_n_threshold}')

    cb_index = create_index(barcodes=cell_barcodes,
                            num_mismatches=cb_num_mismatches)

    fb_index = create_index(barcodes=feature_barcodes.keys(),
                            num_mismatches=fb_num_mismatches)

    logger.info('Matching ...')

    with dnaio.open(file1=read1_file,
                    file2=read2_file,
                    fileformat='fastq',
                    mode='r') as f:

        read_counter = [int(), int()]
        for rec in f:
            read1, read2 = rec

            read_counter[1] += 1
            if read_counter[1] % 10_000_000 == 0:
                logger.info(f'Read pairs processed: {read_counter[1]:,}')

            out = match_barcodes_paired_fastss(
                read_seqs=(read1.sequence, read1.qualities,
                           read2.sequence, read2.qualities),
                cb_index=cb_index,
                fb_index=fb_index,
                feature_barcodes=feature_barcodes,
                read1_coords=read1_coords,
                read2_coords=read2_coords,
                cb_num_mismatches=cb_num_mismatches,
                fb_num_mismatches=fb_num_mismatches,
                cb_num_n_threshold=cb_num_n_threshold,
                fb_num_n_threshold=fb_num_n_threshold
            )
            if out:
                read_counter[0] += 1
                yield '\t'.join(out)

    logger.info(f'Number of read pairs processed: {read_counter[1]:,}')
    logger.info(f'Number of read pairs w/ valid barcodes: {read_counter[0]:,}')


if __name__ == '__main__':
    pass
