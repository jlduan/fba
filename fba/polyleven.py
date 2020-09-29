# polyleven.py

from Bio.SeqIO.QualityIO import FastqGeneralIterator
from polyleven import levenshtein
from itertools import islice, repeat
from multiprocessing import Pool
from fba.utils import open_by_suffix, get_logger


logger = get_logger(logger_name=__name__)


def match_barcodes_polyleven(read_seq,
                             barcodes,
                             read_coords,
                             num_mismatches=3,
                             num_n_threshold=3):
    """Searches one string for known barcode.

    Parameters
    ----------
    read_seq : str
        A DNA string.
    barcodes : list
        A list of barcodes to compare against.
    read_coords : tuple or list
        The positions of read to compare against barcodes.
    num_mismatches : int, optional
        Maximum levenshtein distance allowd.
    num_n_threshold : int, optional
        Maximum Ns allowd.

    Returns
    -------
    str
        The input DNA string.
    str
        Matched barcode.
    int
        The calculated levenshtein distance between input string
        and matched barcode.
    """

    if read_seq.count('N') <= num_n_threshold:
        x, y = read_coords

        for bc in barcodes:
            dist = levenshtein(read_seq[x:y],
                               bc.split('_')[-1],
                               num_mismatches)

            if dist <= num_mismatches:
                break

        if dist <= num_mismatches:
            return read_seq, bc, dist


def match_barcodes_paired_polyleven(read_seqs,
                                    cell_barcodes,
                                    feature_barcodes,
                                    cb_num_mismatches,
                                    fb_num_mismatches,
                                    read1_coords,
                                    read2_coords,
                                    cb_num_n_threshold=3,
                                    fb_num_n_threshold=3):
    """Searches tow strings for known cell and feature barcodes.

    Parameters
    ----------
    read_seqs : tuple or list
        A pair of DNA strings.
    cell_barcodes : list
        A list of cell barcodes (strings) to compare against.
    feature_barcodes : list
        A list of feature barcodes (strings) to compare against.
    cb_num_mismatches : int, optional
        Maximum levenshtein distance allowd for cell barcode matching.
    fb_num_mismatches : int, optional
        Maximum levenshtein distance allowd for feature barcode matching.
    read1_coords : tuple or list
        The positions of read 1 to compare against cell barcodes.
    read2_coords : tuple or list
        The positions of read 2 to compare against feature barcodes.
    cb_num_n_threshold : int, optional
        Maximum Ns allowd for read 1. Read 1 with more Ns than this
        threshold will be skipped.
    fb_num_n_threshold : int, optional
        Maximum Ns allowd for read 2.

    Returns
    -------
    str
        A string of matching result.
    """

    assert len(read_seqs) == 2
    read1_seq, read2_seq = read_seqs

    cb_matched = match_barcodes_polyleven(
        read_seq=read1_seq,
        barcodes=cell_barcodes,
        read_coords=read1_coords,
        num_mismatches=cb_num_mismatches,
        num_n_threshold=cb_num_n_threshold
    )

    if cb_matched:
        fb_matched = match_barcodes_polyleven(
            read_seq=read2_seq,
            barcodes=feature_barcodes,
            read_coords=read2_coords,
            num_mismatches=fb_num_mismatches,
            num_n_threshold=fb_num_n_threshold
        )

    if cb_matched and fb_matched:
        x1, y1 = read1_coords
        x2, y2 = read2_coords

        read1_seq, cb, cb_dist = cb_matched
        read2_seq, fb, fb_dist = fb_matched

        read1_seq = (read1_seq[:x1].lower()
                     + read1_seq[x1:y1]
                     + read1_seq[y1:].lower())

        read2_seq = (read2_seq[:x2].lower()
                     + read2_seq[x2:y2]
                     + read2_seq[y2:].lower())

        out = [
            read1_seq,
            cb,
            cb_dist,
            read2_seq,
            fb,
            fb_dist
        ]

        return '\t'.join([str(i) for i in out])


def extract_feature_barcoding_polyleven(read1_file,
                                        read2_file,
                                        cb_file,
                                        fb_file,
                                        cb_num_mismatches,
                                        fb_num_mismatches,
                                        read1_coords,
                                        read2_coords,
                                        cb_num_n_threshold=3,
                                        fb_num_n_threshold=3,
                                        num_threads=1,
                                        chunk_size=1000):
    """Extracts feature barcodes."""

    logger.info('Using polyleven method ...')

    logger.info(f'Read 1 coordinates to search: {read1_coords}')
    logger.info(f'Read 2 coordinates to search: {read2_coords}')

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

    with open_by_suffix(file_name=cb_file) as f:
        cell_barcodes = [i.split('-')[0].rstrip() for i in f]

    with open_by_suffix(file_name=fb_file) as f:
        feature_barcodes = [i.rstrip().replace('\t', '_') for i in f]

    read1_iter = FastqGeneralIterator(open_by_suffix(file_name=read1_file))
    read2_iter = FastqGeneralIterator(open_by_suffix(file_name=read2_file))

    logger.info('Matching ...')

    def _get_sequence(read1_iter, read2_iter):
        """Gets sequences."""
        for (_, read1_seq, _), (_, read2_seq, _) in zip(read1_iter,
                                                        read2_iter):
            yield read1_seq, read2_seq

    if num_threads == 1:
        for read1_seq, read2_seq in _get_sequence(read1_iter, read2_iter):

            out = match_barcodes_paired_polyleven(
                read_seqs=(read1_seq, read2_seq),
                cell_barcodes=cell_barcodes,
                feature_barcodes=feature_barcodes,
                cb_num_mismatches=cb_num_mismatches,
                fb_num_mismatches=fb_num_mismatches,
                read1_coords=read1_coords,
                read2_coords=read2_coords,
                cb_num_n_threshold=3,
                fb_num_n_threshold=3)

            if out:
                yield out

    else:
        items = list(islice(_get_sequence(read1_iter, read2_iter), chunk_size))

        with Pool(processes=num_threads) as p:

            while items:
                outs = p.starmap(
                    match_barcodes_paired_polyleven,
                    zip(items,
                        repeat(cell_barcodes),
                        repeat(feature_barcodes),
                        repeat(cb_num_mismatches),
                        repeat(fb_num_mismatches),
                        repeat(read1_coords),
                        repeat(read2_coords),
                        repeat(cb_num_n_threshold),
                        repeat(fb_num_n_threshold),
                        )
                )
                outs = '\n'.join([i for i in outs if i])
                yield outs

                items = list(islice(_get_sequence(read1_iter, read2_iter),
                                    chunk_size))


if __name__ == '__main__':
    pass
