# qc.py

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from pathlib import Path
from Bio.SeqIO.QualityIO import FastqGeneralIterator
from itertools import islice, repeat
from multiprocessing import Pool
from fba.utils import open_by_suffix, get_logger
from fba.polyleven import match_barcodes_polyleven


logger = get_logger(logger_name=__name__)


params = {'pdf.fonttype': 42,
          'mathtext.default': 'regular',
          'axes.axisbelow': True}
plt.rcParams.update(params)


def plot_sequence_content(read_composition, title,
                          nucleotide_dict, ax, nucleotides='ACGT'):
    """Plots per base composition.

    Parameters
    ----------
    read_composition : DataFrame
        A DataFraem of per base content. Index is read coordinate. Columns are
        nucleotides in the order of ACGTN.
    title : str
        The title for the generated plot.
    nucleotide_dict : dict
        Color for each base.
    ax : Axes
        Axes for plotting.
    nucleotides : str, optional
        Selected nucleotides to visualize.

    Returns
    -------
    Axes
        Distribution of sequence content per base.
    """

    p_handles = list()
    for i in list(nucleotides):

        p = ax.plot(read_composition.index.values,
                    read_composition[i],
                    c=nucleotide_dict[i],
                    linewidth=1)
        p_handles.append(p[0])

    ax.legend(handles=p_handles,
              labels=list(nucleotides),
              loc='upper left',
              fontsize=6,
              frameon=True,
              shadow=False,
              framealpha=0)

    ax.set_title(label=title, fontsize=7)
    ax.tick_params(labelsize=6, labelcolor='black', direction='out')
    ax.xaxis.set_ticks(range(0, read_composition.shape[0], 2))
    ax.set_yticklabels(labels=[f'{i:3,.1%}' for i in ax.get_yticks()])

    for i in ['top', 'bottom', 'left', 'right']:
        ax.spines[i].set_linewidth(w=0.5)
        ax.spines[i].set_color(c='#333333')

    ax.set_xbound(lower=-1, upper=read_composition.shape[0] + 1)

    return ax


def plot_barcode_startend(s, e, bases, title, ax):
    """Plots barcode starting and ending positions.

    Parameters
    ----------
    s : Series
        The percentage of starting positions on each base. The length equals
        e and bases.
    e : Series
        The percentage of ending positions on each base. The length equals
        s and bases.
    bases : Array
        The bases of reads.
    title : str
        The title for the generated plot.
    ax: Axes
        Axes for plotting.

    Returns
    -------
    Axes
        Distribution of barcode starting and ending positions on reads.
    """

    assert len(s) == len(bases)
    ax.bar(x=bases,
           height=s,
           bottom=0)

    assert len(e) == len(bases)
    ax.bar(x=bases,
           height=e,
           bottom=s)

    ax.set_title(label=title, fontsize=7)
    ax.tick_params(labelsize=6, labelcolor='black', direction='out')
    ax.xaxis.set_ticks(range(0, len(bases), 2))
    ax.set_xbound(lower=-1, upper=len(bases) + 1)
    ax.set_ylim(bottom=0, top=1)
    ax.set_yticklabels(labels=['{:,.1%}'.format(i) for i in ax.get_yticks()])

    for i in ['top', 'bottom', 'left', 'right']:
        ax.spines[i].set_linewidth(w=0.5)
        ax.spines[i].set_color(c='#333333')

    return ax


def summarize_sequence_content(read1_file,
                               read2_file,
                               num_reads=None,
                               output_directory='qc'):
    """Summarizes per base content for reads 1 and reads 2.

    Parameters
    ----------
    read1_file : str
        The path and name of read 1 file.
    read2_file : str
        The path and name of read 2 file.
    num_reads : int, optional
        Number of reads to analyze.
    output_directory : str, optional
        The path and name for the output directory.

    Returns
    -------
    str
        The path and name for the output directory.
    """

    logger.info('Summarizing per base read content ...')
    if num_reads:
        logger.info(f'Number of reads to analyze: {num_reads:,}')
    else:
        logger.info('Number of reads to analyze: all')
    logger.info(f'Output directory: {output_directory}')

    # read1
    Path(output_directory).mkdir(exist_ok=True)
    R1_ACGT_PLOT = \
        Path(output_directory) / 'Pyplot_read1_per_base_seq_content.pdf'
    R1_N_PLOT = \
        Path(output_directory) / 'Pyplot_read1_per_base_seq_content_n.pdf'
    # read2
    R2_ACGT_PLOT = \
        Path(output_directory) / 'Pyplot_read2_per_base_seq_content.pdf'
    R2_N_PLOT = \
        Path(output_directory) / 'Pyplot_read2_per_base_seq_content_n.pdf'

    read1_matrix = []
    read2_matrix = []
    read1_iter = FastqGeneralIterator(open_by_suffix(file_name=read1_file))
    read2_iter = FastqGeneralIterator(open_by_suffix(file_name=read2_file))

    def _get_sequence(read1_iter, read2_iter):
        """Gets sequences."""

        for (_, read1_seq, _), (_, read2_seq, _) in zip(read1_iter,
                                                        read2_iter):
            yield read1_seq, read2_seq

    _reads = islice(
        _get_sequence(read1_iter, read2_iter),
        0,
        num_reads
    )

    counter = 0
    for read1_seq, read2_seq in _reads:
        counter += 1

        read1_matrix.append(read1_seq)
        read2_matrix.append(read2_seq)

    logger.info(f'Number of reads processed: {counter:,}')

    read1_matrix = np.array([list(i) for i in read1_matrix])
    read2_matrix = np.array([list(i) for i in read2_matrix])
    read1_composition = pd.DataFrame(
        {i: (read1_matrix == i).mean(axis=0) for i in list('ACGTN')})
    read2_composition = pd.DataFrame(
        {i: (read2_matrix == i).mean(axis=0) for i in list('ACGTN')})

    nucleotide_dict = {i: j for i, j in zip(
        list('ACGTN'), ['#a0cbe8', '#8cd17d', '#e15759', '#f1ce63', 'black'])}

    read1_length = read1_composition.shape[0]
    read2_length = read2_composition.shape[0]

    # read1
    fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(
        max(2.8,  read1_length / 15), 2.5))
    plot_sequence_content(
        read_composition=read1_composition,
        title='Read 1 per base sequence content',
        nucleotide_dict=nucleotide_dict,
        ax=ax
    )
    plt.tight_layout()
    fig.savefig(fname=R1_ACGT_PLOT,
                transparent=None, bbox_inches='tight')

    fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(
        max(2.8, read1_length / 15), 2.5))
    plot_sequence_content(
        read_composition=read1_composition,
        title='Read 1 per base sequence content',
        nucleotide_dict=nucleotide_dict,
        ax=ax,
        nucleotides='N'
    )
    ax.set_yticklabels(labels=[f'{i:,.3%}' for i in ax.get_yticks()])
    plt.tight_layout()
    fig.savefig(fname=R1_N_PLOT,
                transparent=None, bbox_inches='tight')

    # read2
    fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(
        read2_length / 15, 2.5))
    plot_sequence_content(
        read_composition=read2_composition,
        title='Read 2 per base sequence content',
        nucleotide_dict=nucleotide_dict,
        ax=ax
    )
    plt.tight_layout()
    fig.savefig(fname=R2_ACGT_PLOT,
                transparent=None, bbox_inches='tight')

    fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(
        read2_composition.shape[0] / 15, 2.5))
    plot_sequence_content(
        read_composition=read2_composition,
        title='Read 2 per base sequence content',
        nucleotide_dict=nucleotide_dict,
        ax=ax,
        nucleotides='N'
    )
    ax.set_yticklabels(labels=[f'{i:,.3%}' for i in ax.get_yticks()])
    plt.tight_layout()
    fig.savefig(fname=R2_N_PLOT,
                transparent=None, bbox_inches='tight')

    return output_directory


def summarize_barcode_positions(matching_file, output_directory='qc'):
    """Summarizes barcode positions for reads 1 and reads 2.

    Parameters
    ----------
    matching_file : str
        The path and name of matching result.
    output_directory : str, optional
        The path and name for the output directory.

    Returns
    -------
    str
        The path and name for the output directory.
    """

    logger.info('Summarizing barcode coordinates ...')
    logger.info(f'Output directory: {output_directory}')

    # read1
    Path(output_directory).mkdir(exist_ok=True)
    R1_BC_STARTING_FILE = \
        Path(output_directory) / 'Read1_barcodes_starting.csv'
    R1_BC_ENDING_FILE = \
        Path(output_directory) / 'Read1_barcodes_ending.csv'
    R1_BC_STARTING_ENDING_PLOT = \
        Path(output_directory) / 'Pyplot_read1_barcodes_starting_ending.pdf'
    # read2
    R2_BC_STARTING_FILE = \
        Path(output_directory) / 'Read2_barcodes_starting.csv'
    R2_BC_ENDING_FILE = \
        Path(output_directory) / 'Read2_barcodes_ending.csv'
    R2_BC_STARTING_ENDING_PLOT = \
        Path(output_directory) / 'Pyplot_read2_barcodes_starting_ending.pdf'
    # summary
    CB_MISMATCHES_FILE = \
        Path(output_directory) / 'Read1_barcodes_mismatches.csv'
    FB_MISMATCHES_FILE = \
        Path(output_directory) / 'Read2_barcodes_mismatches.csv'
    MATCHED_BC_RATIO_FILE = Path(
        output_directory) / 'matched_barcode_ratio.csv'

    #
    with open_by_suffix(file_name=matching_file) as f:
        next(f)
        first_line = next(f)

    read1_length = len(first_line.split('\t')[0])
    read2_length = len(first_line.split('\t')[4])

    # barcode starts and ends
    counter = [int(), int()]
    cb_matching_pos = list()
    cb_matching_description = list()
    cb_mismatches = list()
    fb_matching_pos = list()
    fb_matching_description = list()
    fb_mismatches = list()

    with open_by_suffix(file_name=matching_file) as f:
        next(f)
        for line in f:
            i = line.rstrip().split('\t')
            counter[1] += 1

            if (i[2] not in {'no_match', 'n_skipping'}
                    and i[5] not in {'no_match', 'NA'}):
                counter[0] += 1

                cb_matching_pos.append(i[2])
                cb_matching_description.append(i[3])
                _ = [int(ii) for ii in i[2].split(':')]
                cb_mismatches.append(
                    len(i[1]) - (_[1] - _[0]) + sum(
                        [int(ii)
                         for ii in i[3].split(':')])
                )
                fb_matching_pos.append(i[6])
                fb_matching_description.append(i[7])
                _ = [int(ii) for ii in i[6].split(':')]
                fb_mismatches.append(
                    len(i[5]) - (_[1] - _[0]) + sum(
                        [int(ii)
                         for ii in i[7].split(':')])
                )

    counter.append(counter[0] / counter[1])
    with open_by_suffix(file_name=MATCHED_BC_RATIO_FILE, mode='w') as f:
        f.write(
            ','.join(['valid', 'total', 'ratio'])
            + '\n'
            + ','.join([str(i) for i in counter])
            + '\n'
        )

    cb_mismatches = pd.Series(
        cb_mismatches).value_counts().to_frame(name='count')
    cb_mismatches['ratio'] = cb_mismatches['count'] / \
        sum(cb_mismatches['count'])
    cb_mismatches.sort_index().to_csv(CB_MISMATCHES_FILE)

    fb_mismatches = pd.Series(
        fb_mismatches).value_counts().to_frame(name='count')
    fb_mismatches['ratio'] = fb_mismatches['count'] / \
        sum(fb_mismatches['count'])
    fb_mismatches.sort_index().to_csv(FB_MISMATCHES_FILE)

    # cell barcode
    cb_s = [int(i.split(':')[0]) for i in cb_matching_pos]
    cb_e = [int(i.split(':')[1]) - 1 for i in cb_matching_pos]

    cb_start_dist = pd.Series(
        cb_s).value_counts().to_frame(
        name='count').reindex(
            list(range(read1_length))).fillna(0).astype(np.int64)
    cb_start_dist.to_csv(R1_BC_STARTING_FILE)
    cb_end_dist = pd.Series(
        cb_e).value_counts().to_frame(name='count').reindex(
        list(range(read1_length))).fillna(0).astype(np.int64)
    cb_end_dist.to_csv(R1_BC_ENDING_FILE)

    fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(
        max(2.8, read1_length / 15), 2.5))
    plot_barcode_startend(
        s=cb_start_dist['count'] / sum(cb_start_dist['count']),
        e=cb_end_dist['count'] / sum(cb_end_dist['count']),
        bases=cb_start_dist.index.values,
        title='Distribution of cell barcode positions',
        ax=ax
    )
    plt.tight_layout()
    fig.savefig(fname=R1_BC_STARTING_ENDING_PLOT,
                transparent=None, bbox_inches='tight')

    # feature barcode
    fb_s = [int(i.split(':')[0]) for i in fb_matching_pos]
    fb_e = [int(i.split(':')[1]) - 1 for i in fb_matching_pos]

    fb_start_dist = pd.Series(
        fb_s).value_counts().to_frame(
        name='count').reindex(
            list(range(read2_length))).fillna(0).astype(np.int64)
    fb_start_dist.to_csv(R2_BC_STARTING_FILE)
    fb_end_dist = pd.Series(
        fb_e).value_counts().to_frame(name='count').reindex(
        list(range(read2_length))).fillna(0).astype(np.int64)
    fb_end_dist.to_csv(R2_BC_ENDING_FILE)

    fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(
        read2_length / 15, 2.5))
    plot_barcode_startend(
        s=fb_start_dist['count'] / sum(fb_start_dist['count']),
        e=fb_end_dist['count'] / sum(fb_end_dist['count']),
        bases=fb_start_dist.index.values,
        title='Distribution of feature barcode positions',
        ax=ax
    )
    plt.tight_layout()
    fig.savefig(fname=R2_BC_STARTING_ENDING_PLOT,
                transparent=None, bbox_inches='tight')

    return output_directory


def analyze_bulk(read2_file,
                 read2_coords,
                 fb_file,
                 num_mismatches=3,
                 num_n_threshold=3,
                 num_threads=1,
                 chunk_size=10000,
                 num_reads=None):
    """Searches feature barcodes on reads 2 and generates matrix.

    Parameters
    ----------
    read2_file : str
        The path and name of read 2 file.
    read2_coords : tuple or list
        The positions on read 2 to search.
    fb_file : str
        The path and name of feature barcoding file.
    num_mismatches : int, optional
        Maximum levenshtein distance allowed.
    num_n_threshoold : int, optional
        Maximum Ns allowed for reads.
    num_threads : int, optional
        Number of threads to use for searching.
    chunk_size : int, optional
        Chunk size for multiprocessing.
    num_reads ; int, optional
        Number of reads to analyze.

    Returns
    -------
    dict
        Count and frequency of each feature barcode in the provided fastq file.
    """

    logger.info('Using polyleven method ...')
    logger.info(f'Read 2 coordinates to search: {read2_coords}')

    with open_by_suffix(file_name=fb_file) as f:
        feature_barcodes = [i.rstrip().replace('\t', '_') for i in f]
    feature_barcode_dict = {i: int() for i in feature_barcodes}

    logger.info('Number of feature barcodes expected: '
                f'{len(feature_barcode_dict):,}')
    logger.info(
        f'Feature barcode maximum number of mismatches: {num_mismatches}')
    logger.info(
        f'Read 2 maximum number of N allowed: {num_n_threshold}')

    if num_reads:
        logger.info(f'Number of reads to analyze: {num_reads:,}')
    else:
        logger.info('Number of reads to analyze: all')

    logger.info(f'Number of threads: {num_threads}')
    if num_threads > 1:
        logger.info(f'Chunk size: {chunk_size:,}')

    read2_iter = FastqGeneralIterator(open_by_suffix(file_name=read2_file))

    def _get_sequence(read_iter):
        """Gets sequences."""
        for (_, read_seq, _) in read_iter:
            yield read_seq

    _reads = islice(
        _get_sequence(read2_iter),
        0,
        num_reads
    )

    logger.info('Matching ...')

    if num_threads == 1:
        counter = int()

        for read2_seq in _reads:
            counter += 1
            out = match_barcodes_polyleven(read_seq=read2_seq,
                                           barcodes=feature_barcodes,
                                           read_coords=read2_coords,
                                           num_mismatches=num_mismatches,
                                           num_n_threshold=num_n_threshold)
            if out:
                feature_barcode_dict[out[1]] += 1

    else:
        items = list(islice(_reads, chunk_size))

        counter = len(items)
        with Pool(processes=num_threads) as p:

            while items:
                outs = p.starmap(
                    match_barcodes_polyleven,
                    zip(items,
                        repeat(feature_barcodes),
                        repeat(read2_coords),
                        repeat(num_mismatches),
                        repeat(num_n_threshold))
                )

                for i in outs:
                    if i:
                        feature_barcode_dict[i[1]] += 1

                items = list(islice(_reads, chunk_size))
                counter += len(items)

    logger.info(f'Number of reads processed: {counter:,}')
    logger.info(
        'Number of reads with valid feature barcodes: '
        f'{sum(feature_barcode_dict.values()):,}'
    )

    return feature_barcode_dict
