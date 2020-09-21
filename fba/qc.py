# qc.py

from fba.utils import open_by_suffix, get_logger
from Bio.SeqIO.QualityIO import FastqGeneralIterator
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from pathlib import Path


logger = get_logger(logger_name=__name__)


params = {'pdf.fonttype': 42,
          'mathtext.default': 'regular',
          'axes.axisbelow': True
          }
plt.rcParams.update(params)


def plot_sequence_content(read_composition, title,
                          nucleotide_dict, ax, nucleotides='ACGT'):
    """Plots per base composition."""

    p_handles = list()
    for i in list(nucleotides):

        p = ax.plot(
            read_composition.index.values,
            read_composition[i],
            c=nucleotide_dict[i],
            linewidth=1
        )
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
    """Plots barcode starting and ending positions."""

    assert len(s) == len(bases)
    ax.bar(
        x=bases,
        height=s,
        bottom=0
    )

    assert len(e) == len(bases)
    ax.bar(
        x=bases,
        height=e,
        bottom=s
    )

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
                               num_reads=100_000,
                               output_directory='qc'):
    """Summarizes per base content."""

    logger.info('Summarizing per base read content ...')
    if num_reads:
        logger.info(f'Number of reads to analyze: {num_reads:,}')
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

    counter = 0
    for (_, read1_seq, _), (_, read2_seq, _) in zip(read1_iter,
                                                    read2_iter):
        counter += 1

        if counter <= num_reads:
            read1_matrix.append(read1_seq)
            read2_matrix.append(read2_seq)
        else:
            break

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
    """Summarizes barcode positions."""

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
