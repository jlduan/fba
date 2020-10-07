# kallisto.py

import numpy as np
import pandas as pd
import scipy.io
from pathlib import Path
from fba.utils import (
    open_by_suffix,
    get_binary_path,
    get_logger,
    run_executable
)


logger = get_logger(logger_name=__name__)


def fb2fa_kallisto(x, fasta_file, t2g_file):
    """Prepares fasta file, t2g file and returns k-mer length.

    Parameters
    ----------
    x : str
        The path and name of feature barcode file.

        The example content of the file:
        CD3     CTCATTGTAACTCCT
        CD4     TGTTCCCGCTCAACT
        CD8a    GCTGCGCTTTCCATT
        CD11b   GACAAGTGATCTGCA
        CD14    TCTCAGACCTCCGTA
        CD15    TCACCAGTACCTAGT
        CD16    AAGTTCACTCTTTGC
        CD19    CTGGGCAATTACTCG
        CD20    TTCTGGGTCCCTAGA
        CD25    TTTGTCCTGTACGCC

    fasta_file: str
        The path and name of generated fasta file. One mismatch at each
        coordinate.

        The example content of the file:
        >CD3_CTCATTGTAACTCCT_0_A
        ATCATTGTAACTCCT
        >CD3_CTCATTGTAACTCCT_0_C
        CTCATTGTAACTCCT
        >CD3_CTCATTGTAACTCCT_0_G
        GTCATTGTAACTCCT
        >CD3_CTCATTGTAACTCCT_0_T
        TTCATTGTAACTCCT
        >CD3_CTCATTGTAACTCCT_1_A
        CACATTGTAACTCCT
        >CD3_CTCATTGTAACTCCT_1_C
        CCCATTGTAACTCCT
        >CD3_CTCATTGTAACTCCT_1_G
        CGCATTGTAACTCCT
        >CD3_CTCATTGTAACTCCT_1_T
        CTCATTGTAACTCCT
        >CD3_CTCATTGTAACTCCT_2_A
        CTAATTGTAACTCCT
        >CD3_CTCATTGTAACTCCT_2_C
        CTCATTGTAACTCCT

    t2g_file: str
        The path and name of generated t2g file.

        The example content of the file:
        CD3_CTCATTGTAACTCCT_0_A CD3_CTCATTGTAACTCCT     CD3_CTCATTGTAACTCCT
        CD3_CTCATTGTAACTCCT_0_C CD3_CTCATTGTAACTCCT     CD3_CTCATTGTAACTCCT
        CD3_CTCATTGTAACTCCT_0_G CD3_CTCATTGTAACTCCT     CD3_CTCATTGTAACTCCT
        CD3_CTCATTGTAACTCCT_0_T CD3_CTCATTGTAACTCCT     CD3_CTCATTGTAACTCCT
        CD3_CTCATTGTAACTCCT_1_A CD3_CTCATTGTAACTCCT     CD3_CTCATTGTAACTCCT
        CD3_CTCATTGTAACTCCT_1_C CD3_CTCATTGTAACTCCT     CD3_CTCATTGTAACTCCT
        CD3_CTCATTGTAACTCCT_1_G CD3_CTCATTGTAACTCCT     CD3_CTCATTGTAACTCCT
        CD3_CTCATTGTAACTCCT_1_T CD3_CTCATTGTAACTCCT     CD3_CTCATTGTAACTCCT
        CD3_CTCATTGTAACTCCT_2_A CD3_CTCATTGTAACTCCT     CD3_CTCATTGTAACTCCT
        CD3_CTCATTGTAACTCCT_2_C CD3_CTCATTGTAACTCCT     CD3_CTCATTGTAACTCCT

    Returns
    -------
    int
        Largest odd number of the minimal length of the feature barcodes.
    """

    sequence_lengths = list()
    sequence_names = list()

    with open_by_suffix(file_name=x, mode='r') as f:

        with open_by_suffix(file_name=fasta_file, mode='w') as fo:

            for line in f:
                i = line.rstrip().split('\t')

                sequence = i[1]
                sequence_lengths.append(len(sequence))

                for ii in range(len(sequence)):

                    for iii in 'ACGT':
                        sequence_name = '_'.join(
                            [i[0], sequence, str(ii), iii])
                        sequence_mutated = sequence[0:ii] + \
                            iii + \
                            sequence[ii + 1:len(sequence)]

                        sequence_names.append(
                            '\t'.join(
                                [sequence_name] + [i[0] + '_' + sequence] * 2
                            )
                        )

                        fo.write(
                            '>' + sequence_name
                            + '\n'
                            + sequence_mutated
                            + '\n'
                        )

    with open_by_suffix(file_name=t2g_file, mode='w') as foo:
        foo.write('\n'.join(sequence_names))

    num_sequnces = set([i.split('\t')[1] for i in sequence_names])
    logger.info(f'Number of feature barcodes: {len(num_sequnces)}')

    kmer = min(sequence_lengths)
    if kmer % 2 == 0:
        kmer -= 1
    logger.info(f'k-mer length: {kmer}')

    return kmer


def build_kallisto_index(kallisto_index,
                         kmer,
                         fasta_file):
    """Builds kallisto index.

    A wrapper of `kallisto index [arguments] FASTA-files`.

    Parameters
    ----------
    kallisto_index : str
        The path and name of kallisto index.
    kmer : int
        k-mer length, odd number.
    fasta_file : str
        The path and name of generated fasta file.

    Returns
    -------
    str
        The path and name of generated kallisto index.
    """

    cmd = [
        get_binary_path(binary_name='kallisto'),
        'index',
        '-i',
        str(kallisto_index),
        '-k',
        str(kmer),
        str(fasta_file)
    ]
    outs, errs = run_executable(cmd_line=cmd)
    logger.info(errs)

    return kallisto_index


def align_reads_kallisto(read1_file,
                         read2_file,
                         kallisto_index,
                         output_directory,
                         technology,
                         num_threads=1):
    """Aligns reads."""

    cmd = [
        get_binary_path(binary_name='kallisto'),
        'bus',
        '-i',
        str(kallisto_index),
        '-o',
        str(output_directory),
        '-x',
        technology,
        '-t',
        str(num_threads),
        '-n',
        read1_file,
        read2_file
    ]

    outs, errs = run_executable(cmd)
    logger.debug(' '.join(cmd))

    return outs, errs


def correct_cell_barcodes(cb_file,
                          output_directory,
                          bus_file,
                          corrected_bus_file):
    """Corrects cell barcodes."""

    CELL_BARCODE_FILE = output_directory / 'barcodes_no_suffix.tsv'

    with open_by_suffix(file_name=str(cb_file), mode='r') as f:

        with open_by_suffix(
                file_name=str(CELL_BARCODE_FILE),
                mode='w') as fo:
            for line in f:
                i = line.rstrip().split('-')[0]
                fo.write(i + '\n')
    logger.info('Number of whitelisted cell barcodes: '
                + f'{len([i for i in open(CELL_BARCODE_FILE)])}')

    cmd = [
        get_binary_path(binary_name='bustools'),
        'correct',
        '-w',
        str(output_directory / 'barcodes_no_suffix.tsv'),
        '-o',
        str(corrected_bus_file),
        str(bus_file)
    ]

    outs, errs = run_executable(cmd)
    logger.info(errs)

    return corrected_bus_file


def sort_bus_file(bus_file,
                  sorted_bus_file,
                  num_threads):
    """Sorts bus file."""

    # logger.info(f'Number of threads: {num_threads}')
    cmd = [
        get_binary_path(binary_name='bustools'),
        'sort',
        '-t',
        str(num_threads),
        '-o',
        str(sorted_bus_file),
        str(bus_file)
    ]

    outs, errs = run_executable(cmd)
    logger.info(errs)

    return sorted_bus_file


def count_bus_file(bus_file,
                   output_directory,
                   output_prefix,
                   t2g_file,
                   ec_file,
                   transcripts_file):
    """Generates matrix."""

    cmd = [
        get_binary_path(binary_name='bustools'),
        'count',
        '-o',
        str(output_prefix),
        '--genecounts',
        '-g',
        str(t2g_file),
        '-e',
        str(ec_file),
        '-t',
        str(transcripts_file),
        str(bus_file)
    ]

    outs, errs = run_executable(cmd)

    return output_prefix


def run_kallisto(read1_file,
                 read2_file,
                 cb_file,
                 fb_file,
                 technology='10xv3',
                 output_directory='kallisto',
                 num_threads=1):
    """Runs kallisto/bustools."""

    output_directory = Path(output_directory)
    output_directory.mkdir(exist_ok=True)

    FASTA_FILE = output_directory / 'feature_barcode_ref.fasta'
    T2G_FILE = output_directory / 'feature_barcode_ref.t2g'
    KALLISTO_INDEX = output_directory / 'feature_barcode_ref.idx'

    BUS_FILE = output_directory / 'output.bus'
    CORRECTED_BUS_FILE = output_directory / 'output_corrected.bus'
    SORTED_BUS_FILE = output_directory / 'output_corrected_sorted.bus'

    MATRIX_PREFIX = output_directory / 'featurecount'
    EC_FILE = output_directory / 'matrix.ec'
    TRANSCRIPTS_FILE = output_directory / 'transcripts.txt'

    logger.info('Preparing feature barcoding fasta file ...')
    kmer = fb2fa_kallisto(x=fb_file,
                          fasta_file=FASTA_FILE,
                          t2g_file=T2G_FILE)

    logger.info('Creating feature barcoding index reference ...')
    build_kallisto_index(kallisto_index=KALLISTO_INDEX,
                         fasta_file=FASTA_FILE,
                         kmer=kmer)

    logger.info(f'Aligning reads, {num_threads} thread(s) ...')
    _, err = align_reads_kallisto(read1_file=read1_file,
                                  read2_file=read2_file,
                                  kallisto_index=KALLISTO_INDEX,
                                  output_directory=output_directory,
                                  technology=technology,
                                  num_threads=num_threads)
    logger.info(f'{err}')

    logger.info('Preparing whitelisted cell barcodes ...')
    logger.info('Correcting cell barcodes based on provided whitelist ...')
    corrected_bus_file = correct_cell_barcodes(
        cb_file=cb_file,
        output_directory=output_directory,
        bus_file=BUS_FILE,
        corrected_bus_file=CORRECTED_BUS_FILE)

    logger.info(f'Sorting bus file, {num_threads} thread(s) ...')
    sorted_bus_file = sort_bus_file(bus_file=corrected_bus_file,
                                    sorted_bus_file=SORTED_BUS_FILE,
                                    num_threads=num_threads)

    logger.info('Counting bus file ...')
    matrix_prefix = count_bus_file(bus_file=sorted_bus_file,
                                   output_directory=output_directory,
                                   output_prefix=str(MATRIX_PREFIX),
                                   t2g_file=str(T2G_FILE),
                                   ec_file=str(EC_FILE),
                                   transcripts_file=str(TRANSCRIPTS_FILE))

    matrix_featurecount = scipy.io.mmread(
        source=str(matrix_prefix) + '.mtx').astype(dtype=np.int_)

    matrix_featurecount = pd.DataFrame.sparse.from_spmatrix(
        data=matrix_featurecount,
        index=[i.rstrip()
               for i in open(file=str(matrix_prefix) + '.barcodes.txt')],
        columns=[i.rstrip()
                 for i in open(file=str(matrix_prefix) + '.genes.txt')]
    ).T

    logger.info('Done.')

    return matrix_featurecount
