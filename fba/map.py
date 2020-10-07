# map.py

import sys
import pysam
import dnaio
import numpy as np
import pandas as pd
from pathlib import Path
from collections import defaultdict, Counter
from multiprocessing import cpu_count
from tempfile import _get_candidate_names
from umi_tools import UMIClusterer
from umi_tools import __version__ as umi_tools_version
from packaging import version
from fba.levenshtein import (
    create_index,
    query_index,
    select_query
)
from fba.utils import (
    open_by_suffix,
    get_binary_path,
    run_executable,
    get_logger,
    parse_bowtie2_version,
    parse_samtools_version
)
from fba import __version__


logger = get_logger(logger_name=__name__)


def match_cell_barcodes(reads,
                        barcode_index,
                        read_coords,
                        num_mismatches=1,
                        num_n_threshold=3):
    """Matches cell barcodes.

    Parameters
    ----------
    reads : tuple or list
        Read pair info
        (read_name, read1_seq, read1_qual, read2_seq, read2_qual).
    barcode_index : dict
        A FastSS index of barcodes.
    read_coords : tuple or list
        The positions of read to compare against cell barcodes.
    num_mismatches : int, optional
        Maximum levenshtein distance allowd.
    num_n_threshold : int, optional
        Maximum Ns allowd for read. Read with more Ns than this
        threshold will be skipped.

    Returns
    -------
    read_name : str
        The input read name.
    read1_seq : str
        The input read 1 DNA string.
    read1_qual : str
        The input read 1 quality string.
    read2_seq : str
        The input read 2 DNA string.
    read2_qual : str
        The input read 2 quality string.
    bc : str
        Matched barcode.
    dist : int
        The calculated levenshtein distance between input read
        and matched barcode.
    """

    read_name, read1_seq, read1_qual, read2_seq, read2_qual = reads

    if read1_seq.count('N') <= num_n_threshold:
        x1, y1 = read_coords

        cb_queries = query_index(read1_seq[x1: y1],
                                 barcode_index=barcode_index,
                                 num_mismatches=num_mismatches)

        cb_matched = select_query(cb_queries,
                                  read1_seq[x1: y1],
                                  read1_qual[x1: y1])
        if cb_matched:
            bc, dist = cb_matched

            return read_name, read1_seq, read1_qual, \
                read2_seq, read2_qual, bc, dist


def compose_aln(x):
    """Composes unaligned alignment.

    Parameters
    ----------
    x : tuple or list
        A cell barcode matching result.
        The output of \'match_cell_barcodes\' function.

    Returns
    -------
    AlignedSegment
        Unaligned read 2 with cell barcode matching result as tags.
    """

    read_name, read1_seq, read1_qual, read2_seq, read2_qual, bc, dist = x

    a = pysam.AlignedSegment()
    a.query_name = read_name.split(' ')[0]
    a.flag = 0x4
    a.template_length = len(read2_seq)

    a.query_sequence = read2_seq
    a.query_qualities = pysam.qualitystring_to_array(read2_qual)

    tags = [
        ('RG', 'fba'),
        ('R1', read1_seq),
        ('CB', bc),
        ('CM', dist),
    ]
    a.tags = tags

    return a


def generate_unaligned_bam(read1_file,
                           read2_file,
                           cb_file,
                           fb_file,
                           unaligned_bam_file,
                           read1_coords,
                           num_mismatches=1,
                           num_n_threshold=3,
                           num_n_ref=0):
    """Matches cell barcodes and generates unaligned bam.

    Parameters
    ----------
    read1_file : str
        The path and name of read 1 file.
    read2_file : str
        The path and name of read 2 file.
    cb_file : str
        The path and name of cell barcode file.
    fb_file : str
        The path and name of feature barcode file.
    unaligned_bam_file : str
        The path and name of unaligned file.
    read1_coords : tuple or list
        The positions of read 1 to compare against cell barcodes.
    num_mismatches : int, optional
        Maximum levenshtein distance allowd.
    num_n_threshold : int, optional
        Maximum Ns allowd for read 1. Read 1 with more Ns than this
        threshold will be skipped.
    num_n_ref : int, optional
        Number of Ns to use for separating seqeunces belonging to
        the same feature. Needed for correctly constructing bam header.

    Returns
    -------
    str
        The path and name of unaligned file.
    """

    cell_barcodes = [i.rstrip().split('-')[0]
                     for i in open_by_suffix(cb_file, mode='r')]

    cb_index = create_index(barcodes=cell_barcodes,
                            num_mismatches=num_mismatches)

    # create bam header
    feature_barcodes = dict()
    with open_by_suffix(file_name=fb_file, mode='r') as f:
        for line in f:
            i = line.rstrip().split('\t')

            if i[0] not in feature_barcodes:
                feature_barcodes[i[0]] = []
            feature_barcodes[i[0]].append(i[1])

    feature_barcodes = [
        {'LN': len(('N' * num_n_ref).join(feature_barcodes[i])), 'SN': i}
        for i in feature_barcodes
    ]

    rg = {
        'ID': 'fba',
        'LB': 'null',
        'PL': 'illumina',
        'PU': 'null',
        'SM': 'null'
    }

    pg = {
        'ID': 'fba',
        'PN': 'fba',
        'VN': __version__,
        'CL': ' '.join(sys.argv)
    }

    fb_bam_header = {
        'HD': {'VN': '1.6'},
        'SQ': feature_barcodes,
        'RG': [rg],
        'PG': [pg]
    }

    def _get_sequence(read1_file, read2_file):
        """Gets sequences and qualities."""

        with dnaio.open(file1=read1_file,
                        file2=read2_file,
                        fileformat='fastq',
                        mode='r') as f:
            for rec in f:
                read1, read2 = rec

                yield read1.name, read1.sequence, read1.qualities, \
                    read2.sequence, read2.qualities

    read_counter = [int(), int()]
    with pysam.AlignmentFile(
            unaligned_bam_file, 'wb', header=fb_bam_header) as outf:

        for i in _get_sequence(read1_file, read2_file):
            read_counter[1] += 1

            out = match_cell_barcodes(reads=i,
                                      barcode_index=cb_index,
                                      read_coords=read1_coords,
                                      num_mismatches=num_mismatches,
                                      num_n_threshold=num_n_threshold)
            if out:
                read_counter[0] += 1
                outf.write(compose_aln(out))

    return unaligned_bam_file, read_counter


def fb2fa_concatenated(x, fasta_file, num_n=0):
    """Generates feature barcode fasta file.

    Parameters
    ----------
    x : str
        The path and name of feature barcode file.
    fasta_file : str
        The path and name of generated fasta file.
    num_n : int, optional
        Number of Ns to use for separating seqeunces belonging to
        the same feature.

    Returns
    -------
    str
        The path and name of generated fasta file.
    """

    fb = dict()
    with open_by_suffix(file_name=x, mode='r') as f:
        for line in f:
            i = line.rstrip().split('\t')

            if i[0] not in fb:
                fb[i[0]] = []
            fb[i[0]].append(i[1])

    with open_by_suffix(file_name=fasta_file, mode='w') as fo:
        for i in fb:
            fo.write('>' + i + '\n')
            fo.write(('N' * num_n).join(fb[i]) + '\n')

    return fasta_file


def build_bt2_index(fasta_file,
                    bt2_index_base):
    """Builds bowtie2 index."""

    cmd = [
        get_binary_path(binary_name='bowtie2-build'),
        str(fasta_file),
        str(bt2_index_base)
    ]
    outs, errs = run_executable(cmd_line=cmd)

    return bt2_index_base, errs


def align_reads(unaligned_bam_file,
                bt2_index_base,
                alignment_file,
                temp_prefix,
                num_threads=1):
    """Aligns unaligned bam file."""

    bowtie2_align_parameter = (
        '--preserve-tags -D 20 -R 3 -N 1 -L 15 -i S,1,0.50'
    )

    bowtie2_align_parameter = ' '.join(
        ['-p', str(num_threads), bowtie2_align_parameter]
    )

    cmd = [
        get_binary_path(binary_name='bowtie2'),
        bowtie2_align_parameter,
        ' '.join(['-x', str(bt2_index_base),
                  '-b', str(unaligned_bam_file)]),
        '|',
        get_binary_path(binary_name='samtools'),
        'view -uS - | ',
        get_binary_path(binary_name='samtools'),
        'sort -T',
        temp_prefix,
        '-@',
        str(num_threads),
        '-o',
        str(alignment_file),
        '-'
    ]

    outs, errs = run_executable(cmd_line=cmd)

    return alignment_file, errs


def generate_matrix_from_alignment(alignment_file,
                                   umi_pos_start=16,
                                   umi_length=12,
                                   umi_deduplication_method='directional',
                                   umi_deduplication_threshold=1,
                                   mapq=10):
    """Generates matrix from alignments.

    Parameters
    ----------
    alignment_file : str
        The path and name of alignment file.
    umi_pos_start : int, optional
        The starting coordiate of UMI on read 1. If the input matching result
        is from the regex method of extract subcommand, the staring
        coordinate will be auto determined.
    umi_length : int, optional
        The length of UMI on read 1 after cell barcode. The default is 12.
    umi_deduplication_method : str, optional
        The UMI dedupliation method used in UMI-tools
        (Smith, T., et al. (2017). Genome Res. 27, 491–499.).
        See https://cgatoxford.wordpress.com/2015/08/14/unique-molecular-identifiers-the-problem-the-solution-and-the-proof
    umi_deduplication_threshold : int, optional
        The mismatch tolerance for UMI deduplication.
    mapq : int, optional
        The minimal mapping quality threshod. Alignment with mapq less than
        this value will be discarded.

    Returns
    -------
    DataFrame
        A pandas DataFrame of feature count. The columns are cells and
        the rows are features.
    """  # noqa

    matrix_featurecount = defaultdict(dict)
    with pysam.AlignmentFile(alignment_file, mode='rb') as f:

        references = f.references

        for i in references:
            matrix_featurecount[i] = defaultdict(dict)

            for aln in f.fetch(i):

                if aln.mapping_quality >= mapq:
                    cell_barcode = aln.get_tag('CB')

                    if cell_barcode not in matrix_featurecount[i]:
                        matrix_featurecount[i][cell_barcode] = list()

                    if len(aln.get_tag('R1')) >= umi_pos_start + umi_length:
                        matrix_featurecount[i][cell_barcode].append(
                            aln.get_tag('R1')[
                                umi_pos_start:(umi_pos_start + umi_length)
                            ].encode())

            clusterer = UMIClusterer(cluster_method=umi_deduplication_method)

            for ii in matrix_featurecount[i]:
                umis = matrix_featurecount[i][ii]
                matrix_featurecount[i][ii] = len(
                    clusterer(Counter(umis),
                              threshold=umi_deduplication_threshold)
                )
    # cell_barcodes = sorted(set([j for i in matrix_featurecount
    #                             for j in matrix_featurecount[i].keys()]))

    matrix_featurecount = pd.DataFrame.from_dict(
        matrix_featurecount,
        orient='index').fillna(0).astype(dtype=np.int64)

    return matrix_featurecount


def map_feature_barcoding(read1_file,
                          read2_file,
                          cb_file,
                          fb_file,
                          read1_coords,
                          num_mismatches=1,
                          num_n_threshold=3,
                          num_n_ref=0,
                          umi_pos_start=16,
                          umi_length=12,
                          umi_deduplication_method='directional',
                          umi_deduplication_threshold=1,
                          mapq=10,
                          output_directory='barcode_mapping',
                          num_threads=None):
    """Maps feature barcoding. """

    output_directory = Path(output_directory)
    output_directory.mkdir(exist_ok=True)

    FB_FASTA_FILE = str(output_directory / 'feature_ref.fasta')
    FEATURE_BARCODE_REF = str(output_directory / 'feature_ref')
    FEATURE_BARCODE_INDEX_LOG = str(output_directory / 'bowtie2-build.log')

    UNALIGNED_BAM_FILE = str(output_directory / 'unaligned.bam')
    ALIGNMENT_FILE = str(output_directory / 'aligned.bam')
    ALIGNMENT_LOG = str(output_directory / 'bowtie2.log')

    logger.info(f'bowtie2 version: {parse_bowtie2_version()}')
    logger.info(f'samtools version: {parse_samtools_version()}')

    if version.parse(parse_bowtie2_version()) < version.parse('2.4.0'):
        logger.critical('Please use bowtie2 >= 2.4.0')
        sys.exit(1)

    fasta_file = fb2fa_concatenated(
        x=fb_file, fasta_file=FB_FASTA_FILE, num_n=num_n_ref)
    feature_barcode_ref, _ = build_bt2_index(
        fasta_file=fasta_file,
        bt2_index_base=FEATURE_BARCODE_REF)
    with open_by_suffix(file_name=FEATURE_BARCODE_INDEX_LOG, mode='w') as f:
        f.write(_)

    num_cb = len([i for i in open_by_suffix(cb_file)])
    logger.info(f'Number of reference cell barcodes: {num_cb:,}')
    logger.info(f'Read 1 coordinates to search: {read1_coords}')
    logger.info(
        f'Cell barcode maximum number of mismatches: {num_mismatches}')
    logger.info(
        f'Read 1 maximum number of N allowed: {num_n_threshold}')

    logger.info('Matching cell barcodes (read 1) ...')
    unaligned_bam_file, read_counter = generate_unaligned_bam(
        read1_file=read1_file,
        read2_file=read2_file,
        cb_file=cb_file,
        fb_file=fb_file,
        unaligned_bam_file=UNALIGNED_BAM_FILE,
        read1_coords=read1_coords,
        num_mismatches=num_mismatches,
        num_n_threshold=num_n_threshold,
        num_n_ref=num_n_ref
    )

    logger.info(f'Number of read pairs processed: {read_counter[1]:,}')
    logger.info('Number of read pairs w/ valid cell barcodes: '
                f'{read_counter[0]:,}')

    num_fb = len(set([i.split('\t')[0] for i in open_by_suffix(fb_file)]))
    logger.info(f'Number of reference features: {num_fb:,}')

    if not num_threads:
        num_threads = cpu_count()
    logger.info(f'Number of threads: {num_threads}')

    logger.info('Aligning read 2 ...')
    alignment_file, _ = align_reads(
        unaligned_bam_file=unaligned_bam_file,
        bt2_index_base=feature_barcode_ref,
        alignment_file=ALIGNMENT_FILE,
        temp_prefix=next(_get_candidate_names()),
        num_threads=num_threads)
    pysam.index(alignment_file, alignment_file + '.bai')
    with open_by_suffix(file_name=ALIGNMENT_LOG, mode='w') as f:
        f.write(_)
    logger.info(f'\n{_.rstrip()}')

    logger.info('Generating matrix (UMI deduplication) ...')
    logger.info(f'UMI-tools version: {umi_tools_version}')
    logger.info(f'Mapping quality threshold: {mapq}')

    logger.info(f'UMI starting position on read 1: {umi_pos_start}')
    logger.info(f'UMI length: {umi_length}')
    logger.info('UMI-tools deduplication threshold: '
                f'{umi_deduplication_threshold}')
    logger.info('UMI-tools deduplication method: '
                f'{umi_deduplication_method}')

    matrix_featurecount = generate_matrix_from_alignment(
        alignment_file=alignment_file,
        umi_pos_start=umi_pos_start,
        umi_length=umi_length,
        umi_deduplication_method='directional',
        umi_deduplication_threshold=umi_deduplication_threshold
    )

    logger.info(
        f'Number of cell barcodes detected: {matrix_featurecount.shape[1]:,}')
    logger.info(
        f'Number of features detected: {matrix_featurecount.shape[0]:,}')

    logger.info('Total UMIs after deduplication: '
                f'{matrix_featurecount.values.sum():,}')
    logger.info('Median number of UMIs per cell: '
                f'{np.median(matrix_featurecount.sum(axis=0)):,}')

    return matrix_featurecount
