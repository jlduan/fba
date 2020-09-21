# map.py

import sys
import pysam
import numpy as np
import pandas as pd
from pathlib import Path
from multiprocessing import Pool
from collections import Counter
from tempfile import _get_candidate_names
from polyleven import levenshtein
from itertools import islice, repeat
from Bio.SeqIO.QualityIO import FastqGeneralIterator
from umi_tools import UMIClusterer
from umi_tools import __version__ as umi_tools_version
from packaging import version
from fba.utils import (
    open_by_suffix,
    get_binary_path,
    run_executable,
    get_logger,
    parse_bowtie2_version,
    parse_samtools_version
)


logger = get_logger(logger_name=__name__)


def match_cell_barcodes_polyleven(reads,
                                  barcodes,
                                  read1_coords,
                                  num_mismatches=3):
    """Matches cell barcodes.

    Parameters
    ----------
    reads : tuple or list
        Read pair info (read_name, read1_seq, read2_seq, read2_qual).
    barcodes : tuple or list
        A list of known cell barcodes.
    read1_coords : tuple or list
        The part of read 1 to compare against cell barcodes.
    num_mismatches : int, optional
        Maximum levenshtein distance allowd.
    """

    read_name, read1_seq, read2_seq, read2_qual = reads
    x, y = read1_coords

    for bc in barcodes:
        dist = levenshtein(
            read1_seq[x:y],
            bc.split('_')[-1],
            num_mismatches)

        if dist <= num_mismatches:
            break

    if dist <= num_mismatches:
        return read_name, read1_seq, read2_seq, read2_qual, bc, dist


def compose_aln(x):
    """Composes alignment."""

    read_name, read1_seq, read2_seq, read2_qual, bc, dist = x

    a = pysam.AlignedSegment()
    a.query_name = read_name.split(' ')[0]
    a.flag = 4
    a.reference_id = 0
    a.reference_start = 0
    a.mapping_quality = 0
    a.next_reference_id = 0
    a.next_reference_start = 0
    a.template_length = 0

    a.query_sequence = read2_seq
    a.cigar = ((5, len(read2_seq)),)
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
                           unaligned_sorted_bam_file,
                           read1_coords,
                           num_mismatches=3,
                           num_n_threshold=3,
                           num_n=100,
                           num_threads=1,
                           chunk_size=1000):
    """Matches cell barcodes and generates unalignmed bam."""

    read1_iter = FastqGeneralIterator(open_by_suffix(file_name=read1_file))
    read2_iter = FastqGeneralIterator(open_by_suffix(file_name=read2_file))

    cell_barcodes = [i.rstrip().split('-')[0]
                     for i in open_by_suffix(cb_file, mode='r')]

    # create bam header
    feature_barcodes = dict()
    with open_by_suffix(file_name=fb_file, mode='r') as f:
        for line in f:
            i = line.rstrip().split('\t')

            if i[0] not in feature_barcodes:
                feature_barcodes[i[0]] = []
            feature_barcodes[i[0]].append(i[1])

    feature_barcodes = [
        {'LN': len(
            ('N' * num_n).join(feature_barcodes[i])), 'SN': i}
        for i in feature_barcodes]

    fb_bam_header = {
        'HD': {'VN': '1.6'},
        'SQ': feature_barcodes}

    with pysam.AlignmentFile(
            unaligned_bam_file, 'wb', header=fb_bam_header) as outf:

        if num_threads == 1:

            for (read_name, read1_seq, _), \
                (_, read2_seq, read2_qual) in zip(read1_iter,
                                                  read2_iter):

                out = match_cell_barcodes_polyleven(
                    reads=(read_name, read1_seq, read2_seq, read2_qual),
                    barcodes=cell_barcodes,
                    read1_coords=read1_coords,
                    num_mismatches=num_mismatches)

                if out:
                    outf.write(compose_aln(out))

        else:

            def get_sequence(read1_iter, read2_iter):
                """Gets sequences, a generator."""

                for (read_name, read1_seq, _), \
                    (_, read2_seq, read2_qual) in zip(read1_iter,
                                                      read2_iter):
                    yield read_name, read1_seq, read2_seq, read2_qual

            items = list(
                islice(
                    get_sequence(read1_iter, read2_iter),
                    chunk_size
                )
            )

            with Pool(processes=num_threads) as p:

                while items:
                    outs = p.starmap(
                        match_cell_barcodes_polyleven,
                        zip(items,
                            repeat(cell_barcodes),
                            repeat(read1_coords),
                            repeat(num_mismatches))
                    )
                    for i in outs:
                        if i:
                            outf.write(compose_aln(i))

                    items = list(islice(get_sequence(read1_iter, read2_iter),
                                        chunk_size))

    pysam.sort('-n', '-o', unaligned_sorted_bam_file, unaligned_bam_file)
    return unaligned_sorted_bam_file


def fb2fa_concatenated(x, fasta_file, num_n=100):
    """Generates feature barcode fasta file."""

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


def align_reads(unaligned_sorted_bam_file,
                bt2_index_base,
                alignment_file,
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
                  '-b', str(unaligned_sorted_bam_file)]),
        '|',
        get_binary_path(binary_name='samtools'),
        'view -uS - | ',
        get_binary_path(binary_name='samtools'),
        'sort -o',
        str(alignment_file),
        '-'
    ]

    outs, errs = run_executable(cmd_line=cmd)

    return alignment_file, errs


def generate_matrix_from_alignment(alignment_file,
                                   umi_pos_start,
                                   umi_length,
                                   umi_deduplication_method,
                                   umi_deduplication_threshold=1,
                                   mapq=10):
    """Generates matrix from alignments."""

    matrix_featurecount = {}
    with pysam.AlignmentFile(alignment_file,
                             mode='rb') as f:

        references = f.references

        for i in references:

            matrix_featurecount[i] = dict()
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
                          num_n_ref=100,
                          num_mismatches=3,
                          umi_pos_start=16,
                          umi_length=12,
                          umi_deduplication_method='directional',
                          umi_deduplication_threshold=1,
                          mapq=10,
                          output_directory='barcode_mapping',
                          num_threads=1,
                          chunk_size=10000):
    """Maps feature barcoding. """

    output_directory = Path(output_directory)
    output_directory.mkdir(exist_ok=True)

    FB_FASTA_FILE = str(output_directory / 'feature_barcode_ref.fasta')
    FEATURE_BARCODE_REF = str(output_directory / 'feature_barcode_ref')
    FEATURE_BARCODE_INDEX_LOG = str(output_directory / 'bowtie2-build_fb.log')

    UNALIGNED_BAM_FILE = str(next(_get_candidate_names()) + '.bam')
    UNALIGNED_SORTED_BAM_FILE = str(output_directory / 'unaligned_sorted.bam')
    ALIGNMENT_FILE = str(output_directory / 'aligned.bam')
    ALIGNMENT_LOG = str(output_directory / 'bowtie2_fb.log')

    logger.info(f'bowtie2 version: {parse_bowtie2_version()}')
    logger.info(f'samtools version: {parse_samtools_version()}')

    if version.parse(parse_bowtie2_version()) < version.parse('2.4.0'):
        logger.critical('Please use bowtie2 >= 2.4.0')
        sys.exit(1)

    logger.info(f'Read 1 coordinates to search: {read1_coords}')
    logger.info(
        f'Cell barcode maximum number of mismatches: {num_mismatches}')
    logger.info(f'Mapping quality threshold: {mapq}')
    logger.info(f'Number of threads: {num_threads}')
    if num_threads > 1:
        logger.info(f'Chunk size: {chunk_size:,}')

    logger.info(f'UMI-tools version: {umi_tools_version}')
    logger.info('UMI-tools deduplication method: '
                + f'{umi_deduplication_method}')
    logger.info('UMI-tools deduplication threshold: '
                + f'{umi_deduplication_threshold}')
    logger.info(f'UMI length: {umi_length}')
    logger.info(f'UMI starting position on read 1: {umi_pos_start}')

    logger.info('Preparing feature barcodes ...')
    fasta_file = fb2fa_concatenated(
        x=fb_file, fasta_file=FB_FASTA_FILE, num_n=num_n_ref)
    feature_barcode_ref, _ = build_bt2_index(
        fasta_file=fasta_file,
        bt2_index_base=FEATURE_BARCODE_REF)
    with open_by_suffix(file_name=FEATURE_BARCODE_INDEX_LOG, mode='w') as f:
        f.write(_)

    logger.info('Matching cell barcodes ...')
    unaligned_sorted_bam_file = generate_unaligned_bam(
        read1_file=read1_file,
        read2_file=read2_file,
        cb_file=cb_file,
        fb_file=fb_file,
        unaligned_bam_file=UNALIGNED_BAM_FILE,
        unaligned_sorted_bam_file=UNALIGNED_SORTED_BAM_FILE,
        read1_coords=read1_coords,
        num_mismatches=num_mismatches,
        num_n=num_n_ref,
        num_threads=num_threads,
        chunk_size=chunk_size
    )
    pysam.index(unaligned_sorted_bam_file, unaligned_sorted_bam_file + '.bai')

    logger.info('Aligning ...')
    alignment_file, _ = align_reads(
        unaligned_sorted_bam_file=unaligned_sorted_bam_file,
        bt2_index_base=feature_barcode_ref,
        alignment_file=ALIGNMENT_FILE,
        num_threads=num_threads)
    pysam.index(alignment_file, alignment_file + '.bai')
    with open_by_suffix(file_name=ALIGNMENT_LOG, mode='w') as f:
        f.write(_)

    logger.info('Generating matrix (UMI deduplication) ...')
    matrix_featurecount = generate_matrix_from_alignment(
        alignment_file=alignment_file,
        umi_pos_start=umi_pos_start,
        umi_length=umi_length,
        umi_deduplication_method='directional',
        umi_deduplication_threshold=umi_deduplication_threshold)

    logger.info(
        f'Number of cell barcodes detected: {matrix_featurecount.shape[1]:,}')
    logger.info(
        f'Number of features detected: {matrix_featurecount.shape[0]:,}')
    logger.info('Total UMIs after deduplication: '
                + f'{matrix_featurecount.values.sum():,}')
    logger.info('Median number of UMIs per cell: '
                + f'{np.median(matrix_featurecount.sum(axis=0)):,}')

    return matrix_featurecount
