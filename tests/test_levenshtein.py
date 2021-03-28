# test_levenshtein.py

import pytest
from fba.levenshtein import (
    indexkeys,
    bytes2set,
    set2bytes,
    create_index,
    query_index,
    select_query,
    format_one_query,
    match_barcodes_paired_fastss
)


test_data = [
    ('ACGT', 0, {'ACGT'}),
    ('ACGT', 1, {'ACG', 'ACGT', 'ACT', 'AGT', 'CGT'}),
    ('abc', 2, {'a', 'ab', 'abc', 'ac', 'b', 'bc', 'c'}),
    ('ACGTACGT', 0, {'ACGTACGT'}),
    ('A', 2, {'', 'A'})
]


@pytest.mark.parametrize('word, max_dist, z', test_data)
def test_indexkeys(word, max_dist, z):
    assert indexkeys(word=word, max_dist=max_dist) == z


def test_bytes2set():
    assert bytes2set(b'a\x00b\x00c') == {u'a', u'b', u'c'}


def test_set2bytes():
    assert set2bytes({u'a', u'b', u'c'}) == b'a\x00b\x00c'


test_data = [
    (['AC'], 0, {'AC': 'AC'}),
    (['AC'], 1, {'A': 'AC', 'AC': 'AC', 'C': 'AC'}),
    (['AC', 'AG'], 0, {'AC': 'AC', 'AG': 'AG'}),
    (['A'], 3, {'': 'A', 'A': 'A'}),
    (['123'], 1, {'13': '123', '23': '123', '123': '123', '12': '123'}),
    (['ACGT', 'ACGTA'], 1, {'ACG': 'ACGT',
                            'ACT': 'ACGT',
                            'ACGT': 'ACGT\x00ACGTA',
                            'CGT': 'ACGT',
                            'AGT': 'ACGT',
                            'CGTA': 'ACGTA',
                            'ACGA': 'ACGTA',
                            'ACGTA': 'ACGTA',
                            'AGTA': 'ACGTA',
                            'ACTA': 'ACGTA'}

     )
]


@pytest.mark.parametrize('barcodes, num_mismatches, z', test_data)
def test_create_index(barcodes, num_mismatches, z):
    assert create_index(barcodes=barcodes, num_mismatches=num_mismatches) == z


test_data = [
    ('AC', {'AC': 'AC'}, 0, {0: ['AC']}),
    ('ACGT', {'ACGT': 'ACGT'}, 0, {0: ['ACGT']}),
    ('ACGT', {'ACG': 'ACGT',
              'ACT': 'ACGT',
              'ACGT': 'ACGT\x00ACGTA',
              'CGT': 'ACGT',
              'AGT': 'ACGT',
              'CGTA': 'ACGTA',
              'ACGA': 'ACGTA',
              'ACGTA': 'ACGTA',
              'AGTA': 'ACGTA',
              'ACTA': 'ACGTA'}, 1, {0: ['ACGT'], 1:[]}),
    ('ACGT', {'ACG': 'ACGT',
              'ACT': 'ACGT',
              'ACGT': 'ACGT\x00ACGTA',
              'CGT': 'ACGT',
              'AGT': 'ACGT',
              'CGTA': 'ACGTA',
              'ACGA': 'ACGTA',
              'ACGTA': 'ACGTA',
              'AGTA': 'ACGTA',
              'ACTA': 'ACGTA'}, 0, {0: []})

]


@pytest.mark.parametrize('seq, barcode_index, num_mismatches, z', test_data)
def test_query_index(seq, barcode_index, num_mismatches, z):
    assert query_index(
        seq=seq,
        barcode_index=barcode_index,
        num_mismatches=num_mismatches) == z


test_data = [
    ({0: ['ACGT'], 1:[]}, 'ACGT', '$$$$', ('ACGT', 0)),
    ({0: ['ACGT'], 1:[1111]}, 'ACGT', '$$$$', ('ACGT', 0)),
    ({0: [], 1:['ACTT', 'ACCT']}, 'ACGT', '$$$$', ('ACTT', 1))
]


@pytest.mark.parametrize('x, read_seq, read_qual, z', test_data)
def test_select_query(x, read_seq, read_qual, z):
    assert select_query(x=x, read_seq=read_seq, read_qual=read_qual) == z


test_data = [
    (('ACGT', 0), 'ACGTACGT', (0, 4), None, ('ACGTacgt', 'ACGT', '0')),
    (('ACGT', 1), 'ACGTACGT', (0, 4), None, ('ACGTacgt', 'ACGT', '1')),
    (('ACGT', 1), 'ACGTACGT', (0, 4), {'ACGT': 'a'}, ('ACGTacgt', 'a', '1')),
    (('ACGT', 1), 'ACGTACGT', (1, 5), {'ACGT': 'a'}, ('aCGTAcgt', 'a', '1'))
]


@pytest.mark.parametrize('q, read_seq, read_coords, barcode_dict, z',
                         test_data)
def test_format_one_query(q, read_seq, read_coords, barcode_dict, z):
    assert format_one_query(q=q,
                            read_seq=read_seq,
                            read_coords=read_coords,
                            barcode_dict=barcode_dict) == z


test_data = [
    (
        ('TCACTCGAGGTTACCCTCTTATTACGCC', 'FFFFFFFFFFFFFFFFFFFFFFFFFFFF', 'AAGCAGTGGTATCAACGCAGAGTACATGGGGGCCGGCGAACCAGGAAATAGTTTAAGAGCTAAGCTGGAAACAGCATAGCAAGTTTAAAT',
         'FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF:FFFFFF:FFFFFFFFFFFFFF:F,FFFFF,FFFFFFFFFFFFFFFFFFFF'),
        create_index(['TCACTCGAGGTTGCCC'], 1),
        create_index(['GCCGGCGAACCAGGAAATAG'], 0),
        {'GCCGGCGAACCAGGAAATAG': 'RAB1A-2_GCCGGCGAACCAGGAAATAG'},
        (0, 16),
        (31, 51),
        1,
        0,
        3,
        3,
        ('TCACTCGAGGTTACCCtcttattacgcc', 'TCACTCGAGGTTGCCC', '1',
         'aagcagtggtatcaacgcagagtacatggggGCCGGCGAACCAGGAAATAGtttaagagctaagctggaaacagcatagcaagtttaaat', 'RAB1A-2_GCCGGCGAACCAGGAAATAG', '0')
    )
] # noqa


@pytest.mark.parametrize(', '.join(['read_seqs',
                                    'cb_index',
                                    'fb_index',
                                    'feature_barcodes',
                                    'read1_coords',
                                    'read2_coords',
                                    'cb_num_mismatches',
                                    'fb_num_mismatches',
                                    'cb_num_n_threshold',
                                    'fb_num_n_threshold',
                                    'z']), test_data)
def test_match_barcodes_paired_fastss(read_seqs,
                                      cb_index,
                                      fb_index,
                                      feature_barcodes,
                                      read1_coords,
                                      read2_coords,
                                      cb_num_mismatches,
                                      fb_num_mismatches,
                                      cb_num_n_threshold,
                                      fb_num_n_threshold,
                                      z):
    assert match_barcodes_paired_fastss(
        read_seqs=read_seqs,
        cb_index=cb_index,
        fb_index=fb_index,
        feature_barcodes=feature_barcodes,
        read1_coords=read1_coords,
        read2_coords=read2_coords,
        cb_num_mismatches=cb_num_mismatches,
        fb_num_mismatches=fb_num_mismatches,
        cb_num_n_threshold=cb_num_n_threshold,
        fb_num_n_threshold=fb_num_n_threshold) == z
