# test_levenshtein.py

import pytest
from fba.levenshtein import (
    indexkeys,
    bytes2set,
    set2bytes,
    create_index,
    query_index,
    # select_query,
    # format_one_query,
    # match_barcodes_paired_fastss,
    # extract_feature_barcoding_fastss
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
    (['123'], 1, {'13': '123', '23': '123', '123': '123', '12': '123'})
]


@pytest.mark.parametrize('barcodes, num_mismatches, z', test_data)
def test_create_index(barcodes, num_mismatches, z):
    assert create_index(barcodes=barcodes, num_mismatches=num_mismatches) == z


test_data = [
    ('AC', {'AC': 'AC'}, 0, {0: ['AC']})
]


@pytest.mark.parametrize('seq, barcode_index, num_mismatches, z', test_data)
def test_query_index(seq, barcode_index, num_mismatches, z):
    assert query_index(
        seq=seq,
        barcode_index=barcode_index,
        num_mismatches=num_mismatches) == z
