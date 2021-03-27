# test_levenshtein.py

from fba.levenshtein import (
    indexkeys,
    bytes2set,
    set2bytes,
    create_index,
    query_index,
    select_query,
    format_one_query,
    match_barcodes_paired_fastss,
    extract_feature_barcoding_fastss
)


def test_indexkeys():
    words = {'ACG', 'ACGT', 'ACT', 'AGT', 'CGT'}
    assert indexkeys(word='ACGT', max_dist=1) == words


def test_bytes2set():
    assert bytes2set(b'a\x00b\x00c') == {u'a', u'b', u'c'}


def test_set2bytes():
    assert set2bytes({u'a', u'b', u'c'}) == b'a\x00b\x00c'


def test_create_index():
    a = create_index(barcodes=['ACGT'], num_mismatches=1)
    b = {'CGT': 'ACGT', 'ACT': 'ACGT',
         'ACGT': 'ACGT', 'AGT': 'ACGT', 'ACG': 'ACGT'}
    assert a == b
