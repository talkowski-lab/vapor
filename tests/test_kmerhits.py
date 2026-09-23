"""The vectorised kmerhits must return exactly the list the original pure-Python code returned."""
import os
import random

import pytest

import vapor_vali.Simple_function as sf

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
REF = os.path.join(REPO, 'test_data', 'Homo_sapiens_assembly38.fasta')
BAM = os.path.join(REPO, 'test_data', 'chr1.bam')


def check(s1, s2, k):
    try:
        expected = sf.kmerhits_python(s1, s2, k, 1, True)
    except KeyError:
        with pytest.raises(KeyError):
            sf.kmerhits(s1, s2, k, 1, True)
        return
    got = sf.kmerhits(s1, s2, k, 1, True)
    assert got == expected, (s1[:50], s2[:50], k)
    assert all(type(a) is int and type(b) is int for a, b in got[:50])


@pytest.mark.parametrize('k', [1, 2, 3, 10, 16, 17, 20, 30, 32, 33, 40])
def test_random_sequences(k):
    rnd = random.Random(k)
    alphabets = ['ACGT', 'ACGTN', 'ACGTacgtNn', 'AT', 'ACGTRYSWKMBDHVryswkmbdhv', 'ACGTX']
    for _ in range(150):
        a1 = rnd.choice(alphabets[:5])
        a2 = rnd.choice(alphabets)
        s1 = ''.join(rnd.choice(a1) for _ in range(rnd.randint(0, 300)))
        s2 = ''.join(rnd.choice(a2) for _ in range(rnd.randint(0, 300)))
        if rnd.random() < 0.3:   # repeats and self comparisons
            s2 = s1
        check(s1, s2, k)


def test_palindromes_and_homopolymers():
    for k in (2, 4, 10, 20):
        check('ACGT' * 50, 'ACGT' * 50, k)
        check('A' * 100 + 'T' * 100, 'T' * 100 + 'A' * 100, k)
        check('N' * 60 + 'ACGTTGCA' * 10, 'N' * 30 + 'acgt' * 10 + 'ACGTTGCA' * 5, k)


def test_invalid_characters_in_seq1_raise_like_original():
    check('ACGTXACGTACGT', 'ACGTACGT', 4)
    check('ACG*', 'ACGT', 3)
    check('ACGU', 'ACGT', 3)
    check('ACGX', 'ACGT', 5)   # shorter than k: no windows, no error


def test_unsupported_configurations_fall_back():
    assert sf.kmerhits('ACGTACGT', 'ACGTACGT', 3, 1, False) == sf.kmerhits_python('ACGTACGT', 'ACGTACGT', 3, 1, False)
    assert sf.kmerhits('ACGTACGT', 'ACGTACGT', 3, 2, True) == sf.kmerhits_python('ACGTACGT', 'ACGTACGT', 3, 2, True)


@pytest.mark.skipif(not os.path.exists(REF), reason='test data missing')
def test_real_reference_and_reads():
    rnd = random.Random(3)
    for _ in range(25):
        s = rnd.randint(1000000, 240000000)
        ref = sf.ref_seq_readin(REF, 'chr1', s, s + rnd.choice([150, 600, 2500]))
        reads = sf.chop_pacbio_read_by_pos(BAM, 'chr1', s + 20, s + 100, 50)
        for k in (10, 20, 30, 40):
            check(ref, ref, k)
            for r in reads[:3]:
                check(r[0], ref[r[1]:], k)
    # a repetitive, N-rich region (chromosome start and a satellite-like region)
    for s, e in ((10000, 13000), (125000000, 125003000)):
        ref = sf.ref_seq_readin(REF, 'chr1', s, e)
        for k in (10, 20, 30, 40):
            check(ref, ref, k)
