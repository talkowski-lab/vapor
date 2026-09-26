"""pysam-based fetches must return exactly what the original samtools command-line code returned.

Run from the repo root:  pytest tests/test_fetch_equivalence.py
Needs samtools on PATH and the test data in test_data/.
"""
import os
import random

import pytest

import vapor_vali.Simple_function as sf

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
REF = os.path.join(REPO, 'test_data', 'Homo_sapiens_assembly38.fasta')
BAM = os.path.join(REPO, 'test_data', 'chr1.bam')
BED = os.path.join(REPO, 'test_data', 'out.bed')
CHR1_LEN = 248956422

pytestmark = pytest.mark.skipif(not (os.path.exists(REF) and os.path.exists(BAM)), reason='test data missing')


def bed_regions(n=150, seed=1):
    rows = [l.split('\t') for l in open(BED)]
    random.seed(seed)
    out = []
    for r in random.sample(rows, min(n, len(rows))):
        s, e = int(r[1]), int(r[2])
        fl = random.choice([20, 60, 100, 250, 500])
        out.append(('chr1', s - fl, s + fl))
        out.append(('chr1', s - fl, min(e, s + 1500) + fl))
    return out


EDGE = [('chr1', 1, 10), ('chr1', 0, 10), ('chr1', -5, 10), ('chr1', 10, 5), ('chr1', 10, 10),
        ('chr1', 1, 0), ('chr1', 0, 0), ('chr1', CHR1_LEN - 20, CHR1_LEN + 50), ('chr1', CHR1_LEN, CHR1_LEN),
        ('chr1', CHR1_LEN + 1, CHR1_LEN + 50), ('chrFOO', 1, 10), ('chr2', 100000, 100500),
        ('chr1', 10000, 10009), ('chr1', 1000000, 1000000)]


@pytest.mark.parametrize('rev', ['FALSE', 'TRUE'])
def test_ref_seq_readin(rev):
    for chrom, s, e in EDGE + bed_regions(60):
        a = sf.ref_seq_readin_samtools(REF, chrom, s, e, rev)
        b = sf.ref_seq_readin(REF, chrom, s, e, rev)
        assert a == b, (chrom, s, e, rev)


def test_ref_seq_case_kept():
    # sequence bytes (including case, should the reference be soft-masked) are returned unchanged
    seq = sf.ref_seq_readin(REF, 'chr1', 10000, 12000)
    assert len(seq) == 2001 and seq == sf.ref_seq_readin_samtools(REF, 'chr1', 10000, 12000)


def test_chop_pacbio_read_by_pos():
    nonempty = 0
    for chrom, s, e in EDGE + bed_regions(150):
        for fl in (s and [min(500, max(1, (e - s) // 2))] or [1]):
            a = sf.chop_pacbio_read_by_pos_samtools(BAM, chrom, s, e, fl)
            b = sf.chop_pacbio_read_by_pos(BAM, chrom, s, e, fl)
            assert a == b, (chrom, s, e, fl)
            nonempty += bool(a)
    assert nonempty > 20  # the comparison must actually exercise reads


def test_chop_short_windows_many_reads():
    # small windows spanned by many short reads: exercises CIGAR parsing and ordering
    random.seed(7)
    for _ in range(80):
        s = random.randint(1000000, 240000000)
        e = s + random.choice([10, 40, 80, 120])
        for fl in (10, 50, 100):
            assert sf.chop_pacbio_read_by_pos_samtools(BAM, 'chr1', s, e, fl) == \
                sf.chop_pacbio_read_by_pos(BAM, 'chr1', s, e, fl), (s, e, fl)
