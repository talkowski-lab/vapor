"""Per-SV driver for `vapor bed`, optionally spread over several worker processes.

Every SV is scored independently (clustering is seeded, file handles are per process), so the
records can be processed in any process and in any order; results are written back strictly in
input order, and each record's stdout is replayed in that order too, so the output files and the
log are the same as a single-process run.
"""
from __future__ import print_function

import io
import multiprocessing
import os
import sys
from contextlib import redirect_stdout

_THREAD_ENV = ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS',
               'VECLIB_MAXIMUM_THREADS', 'NUMEXPR_NUM_THREADS')


def limit_native_threads():
    """One native thread per process unless the user asked otherwise.

    Must run before numpy/scipy/sklearn are imported. Parallelism comes from processes;
    extra BLAS/OpenMP threads per process would only oversubscribe the CPUs.
    """
    for k in _THREAD_ENV:
        os.environ.setdefault(k, '1')


def _score_bed_record(task):
    """Score one bed record exactly as the original `vapor bed` loop did.

    Returns (row for write_output_main or None, text the original printed to stdout).
    """
    from vapor_vali import Simple_function as sf
    plt_li, x, num_reads_cff, bam_in, ref, out_path, sample_name = task
    buf = io.StringIO()
    row = None
    with redirect_stdout(buf):
        if x[-1] in ['a/', '/a', '/', 'DEL']:
            key_event = ':'.join([str(i) for i in x[:-3]] + ['DEL'])
            vapor_score_event = sf.vapor_simple_del_Vapor(num_reads_cff, plt_li, bam_in, ref, x[:-3], out_path + sample_name + '.DEL.' + key_event.replace(':', '__') + '.png')
        elif x[-1] in ['a/a^', 'a^/a', 'a^/a^', 'INV']:
            key_event = ':'.join([str(i) for i in x[:-3]] + ['INV'])
            vapor_score_event = sf.vapor_simple_inv_Vapor(num_reads_cff, plt_li, bam_in, ref, x[:-3], out_path + sample_name + '.INV.' + key_event.replace(':', '__') + '.png')
        elif x[-1] in ['INS']:
            key_event = ':'.join([str(i) for i in x[:-3] + ['INS']])
            ins_pos = '_'.join([str(i) for i in x[:2]])
            POLARITY = '+'
            if type(x[4]) == type(4):
                ins_seq = ''.join(['X' for i in range(x[4])])
            else:
                ins_seq = x[4]
            vapor_score_event = sf.vapor_simple_ins_Vapor(num_reads_cff, plt_li, bam_in, ref, ins_pos, ins_seq, out_path + sample_name + '.INS.' + key_event.replace(':', '__') + '.png', POLARITY)
        elif x[-1] in ['a/aa', 'aa/a', 'aa/aa', 'DUP', 'TANDUP']:
            key_event = ':'.join([str(i) for i in x[:-3]] + ['TANDUP'])
            vapor_score_event = sf.vapor_simple_tandup_Vapor(num_reads_cff, plt_li, bam_in, ref, x[:-3], out_path + sample_name + '.TANDUP.' + key_event.replace(':', '__') + '.png')
        else:
            print(x)
            return None, buf.getvalue()
        vapor_result = sf.result_organize_ins([key_event, vapor_score_event])
        row = vapor_result[0].split(':') + [x[3]] + vapor_result[1:]
        # the original printed after writing the row; stdout and the output file are separate
        # streams, so replaying the text after the row is written reproduces both exactly
        print(vapor_result)
    return row, buf.getvalue()


def _bed_tasks(bed_info, num_reads_cff, bam_in, ref, out_path, sample_name):
    plt_li = 0
    for x in bed_info:
        known = x[-1] in ['a/', '/a', '/', 'DEL', 'a/a^', 'a^/a', 'a^/a^', 'INV', 'INS',
                          'a/aa', 'aa/a', 'aa/aa', 'DUP', 'TANDUP']
        if known:
            plt_li += 1   # figure number; the original only counted recognised records
        yield (plt_li, x, num_reads_cff, bam_in, ref, out_path, sample_name)


def run_bed(bed_info, out_name, num_reads_cff, bam_in, ref, out_path, sample_name, threads=1):
    from vapor_vali import Simple_function as sf
    tasks = _bed_tasks(bed_info, num_reads_cff, bam_in, ref, out_path, sample_name)
    with open(out_name, 'a') as fo:
        def emit(result):
            row, text = result
            if row is not None:
                sf.write_output_row(fo, row)
                fo.flush()
            sys.stdout.write(text)

        if threads <= 1:
            for t in tasks:
                emit(_score_bed_record(t))
            return
        ctx = multiprocessing.get_context(os.environ.get('VAPOR_MP_START') or None)
        with ctx.Pool(threads) as pool:
            # imap keeps input order; chunksize 1 balances SVs of very different cost
            for result in pool.imap(_score_bed_record, tasks, chunksize=1):
                emit(result)
