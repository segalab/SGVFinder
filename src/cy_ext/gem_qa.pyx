cimport cython
from libc.math cimport pow, log10
from cpython.mem cimport PyMem_Malloc, PyMem_Free
import numpy as np
cimport numpy as np
cdef extern from "stdlib.h":
    int isdigit (int c)
from libc.string cimport memset
__version__ = '0.0.1'
cdef int TEN = 10
cdef int M1 = -1
cdef int P1 = 1
cdef int ZERO = 0
cdef int BASE_10 = 10
cdef int ASCII_0 = 48
cdef char N = 'N'
cdef char Y = 'Y'
cdef char R = 'R'
cdef double * ps = < double *> PyMem_Malloc(42 * cython.sizeof(double))
cdef double * qs = < double *> PyMem_Malloc(42 * cython.sizeof(double))
if ps is NULL or qs is NULL:
    raise MemoryError()
cdef int x = 0
while x < 42:
    ps[x] = log10(1 - pow(10, -0.1 * x))
    qs[x] = -0.1 * x
    x += 1

np.import_array()

ctypedef np.int32_t DTYPE_t

cdef extern from "numpy/arrayobject.h":
    void PyArray_ENABLEFLAGS(np.ndarray arr, int flags)

cdef data_to_numpy_array_with_spec(void * ptr, np.npy_intp N, int t):
    cdef np.ndarray[DTYPE_t, ndim = 1] arr = np.PyArray_SimpleNewFromData(1, & N, t, ptr)
    PyArray_ENABLEFLAGS(arr, np.NPY_OWNDATA)
    return arr

cdef int * _help = < int *> PyMem_Malloc(1000 * cython.sizeof(int))

cdef int * create_gigar_array_ptr(char * gigar, int * size):
    cdef int h_i = 0
    cdef int gigar_len = len(gigar)
    cdef int gigar_i = 0
    cdef int curchr, value, _v
    while gigar_i < gigar_len:
        curchr = gigar[gigar_i]
        if not isdigit(curchr):
            if curchr != '>':
                _help[h_i] = TEN
                h_i += 1
            gigar_i += 1
        else:
            value = 0
            while isdigit(curchr):
                value = value * BASE_10 + (curchr - ASCII_0)
                gigar_i += 1
                if (gigar_i == gigar_len):
                    break
                curchr = gigar[gigar_i]
            if curchr == '-':
                _v = M1
                gigar_i += 1
            elif curchr == '+':
                _v = P1
                gigar_i += 1
            else:
                _v = 0
            h_i += value
            while value > 0:
                _help[h_i - value] = _v
                value -= 1
    cdef int * res = < int *> PyMem_Malloc(h_i * cython.sizeof(int))
    if res == NULL:
        print 'memory in create_gi_array'
    _v = 0
    while _v < h_i:
        res[_v] = _help[_v]
        _v += 1
    size[0] = h_i
    return res

cdef int opt_hamming(char * seq, int lseq, char * genesq, int lgseq, int rstart, int rend) except -100:
    cdef int h = 0
    cdef int o
    if rstart < rend:
        o = rstart
        while rstart < rend:
            if rstart - o < 0 or rstart - o > lseq or rstart < 0 or \
                    rstart > lgseq or seq[rstart - o] != genesq[rstart]:
                h += 1
            rstart += 1
        return h    
    else:
        if lgseq-rstart > rend:
            return opt_hamming(seq, lseq, genesq, lgseq, rstart, lgseq + rend)
        else:
            return opt_hamming(seq, lseq, genesq, lgseq, 0 - (lgseq - rstart), rend)

cdef int * opt_hamming_gigar_arr(char * seq, char * genesq, int rstart, int rend, int * size, int lgsq):
    if rstart > rend:
        if lgsq-rstart > rend:
            rend = lgsq + rend
        else:
            rstart = 0 - (lgsq - rstart)
    cdef int i = 0
    cdef int * res = < int *> PyMem_Malloc((rend - rstart) * cython.sizeof(int))
    if res == NULL:
        print 'memory in opt hamming'
    while i + rstart < rend:
        if i + rstart < 0 or i + rstart > lgsq or seq[i] != genesq[i + rstart]:
            res[i] = TEN
        else:
            res[i] = ZERO
        i += 1
    size[0] = rend - rstart
    return res

cdef int * get_opt_gigar_array_ptr(int lg, char * genesq, int start, int end,
                        int strand, char * gigar, int lseq,
                         char * seq, int mismatches, int * size, int * offst):
    cdef int len_ga;
    cdef int * ga = create_gigar_array_ptr(gigar, & len_ga)
    cdef int curroffset = 0
    cdef int prev = 0
    cdef int g
    cdef int i = 0
    cdef int best_of = 0
    cdef int best_of_dst = 1000
    if (end - start) == lseq or end + lg - start == lseq:
        best_of_dst = opt_hamming(seq, lseq, genesq, lg, start, end)
    cdef int cur_dst
    
    while i < len_ga:
        g = ga[i]
        if g != prev:
            if (prev == 1 or prev == -1) and \
               (end + strand * curroffset) % lg - (start + strand * curroffset) == lseq:
                cur_dst = opt_hamming(seq, lseq, genesq, lg, start + strand * curroffset, end + strand * curroffset)
                if cur_dst < best_of_dst:
                    best_of_dst = cur_dst
                    best_of = curroffset
            prev = g
        if g == 1 or g == -1:
            curroffset += g
        i += 1
    if (prev == 1 or prev == -1) and \
       (end + strand * curroffset) % lg - (start + strand * curroffset) == lseq:
        cur_dst = opt_hamming(seq, lseq, genesq, lg, start + strand * curroffset, end + strand * curroffset)
        if cur_dst < best_of_dst:
            best_of_dst = cur_dst
            best_of = curroffset
    if best_of_dst <= mismatches:
        PyMem_Free(ga)
        offst[0] = best_of
        return opt_hamming_gigar_arr(seq, genesq, start + strand * best_of, (end + strand * best_of) % lg, size, lg)
    else:
        size[0] = len_ga
        offst[0] = 0
        return ga

cpdef get_opt_gigar_array(int lg, char * genesq, int start, int end,
                        int strand, char * gigar, int lseq,
                         char * seq, int mismatches):
    cdef int size, ofst
    cdef int * ga = get_opt_gigar_array_ptr(lg, genesq, start, end,
                                           strand, gigar, lseq, seq,
                                           mismatches, & size, & ofst)
    return data_to_numpy_array_with_spec(ga, size, np.NPY_INT32), ofst

cpdef create_gigar_array(char * gigar):
    cdef int h_i = 0
    cdef int gigar_len = len(gigar)
    cdef int gigar_i = 0
    cdef int curchr, value, _v
    while gigar_i < gigar_len:
        curchr = gigar[gigar_i]
        if not isdigit(curchr):
            if curchr != '>':
                _help[h_i] = TEN
                h_i += 1
            gigar_i += 1
        else:
            value = 0
            while isdigit(curchr):
                value = value * BASE_10 + (curchr - ASCII_0)
                gigar_i += 1
                if (gigar_i == gigar_len):
                    break
                curchr = gigar[gigar_i]
            if curchr == '-':
                _v = M1
                gigar_i += 1
            elif curchr == '+':
                _v = P1
                gigar_i += 1
            else:
                _v = 0
            h_i += value
            while value > 0:
                _help[h_i - value] = _v
                value -= 1
    cdef int * res = < int *> PyMem_Malloc(h_i * cython.sizeof(int))
    _v = 0
    while _v < h_i:
        res[_v] = _help[_v]
        _v += 1
    return data_to_numpy_array_with_spec(res, h_i, np.NPY_INT32)

@cython.boundscheck(False)
cpdef double calc_quality(quals, gigar):
    cdef int ql = len(quals)
    cdef int * qa = < int *> PyMem_Malloc(ql * cython.sizeof(int))
    cdef int * ga = < int *> PyMem_Malloc(ql * cython.sizeof(int))
    if qa is NULL or ga is NULL:
        raise MemoryError()
    cdef int i = 0
    while i < ql:
        qa[i] = quals[i]
        ga[i] = gigar[i]
        i += 1
    return _calc_quality(qa, ga, ql)

@cython.boundscheck(False)
cdef double _calc_quality(int * quals, int * gigar, int ql):
    cdef int i = 0
    cdef double res = 0
    while i < ql:
        if gigar[i] == 0:
            res += ps[quals[i]]
        elif gigar[i] == -1 or gigar[i] == 10:
            res += qs[quals[i]]
        i += 1
    PyMem_Free(quals)
    PyMem_Free(gigar)
    return res

@cython.boundscheck(False)
cpdef int sum_skips(int[:] gigar):
    cdef int sum = 0
    cdef int i = 0
    cdef int gigar_len = len(gigar)
    cdef int g
    while i < gigar_len:
        g = gigar[i]
        if g != TEN:
            sum += g
        i += 1
    return sum

cpdef int validate_match3(char * seq, char * destseq, int[:] gigar):
    cdef int seq_len = len(seq)
    cdef int destseq_len = len(destseq)
    cdef int seq_i = 0
    cdef int destseq_i = 0
    cdef int gigar_i = 0
    cdef int g
    
    while gigar_i < seq_len:
        g = gigar[gigar_i] 
        if g == ZERO and seq[seq_i] != destseq[destseq_i] and not (seq[seq_i] == N and (destseq[destseq_i] == R or destseq[destseq_i] == Y)):
            return 0
        elif g == M1:
            destseq_i -= 1
        elif g == P1:
            seq_i -= 1
        seq_i += 1
        destseq_i += 1
        gigar_i += 1
    return 1
