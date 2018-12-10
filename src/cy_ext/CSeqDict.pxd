from libc.stdio cimport *
from libcpp.string cimport string
from libcpp.unordered_map cimport unordered_map
from libcpp.pair cimport pair
from libc.stdlib cimport malloc, free
cimport cython
from cpython.mem cimport PyMem_Malloc, PyMem_Free
import ujson
from os.path import splitext

cpdef cppclass SeqIDLL:
    string name
    SeqIDLL *newer
    SeqIDLL *older
    
cdef class CSeqDict:
    cdef unordered_map[string, char *] _dct
    cdef unordered_map[string, char *] _revdct
    cdef unordered_map[string, long] _lendct
    cdef unordered_map[string, long] _idx 
    cdef unordered_map[string, SeqIDLL *] _cppcq
    cdef long _cursum, _lim
    cdef FILE* _file
    cdef SeqIDLL *newest
    cdef SeqIDLL *oldest
    cdef int _add_to_cacheq(self, string seqid) except -1
    cdef char * _load_seq(self, string seqid)
    cdef int _clean_mem(self) except -1
    cdef char * get_seq(self, string seqid, int callcache) except NULL
    cdef char * get_rev(self, string seqid) except NULL
    cpdef long get_len(self, string seqid)
    cpdef long get_cursum(self)
    