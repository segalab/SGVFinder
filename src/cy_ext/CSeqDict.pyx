from libc.stdio cimport *
from libcpp.string cimport string
from libcpp.unordered_map cimport unordered_map
from libcpp.pair cimport pair
from libc.stdlib cimport malloc, free
cimport cython
from cpython.mem cimport PyMem_Malloc, PyMem_Free
import ujson
from os.path import splitext

cdef char *basemap = <char *>PyMem_Malloc(255 * cython.sizeof(char))
cdef int c_i = 0
while c_i < 255:
    basemap[c_i] = 'N'
    c_i += 1
basemap[65] = 'T'
basemap[67] = 'G'
basemap[71] = 'C'
basemap[84] = 'A'
basemap[97] = 'T'
basemap[99] = 'G'
basemap[103] = 'C'
basemap[116] = 'A'
basemap[78] = 'N'

PLUS = '+'

cdef SeqIDLL * newSeqID(string seqid):
    cdef SeqIDLL *new_ = new SeqIDLL();
    new_.name = seqid
    new_.newer = EMPTY
    new_.older = EMPTY
    return new_
    
cdef struct TestS:
    string name
 
cpdef int test(string n) except -1:
    cdef TestS cur
    cur.name = n
    return 0
    

cdef char * revcomp(unsigned char *seq, long seqlen):
    cdef char *ret = <char *> PyMem_Malloc(seqlen * cython.sizeof(char) + 1) 
    if ret == NULL:
        print 'mem in revcomp'
    cdef long i = 0
    while i < seqlen:
        ret[i] = basemap[seq[seqlen-i-1]]
        i+=1
    ret[seqlen] = '\0'
    return ret

cdef SeqIDLL *EMPTY

cdef class CSeqDict:
    
    def __init__(self, pth, lim = None):
        with open(splitext(pth)[0] + '.lens') as ol:
            for k, v in ujson.load(ol).iteritems():
                self._lendct[k] = v 
        with open(splitext(pth)[0] + '.idx') as o:
            for k, v in ujson.load(o).iteritems():
                self._idx[k] = v
        self._cursum = 0
        self._lim = lim
        self._file = fopen(pth, "rb")
        if self._file == NULL:
            raise Exception(2, "No such file or directory: '%s'" % pth)
    
    def get_lendct(self):
        return self._lendct
    
    def __dealloc__(self):
        fclose(self._file)
    
    def _get_rev(self, seqid):
        return self.get_rev(seqid)
    
    def _get_seq(self, seqid):
        return self.get_seq(seqid, 1)
    
    def get_sequence(self, gid, strand, pos, lngt):
        ln = self.get_len(gid)
        if strand == PLUS:
            sq, st, en, strnd =  self.get_seq(gid, 1), pos-1, pos+lngt-1, 1 
        else:
            sq = self.get_rev(gid)
            st, en, strnd =  ln-pos-lngt+1, ln-pos+1, -1
        return sq, st%ln, en%ln, strnd
    
    cdef char * get_seq(self, string seqid, int callcache) except NULL:
        cdef char *sq
        if callcache:
            self._add_to_cacheq(seqid)
        if self._dct.count(seqid) > 0:
            return self._dct[seqid]
        else:
            sq = self._load_seq(seqid)
            self._dct[seqid] = sq
            self._cursum += self._lendct[seqid]
            self._clean_mem()
            return self._dct[seqid]
    
    cdef char * get_rev(self, string seqid) except NULL:
        cdef long slen
        cdef char *sq
        self._add_to_cacheq(seqid)
        if self._revdct.count(seqid) > 0:
            return self._revdct[seqid]
        else:
            slen = self._lendct[seqid]
            sq = revcomp(<unsigned char *>self.get_seq(seqid, 0), slen)
            self._revdct[seqid] = sq
            self._cursum += slen
            self._clean_mem()
            return self._revdct[seqid]
        
    cpdef long get_len(self, string seqid):
        return self._lendct[seqid]
    
    cpdef long get_cursum(self):
        return self._cursum

    cdef int _add_to_cacheq(self, string seqid) except -1:
        cdef SeqIDLL *cur
        if self._cppcq.count(seqid) > 0:
            cur = self._cppcq[seqid]
            if cur.newer == EMPTY:
                return 0
            if cur.older == EMPTY:
                cur.newer.older = EMPTY
                self.oldest = cur.newer
            else:
                cur.older.newer = cur.newer
                cur.newer.older = cur.older
            self.newest.newer = cur
            cur.newer = EMPTY
            cur.older = self.newest
            self.newest = cur
        else:
            cur = newSeqID(seqid)
            cur.newer = EMPTY
            cur.older = EMPTY
            if self.newest == NULL:
                self.newest = cur
                self.oldest = cur
            else:
                self.newest.newer = cur
                cur.newer = EMPTY
                cur.older = self.newest
                self.newest = cur
            self._cppcq[seqid] = cur
        return 0
    
    cdef char * _load_seq(self, string seqid):
        cdef long pos = self._idx[seqid]
        cdef char * line = NULL
        cdef size_t l = 0
        fseek(self._file, pos, SEEK_SET)
        getline(&line, &l, self._file)
         
        cdef long lseq = self._lendct[seqid]
        cdef long ind = lseq - 2
        while line[ind] != 10 and line[ind] != 0:
            ind += 1
        line[ind] = '\0' 
        return line 
    
    cdef int _clean_mem(self) except -1:
        cdef string seqtodl
        cdef SeqIDLL *new_oldest = NULL
        while self._cursum > self._lim:
            seqtodl = self.oldest.name
            self._cppcq.erase(seqtodl)
            PyMem_Free(self._dct[seqtodl])
            self._dct.erase(seqtodl)
            self._cursum -= self._lendct[seqtodl]
            if self._revdct.count(seqtodl) > 0:
                PyMem_Free(self._revdct[seqtodl])
                self._revdct.erase(seqtodl)
                self._cursum -= self._lendct[seqtodl]
            if self.oldest.newer != NULL:
                new_oldest = self.oldest.newer
                new_oldest.older = EMPTY
            PyMem_Free(self.oldest)
            self.oldest = new_oldest
            if new_oldest == NULL:
                self.newest = NULL
        return 0
