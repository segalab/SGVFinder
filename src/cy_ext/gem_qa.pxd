cimport cython
from libc.math cimport pow, log10
from cpython.mem cimport PyMem_Malloc, PyMem_Free
import numpy as np
cimport numpy as np
cdef extern from "stdlib.h":
    int isdigit (int c)
from libc.string cimport memset
cdef int * create_gigar_array_ptr(char *gigar, int *size)
cdef int * get_opt_gigar_array_ptr(int lg, char *genesq, int start, int end, 
                        int strand, char *gigar, int lseq, 
                         char *seq, int mismatches, int * size, int * offst)
