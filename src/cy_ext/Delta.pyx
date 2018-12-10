from __future__ import division
cimport cython
from itertools import groupby
from operator import itemgetter
from collections import defaultdict

cdef: 
    int r_DEST      = 0
    int r_QUALITY   = 1
    int r_POS1      = 2
    int r_POS2      = 3

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.initializedcheck(False)
@cython.cdivision(True)
def do_delta(pi, theta1, read_container, includedbacs, epsilon_method, max_mismatch, b_thresh):
    cdef:
        char* read_id
        double thetamul, readprod, denom_delta_for_read, pidest, qual, delrid
        double thetamin = theta1['_Minimum']
        double epsilon = 1e-18 
        double ep_factor = 1
        int i_, pos1, pos2, lreads
        set includedbacs_ = includedbacs
        set thetakeys = set(theta1.keys())
    if (epsilon_method == "genomes"):
        epsilon = 1e-13
    elif (epsilon_method == "genes"):
        epsilon = 1e-16 
    elif (epsilon_method == "dynamic"):
        ep_factor = (len(pi)**2) * (10**max_mismatch)  
        epsilon = 1./ep_factor
    delta = defaultdict(lambda: defaultdict(float))
    winning_hard = defaultdict(int)
    for read_id, reads in read_container:
        lreads = len(reads)
        if lreads == 0:
            continue
        if lreads == 1 or reads[0][r_QUALITY] > reads[1][r_QUALITY] + b_thresh:
            winning_hard[reads[0][r_DEST]] += 1.
        denom_delta_for_read = 0.
        i_ = 0
        while i_ < lreads:
            read = reads[i_]
            dest = read[r_DEST]
            if dest in includedbacs_:
                qual = read[r_QUALITY]
                pidest = pi[dest]
                if dest in thetakeys:
                    thetamul = theta1[dest]
                else:
                    thetamul = thetamin
                readprod = pidest * thetamul * qual
                delta[read_id][(dest,read[r_POS1],read[r_POS2])] += readprod
                denom_delta_for_read += readprod
            i_ += 1
        for germ in delta[read_id].iterkeys():
            delrid = delta[read_id][germ] 
            if delrid < epsilon:
                denom_delta_for_read -= delrid
        for germ in delta[read_id].iterkeys():
            delrid = delta[read_id][germ] 
            if delrid < epsilon: 
                delta[read_id][germ] = 0
            else:
                delta[read_id][germ] = delrid / denom_delta_for_read
    return delta, winning_hard
