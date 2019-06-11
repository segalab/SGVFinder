from glob import glob
from os.path import join, split, basename, realpath, exists, splitext
import logging
import numpy as np
from datetime import datetime
from collections import defaultdict
from ICRAUtils import timeit
from pandas import to_pickle, read_pickle
import cy_ext.gem_qa as gqa
from GEMMapping import gemmap_single_file, _load_from_file
from cy_ext import CSeqDict
import ujson
from ReadContainer import ReadContainer
log_ = logging.getLogger('ICRA')

GENOMES = 'genomes'
CAMI = 'cami'

def single_file(fq1, fq2, outfol, max_mismatch, consider_lengths, epsilon, \
                max_iterations, min_bins = 4, max_bins = 100, \
                min_reads = 5, dense_region_coverage = 60, length_minimum = 300, \
                length_maximum = 1e6, usage = 'cami', use_theta = False):
    assert usage in {CAMI, GENOMES}
    log_.debug('Loading database...')
    seqdict = C_CamidDBDict() if usage == CAMI else C_EmblDB() 
    indexl = cami_index if usage == CAMI else genomes_index 
    dest_dictf =  cami_destdictf if usage == CAMI else genomes_destdictf
    pmpf = join(outfol, basename(fq1).replace('_1.fastq', '.pmp').replace('.fastq', '.pmp'))
    log_.debug('Loaded. Mapping file...')
    gemmap_single_file(fq1, fq2, pmpf, indexl, seqdict, 0.05, 0.05, 2)
    log_.debug('Mapped. Running ICRA...') 
    read_container, pi, theta1, average_read_length, lengthdb = \
        _initialize(fq1, fq2, pmpf, seqdict, max_mismatch, consider_lengths, \
                    length_minimum, length_maximum, min_reads, min_bins, max_bins, \
                    dense_region_coverage, dest_dictf, use_theta)
    del seqdict
    if len(pi) == 0:
        delta = {}
    else:
        delta, pi = _runIterative(pi, theta1, read_container, min_bins, max_bins, \
                              min_reads, average_read_length, max_mismatch, \
                              lengthdb, dense_region_coverage, consider_lengths, \
                              length_minimum, length_maximum, epsilon, max_iterations, use_theta)
    with open(join(outfol, basename(pmpf).replace('.pmp', '.jsdel')), 'wb') as of:
        ujson.dump(delta, of, double_precision=100)
    return pmpf.replace('.pmp', '.jsdel')

@timeit(log_, "%s: Parameter init complete. Time: %s", logging.INFO)    
def _initialize(fq1, fq2, pmpf, seqdict, max_mismatch, consider_lengths, \
               length_minimum, length_maximum, min_reads, min_bins,  
               max_bins, dense_region_coverage, dest_dictf, use_theta): 
    dest_dict = read_pickle(dest_dictf) if dest_dictf is not None else None
    average_read_length = (_getaveragereadlength(fq1)*2+500) if fq2 is not None else _getaveragereadlength(fq1)
    lengthdb = _get_len_dct(seqdict, dest_dict)
    read_container = ReadContainer()
    pe = fq2 is not None 
    allowed_dests = {k for k,v in lengthdb.iteritems() if length_minimum <= v <= length_maximum}
    for sr in _load_from_file(pmpf): 
        read_container.add_pmp_sr(sr, pe, allowed_dests, max_mismatch, None, seqdict, dest_dict)
    del dest_dict, allowed_dests, seqdict
    pi = read_container.get_pi()
    
    read_container.remove_dests([k for k,val in pi.iteritems() if val<=min_reads])
    pi = {k:val for k, val in pi.iteritems() if  val > min_reads}
    
    read_container.init_cov_dict(pi, min_bins, max_bins, min_reads, average_read_length, lengthdb)
    #note: this will create a fully updated coverage dict in the process
    #(in the original implementation this was explicit and not implicit)
    pi = read_container.get_dense_region_coverage_pi(dense_region_coverage) 
    
    read_container.remove_dests([k for k,val in pi.iteritems() if val<=min_reads])
    pi =  {k:v for k,v in pi.iteritems() if v>min_reads}
    
    if consider_lengths:
        pi = {k:(v/(lengthdb[k] - (2 * average_read_length))) for k, v in pi.iteritems()}
    pisum = sum(pi.itervalues())
    pi = {k:v/float(pisum) for k,v in pi.iteritems()}
    
    theta1 = read_container.get_theta() if use_theta else None
    
    return read_container, pi, theta1, average_read_length, lengthdb
    
def _runIterative(pi, theta1, read_container, min_bins, max_bins, min_reads, average_read_length, \
                  max_mismatch, lengthdb, dense_region_coverage, consider_lengths, length_minimum, 
                  length_maximum, epsilon, max_iterations, use_theta):
    starttime = datetime.now()
    pi_dict = {}
    i = 1
    
    while True:
        pi_dict[i] = pi
        delta = read_container.do_delta(pi, theta1, 1. / ((len(pi)**2) * (10**max_mismatch)), use_theta)
        read_container.recalc_cov_bins()
        prevPi = pi
        #reminder: implicitly creates an updated covdic
        pi = read_container.get_dense_region_coverage_pi(dense_region_coverage, delta)
        read_container.remove_dests([k for k,val in pi.iteritems() if val<min_reads])
        pi =  {k:v for k,v in pi.iteritems() if v>=min_reads}
        
        if len(pi) == 0:
            log_.info("No adequately covered strains found")
            return {}, {}
        
        theta1 = read_container.get_theta() if use_theta else None
        
        if consider_lengths:
            pi = {k:(v/(lengthdb[k] - (2 * average_read_length))) for k, v in pi.iteritems()}
        pisum = sum(pi.itervalues())
        pi = {k:v/float(pisum) for k,v in pi.iteritems()}
        
        prevPi = {k:v for k, v in prevPi.iteritems() if k in pi}
        dPi = _LogDistDict(pi, prevPi)
         
        i += 1
        log_.info("Iteration {} - Time: {}, dPi = {:.2e}, nPi = {}".format(i, datetime.now() - starttime, dPi, 
                                                                            len(pi)))
         
        if  dPi < epsilon or i > max_iterations:
            break
    log_.info("Final result - Time: {}".format(datetime.now() - starttime))
    return read_container.get_full_delta(pi, theta1, 1. / ((len(pi)**2) * (10**max_mismatch)), use_theta), pi

def _write_if_vb(vb, pth, obj):
    if vb:
        to_pickle(obj, pth)

def _calc_qual(sr, pmap_map, seqdict):
        if hasattr(pmap_map, 'mismatches'):
            ga1, off1 = _get_opt_gigar_array(seqdict, pmap_map.dest_id, pmap_map.strand, pmap_map.pos, pmap_map.mismatches, sr.seq[pmap_map.end], pmap_map.gigar)
            return 10**gqa.calc_quality(sr.quals[pmap_map.end], ga1), off1, None
        else:
            ga1, off1 = _get_opt_gigar_array(seqdict, pmap_map.dest_id, pmap_map.strand1, pmap_map.pos1, pmap_map.mismatches1, sr.seq[0], pmap_map.gigar1)
            ga2, off2 = _get_opt_gigar_array(seqdict, pmap_map.dest_id, pmap_map.strand2, pmap_map.pos2, pmap_map.mismatches2, sr.seq[1], pmap_map.gigar2)
            return 10**(gqa.calc_quality(sr.quals[0], ga1) + gqa.calc_quality(sr.quals[1], ga2)), off1, off2

def _logdist(u, v):
    x = np.abs(np.log10(np.asarray(u)) - np.log10(np.asarray(v)))
    return np.percentile(x, 90)

def _LogDistDict(u, v):
    x = []
    y = []
    if not (set(u.keys()).intersection(set(v.keys())) == set(u.keys()) == set(v.keys())):
        raise KeyError("Keys not overlapping between the two dictionaries")
    for key in u:
        if u[key]>1e-5 or v[key]>1e-5:
            x.append(u[key])
            y.append(v[key])
    return _logdist(x,y)

def _getaveragereadlength(fq):
    with open(fq) as inf:
        sum_bps = 0
        for i, read in enumerate(inf):
            if i % 4 == 1:
                sum_bps += len(read) - 1 #has a /n at the end
            if i >= 4000: break
    arl = int(np.round(sum_bps/1000))
    log_.info("Average read length for {} set to {}".format(basename(fq), arl))
    return arl

def _get_len_dct(sequence_dict = None, read_dest_dict = None):
    if hasattr(sequence_dict, 'get_lendct'):
        lengthdb = sequence_dict.get_lendct()
    else:
        lengthdb = {k:sequence_dict.get_len(k) for k in sequence_dict._dct.keys()}
    if read_dest_dict is not None:
        newlengthdb = defaultdict(int)
        for l in lengthdb:
            newlengthdb[read_dest_dict[l][0]] += lengthdb[l]
        return dict(newlengthdb)
    else:
        return lengthdb

def _get_opt_gigar_array(seqdict, dest_id, strand, pos, mismatches, seq, gigar):
    if '>' in gigar:
        lseq = len(seq)
        genesq, start, end, strand = seqdict.get_sequence(dest_id, strand, pos, lseq)
        return gqa.get_opt_gigar_array(seqdict.get_len(dest_id),genesq, start, end,
                                                                 strand, gigar, lseq, seq, mismatches)
    else:
        return gqa.create_gigar_array(gigar), 0
  
def C_CamidDBDict():
    pth = join(split(realpath(__file__))[0], '../DataFiles/2017-10-CAMI_genomes.filtered.if_ujsn')
    assert exists(pth), 'Could not find database file at ' + pth
    return CSeqDict.CSeqDict(pth, 9e9)

def C_EmblDB():
    pth = join(split(realpath(__file__))[0], '../DataFiles/representatives.contigs.drepped.if_ujsn')
    assert exists(pth), 'Could not find database file at ' + pth
    assert exists(splitext(pth)[0] + '.lens'), 'Could not find database file at ' + splitext(pth)[0] + '.lens'
    assert exists(splitext(pth)[0] + '.idx'), 'Could not find database file at ' + splitext(pth)[0] + '.idx'
    return CSeqDict.CSeqDict(pth, 18e9)

cami_index = [join(split(realpath(__file__))[0], '../DataFiles/2017-10-CAMI_genomes.filtered.gem')]

cami_destdictf = join(split(realpath(__file__))[0], '../DataFiles/2017-10-CAMI_genomes.filtered.dests')

genomes_index = sorted(glob(join(split(realpath(__file__))[0], '../DataFiles/GEMSplit/*.gem')))

genomes_destdictf = join(split(realpath(__file__))[0], '../DataFiles/representatives.contigs.drepped.dests')

    