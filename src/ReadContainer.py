import cy_ext.gem_qa as gqa
import numpy as np
from ICRAUtils import timeit
import logging

log_ = logging.getLogger('ICRA')

r_QUALITY   = 0
r_POS1      = 1
r_POS2      = 2
_r_DEST     = 3 

class ReadContainer(object):
    def __init__(self):
        #simple conversion between a running number and the string
        self._did_to_dest = dict()
        self._dest_to_did = dict()
        self._next_did = 0 #did => dest id
        
        #holds all the ambiguous reads, accessed by their rrid
        self._ambig_rds = dict()
        #holds the ambiguous rrids belonging to each did. NOTE: This is not guaranteed to be updated!
        self._did_to_ambig_rrids = dict()
        #holding the unique reads (not rids) belonging to each did
        self._did_to_uniq_rds = dict()
        #hold the original rid of each rrid
        self._rrid_to_rid = dict()
        self._next_rrid = 0 #rrid => running read id
        
        #the pi dict composed of the unique reads - only an entire did can be removed from it
        self._uniq_pi = dict()
        
        self._unique_cov_dic = None
    
    def _id_from_dest(self, dest):
        try:
            return self._dest_to_did[dest]
        except KeyError:
            self._did_to_dest[self._next_did] = dest
            self._dest_to_did[dest] = self._next_did
            self._next_did += 1
            return self._next_did - 1
    
    def _get_rrid(self, rid): #make sure to call only once per rid, doesn't retrieve existing ones
        self._rrid_to_rid[self._next_rrid] = rid
        self._next_rrid += 1
        return self._next_rrid - 1
    
    def _add_unique(self, pos1, pos2, did, rid):
        rrid = self._get_rrid(rid)
        self._add_unique_rrid(pos1, pos2, did, rrid)
        
    def _add_unique_rrid(self, pos1, pos2, did, rrid):
        try:
            self._did_to_uniq_rds[did].append((rrid, pos1, pos2))
        except KeyError:
            self._did_to_uniq_rds[did] = [(rrid, pos1, pos2)]
        try:
            self._uniq_pi[did] += 1
        except KeyError:
            self._uniq_pi[did] = 1
        if self._unique_cov_dic is not None:
            _add_read(self._unique_cov_dic[did], pos1, pos2, self.avg_read_length, 1, self._reduced_len_dct[did])
            
    def _add_ambig(self, rid, mpings):
        rrid = self._get_rrid(rid)
        self._ambig_rds[rrid] = mpings
        for m in mpings:
            try:
                self._did_to_ambig_rrids[m[_r_DEST]].append(rrid)
            except KeyError:
                self._did_to_ambig_rrids[m[_r_DEST]] = [rrid]
    
    def add_pmp_sr(self, sr, pe, allowed_dests, max_mismatch, filtermaps, seqdict, dest_dict):
        cur = []
        if pe is None or pe: 
            for pemap in sr.pe_maps:
                if pemap.mismatches1 + pemap.mismatches2 <= max_mismatch:
                    newdest, addpos = dest_dict[pemap.dest_id] if dest_dict is not None else (pemap.dest_id, 0)
                    if newdest in allowed_dests:
                        cur.append((pemap, newdest, addpos))
        if pe is None or not pe:
            for semap in sr.se_maps:
                if semap.mismatches <= max_mismatch:
                    newdest, addpos = dest_dict[semap.dest_id] if dest_dict is not None else (semap.dest_id, 0)
                    if newdest in allowed_dests:
                        cur.append((semap, newdest, addpos))
                    
        if len(cur) == 1:
            mp, newdest, addpos = cur[0]
            self._add_unique(addpos + (mp.pos if hasattr(mp, 'pos') else mp.pos1), \
                             -1 if hasattr(mp, 'pos') else (addpos + mp.pos2), self._id_from_dest(newdest), sr.rid)
        elif len(cur) > 1:
            cur2 = []
            for mp, newdest, addpos in cur:
                qual, off1, off2 = _calc_qual(sr, mp, seqdict)
                cur2.append((qual, addpos + off1 + (mp.pos if hasattr(mp, 'pos') else mp.pos1), \
                             -1 if hasattr(mp, 'pos') else (addpos + off2 + mp.pos2), self._id_from_dest(newdest)))
            cur2.sort(key = lambda m:m[0], reverse = True)
            if filtermaps:
                cur2 = [m for m in cur2 if m[0]*filtermaps > cur2[0][0]]
            if len(cur2) == 1: #filtermaps has to run after quals are calculated, and it could result in a unique rd 
                self._add_unique(cur2[0][1], cur2[0][2], cur2[0][3], sr.rid)
            else:
                self._add_ambig(sr.rid, cur2)
    
    def get_pi(self):
        pi = self._uniq_pi.copy()
        for temp_reads in self._ambig_rds.itervalues():
            qual_sum = sum(r[r_QUALITY] for r in temp_reads)
            for r in temp_reads:
                try:
                    pi[r[_r_DEST]] += r[r_QUALITY] / qual_sum
                except KeyError:
                    pi[r[_r_DEST]] = r[r_QUALITY] / qual_sum
        return {self._did_to_dest[k]:v for k,v in pi.iteritems()}
    
    def get_theta(self):
        win_hard = self._uniq_pi.copy() #win hard differs from pi only in its handling of ambig reads
        winning_hard = {k:v for k,v in win_hard.iteritems()}
        winning_hard['_Minimum'] = 1
        theta1sum = sum(winning_hard.values())
        return {k:float(v)/theta1sum for k,v in winning_hard.iteritems()}
    
    def remove_dests(self, dests):
        #some of the dels here are not required, but make the data structure more tidy and slightly smaller (memory-wise)
        dids = [self._dest_to_did[dest] for dest in dests]
        for did in dids:
            if did not in self._did_to_ambig_rrids: continue
            for rrid in self._did_to_ambig_rrids[did]:
                if rrid not in self._ambig_rds:
                    continue
                new_mping = [m for m in self._ambig_rds[rrid] if m[_r_DEST] != did]
                #note that if that single mapping is to a dest that should also be removed - it will 
                #be when we handle unique mappings - and this is also why we loop through the dests 
                #twice, first going through all ambiguous and then through all unique
                if len(new_mping) == 1:
                    self._add_unique_rrid(new_mping[0][r_POS1], new_mping[0][r_POS2], new_mping[0][_r_DEST], rrid)
                    del self._ambig_rds[rrid]
                elif len(new_mping) > 1:
                    self._ambig_rds[rrid] = new_mping
                else:
                    del self._ambig_rds[rrid]
            del self._did_to_ambig_rrids[did]
        for did in dids:
            if did in self._uniq_pi:
                del self._uniq_pi[did]
                for rd in self._did_to_uniq_rds[did]:
                    del self._rrid_to_rid[rd[0]]
                del self._did_to_uniq_rds[did] 
            del self._did_to_dest[did]
            if self._unique_cov_dic is not None:
                del self._unique_cov_dic[did]
        for dest in dests:
            del self._dest_to_did[dest]
    
    def init_cov_dict(self, pi, min_bins, max_bins, min_reads, avg_read_length, lengthdict):
        #pi is the entire pi, which we can of course get through self.get_pi, but at this stage was already generated
        self.min_bins = min_bins
        self.max_bins = max_bins
        self.avg_read_length = avg_read_length
        self._unique_cov_dic, self._reduced_len_dct = self._get_unique_cov_dic(pi, lengthdict)
    
    @timeit(log_) 
    def _get_unique_cov_dic(self, pi, lengthdic):
        unique_cov_dic = self._get_empty_cov_dic(pi)
        reduced_len_dct = {self._dest_to_did[k]:lengthdic[k] for k in pi}
        for did, rdlst in self._did_to_uniq_rds.iteritems():
            covmap = unique_cov_dic[did]
            for _, pos1, pos2 in rdlst:
                _add_read(covmap, pos1, pos2, self.avg_read_length, 1, reduced_len_dct[did])
        return unique_cov_dic, reduced_len_dct
    
    def recalc_cov_bins(self):
        dids = self._unique_cov_dic.keys()
        for did in dids:
            new_bin_num = self._define_number_of_bins(self.curcovdic[did][1:-1].sum())
            if new_bin_num != len(self._unique_cov_dic[did]):
                covmap = np.zeros(new_bin_num)
                dlen = self._reduced_len_dct[did]
                if did in self._did_to_uniq_rds: #might not have any...
                    for _, pos1, pos2 in self._did_to_uniq_rds[did]:
                        _add_read(covmap, pos1, pos2, self.avg_read_length, 1, dlen)
                self._unique_cov_dic[did] = covmap
    
    def _get_empty_cov_dic(self, pi):
        ret = dict()
        for k_, numreads_ in pi.iteritems():
            numbins = self._define_number_of_bins(numreads_)
            ret[self._dest_to_did[k_]] = np.zeros(numbins)
        return ret
    
    def _define_number_of_bins(self, num_reads):
        return max(self.min_bins, min(self.max_bins,int(num_reads/100))) + 2
    
    def _get_cur_cov_dic(self, delta = None):
        cov_dic = {k:v.copy() for k,v in self._unique_cov_dic.iteritems()}
        if delta is None:
            for mpings in self._ambig_rds.itervalues():
                qualsum = sum(m[r_QUALITY] for m in mpings)
                for m in mpings:
                    _add_read(cov_dic[m[_r_DEST]], m[r_POS1], m[r_POS2], self.avg_read_length, \
                              m[r_QUALITY] / qualsum, self._reduced_len_dct[m[_r_DEST]])
        else:
            for mpings in delta.itervalues():
                for (did, pos1, pos2), koef in mpings.iteritems():
                    _add_read(cov_dic[did], pos1, pos2, self.avg_read_length, koef, self._reduced_len_dct[did])
        return cov_dic
    
    def get_dense_region_coverage_pi(self, perc, delta = None):
        newpi = {}
        self.curcovdic = self._get_cur_cov_dic(delta)
        for did, curcovmap in self.curcovdic.iteritems():
            arr = curcovmap[1:-1]
            arr.sort()
            region_size = int(np.ceil(len(arr)* perc / 100.0)) 
            st = np.argmin(arr[region_size-1:] - arr[:-region_size+1])
            newpi[self._did_to_dest[did]] = np.median(arr[st:st+region_size]) * len(arr)
        return newpi

    def do_delta(self, pi, theta1, epsilon, use_theta):
        METHOD_OLD = False
        delta = dict() 
        if use_theta:
            for did in self._did_to_dest.iterkeys():
                if did not in theta1:
                    theta1[did] = theta1['_Minimum']
        pi = {self._dest_to_did[k]:v for k,v in pi.iteritems()}
        for rrid, mpngs in self._ambig_rds.iteritems():
            if use_theta:
                noms = {(did, pos1, pos2): pi[did] * theta1[did] * qual for qual, pos1, pos2, did in mpngs}
            else:
                noms = {(did, pos1, pos2): pi[did] * qual for qual, pos1, pos2, did in mpngs}
            if METHOD_OLD:
                noms = {k:v for k,v in noms.iteritems() if v >= epsilon}
            else:
                denom = sum(noms.itervalues())
                noms = {k:v for k,v in noms.iteritems() if v/denom >= epsilon}
            denom = sum(noms.itervalues()) 
            delta[rrid] = {k:v/denom for k,v in noms.iteritems()}
        return delta
    
    def get_full_delta(self, pi, theta1, epsilon, use_theta):
        delta = [(self._rrid_to_rid[rrid],((self._did_to_dest[did], pos1, pos2, 1, 1),)) \
                for did, uniq_rds in self._did_to_uniq_rds.iteritems() \
                for rrid, pos1, pos2 in uniq_rds]
        if use_theta:
            for did in self._did_to_dest.iterkeys():
                if did not in theta1:
                    theta1[did] = theta1['_Minimum']
        pi = {self._dest_to_did[k]:v for k,v in pi.iteritems()}
        for rrid, mpngs in self._ambig_rds.iteritems():
            if use_theta:
                noms = tuple((self._did_to_dest[did], pos1, pos2, pi[did] * theta1[did] * qual, qual) for qual, pos1, pos2, did in mpngs)
            else:
                noms = tuple((self._did_to_dest[did], pos1, pos2, pi[did] * qual, qual) for qual, pos1, pos2, did in mpngs)
            denom = sum(v[3] for v in noms)
            noms = tuple(v for v in noms if v[3] / denom >= epsilon)
            denom = sum(v[3] for v in noms)
            denom2 = sum(v[4] for v in noms)
            delta.append((self._rrid_to_rid[rrid], tuple((a,b,c,d/denom, e/denom2) for a,b,c,d,e in noms)))
        return delta
            
def _add_read(covmap, pos1, pos2, avg_read_length, koef, length):
    numbins = len(covmap) 
    binlength = (length - (avg_read_length * 2.))/(numbins - 2)
    ind1 = int((pos1 / binlength) - (avg_read_length/binlength) + 1) 
    ind2 = int((pos2 / binlength) - (avg_read_length/binlength) + 1)
    pos1_end = pos1 + avg_read_length
    pos2_end = pos2 + avg_read_length
    ind1_end = int((pos1_end / binlength) - (avg_read_length/binlength) + 1) 
    ind2_end = int((pos2_end / binlength) - (avg_read_length/binlength) + 1)
    
    if (ind1 < 0):
        ind1 = 0
    if (ind2 < 0):
        ind2 = 0
    if (ind1 > numbins - 1):
        ind1 = numbins - 1
    if (ind2 > numbins - 1):
        ind2 = numbins - 1
    if (ind1_end > numbins - 1):
        ind1_end = numbins - 1
    if (ind2_end > numbins - 1):
        ind2_end = numbins - 1
    if (pos2 > -1): #PE
        koef = koef / 2.
        if (ind2 < ind2_end):
            if (ind2 > 0 and ind2 < numbins - 1): 
                covmap[ind2] += koef * ((binlength - ((pos2 - avg_read_length)%binlength))/avg_read_length)
            if (ind2_end < numbins - 1):
                covmap[ind2_end] += koef * (((pos2_end - avg_read_length)%binlength)/avg_read_length)
            for i in range(ind2+1, ind2_end):
                covmap[i] += (koef * binlength / avg_read_length)
        else:
            covmap[ind2] += koef
    if (ind1 < ind1_end):
        if (ind1 > 0 and ind1 < numbins - 1):
            covmap[ind1] += koef * ((binlength - ((pos1 - avg_read_length)%binlength))/avg_read_length)
        if (ind1_end < numbins - 1):
            covmap[ind1_end] += koef * (((pos1_end - avg_read_length)%binlength)/avg_read_length)
        for i in range(ind1+1, ind1_end):
            covmap[i] += (koef * binlength / avg_read_length)
    else:
        covmap[ind1] += koef 
              
def _calc_qual(sr, pmap_map, seqdict):
    if hasattr(pmap_map, 'mismatches'):
        ga1, off1 = _get_opt_gigar_array(seqdict, pmap_map.dest_id, pmap_map.strand, pmap_map.pos, pmap_map.mismatches, sr.seq[pmap_map.end], pmap_map.gigar)
        return 10**gqa.calc_quality(sr.quals[pmap_map.end], ga1), off1, None
    else:
        ga1, off1 = _get_opt_gigar_array(seqdict, pmap_map.dest_id, pmap_map.strand1, pmap_map.pos1, pmap_map.mismatches1, sr.seq[0], pmap_map.gigar1)
        ga2, off2 = _get_opt_gigar_array(seqdict, pmap_map.dest_id, pmap_map.strand2, pmap_map.pos2, pmap_map.mismatches2, sr.seq[1], pmap_map.gigar2)
        return 10**(gqa.calc_quality(sr.quals[0], ga1) + gqa.calc_quality(sr.quals[1], ga2)), off1, off2

def _get_opt_gigar_array(seqdict, dest_id, strand, pos, mismatches, seq, gigar):
    if '>' in gigar:
        lseq = len(seq)
        genesq, start, end, strand = seqdict.get_sequence(dest_id, strand, pos, lseq)
        return gqa.get_opt_gigar_array(seqdict.get_len(dest_id),genesq, start, end,
                                                                 strand, gigar, lseq, seq, mismatches)
    else:
        return gqa.create_gigar_array(gigar), 0