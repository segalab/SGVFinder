import re
from Bio import SeqIO
import os
import logging
from ICRAUtils import _shell_command, _tryrm
log_ = logging.getLogger(__name__)

GEM_BIN = 'gem-mapper'
statslinergx = re.compile('.* -- #(\d+) sequences processed$')

def _dope(input1, input2, output_prefix, indexf, qual_format="offset-33", qual_thresh=26,
         mismatch = 0.04, editdistance = 0.04, min_matched_bases=0.8, bigindellength=15,
         strataafterbest = 0, uniquemapping = False, max_decoded_matches = 20, 
         min_decoded_strata = 1, peeditdistance = 0.3, max_extendable_matches = 20,
         max_extensions_per_match=1, verbose = True, threads = 1, maxinsertsize = 1000):
    uniquemapping = ' --unique-mapping' if uniquemapping else ''
    uniquepairing = ' --unique-pairing' if uniquemapping else ''
    verbose = ' -v' if verbose else ''
    out = _shell_command(('{bin} -I {indexf} -1 {input1} -2 {input2} -o {output_prefix} -q {qual_format} ' +\
                   '--gem-quality-threshold {qual_thresh} -m {mismatch} -e {editdistance} ' +\
                   '--min-matched-bases {min_matched_bases} --max-big-indel-length {bigindellength} ' +\
                   '-s {strataafterbest}{uniquemapping} -d {max_decoded_matches} -D {min_decoded_strata} ' +\
                   '-p -E {peeditdistance} --max-extendable-matches {max_extendable_matches} ' +\
                   '--max-extensions-per-match {max_extensions_per_match}{uniquepairing}{verbose} ' +\
                   '--max-insert-size {maxinsertsize} -T {threads}').\
                   format(bin = GEM_BIN, **locals()), loglevel = logging.DEBUG)
    return not out.endswith('Wrong alignment\n'), int(statslinergx.match(out.splitlines()[-2]).group(1)) 

def _dose(input1, output_prefix, indexf, qual_format="offset-33", qual_thresh=26,
         mismatch = 0.04, editdistance = 0.04, min_matched_bases=0.8, bigindellength=15,
         strataafterbest = 0, uniquemapping = False, max_decoded_matches = 20, 
         min_decoded_strata = 1, verbose = True, threads = 1):
    uniquemapping = ' --unique-mapping' if uniquemapping else ''
    verbose = ' -v' if verbose else ''
    out = _shell_command(('{bin} -I {indexf} -i {input1} -o {output_prefix} -q {qual_format} ' +\
                   '--gem-quality-threshold {qual_thresh} -m {mismatch} -e {editdistance} ' +\
                   '--min-matched-bases {min_matched_bases} --max-big-indel-length {bigindellength} ' +\
                   '-s {strataafterbest}{uniquemapping} -d {max_decoded_matches} -D {min_decoded_strata}' +\
                   '{verbose} -T {threads}').\
                   format(bin = GEM_BIN, **locals()), loglevel = logging.DEBUG)
    return not out.endswith('Wrong alignment\n'), int(statslinergx.match(out.splitlines()[-2]).group(1)) 

def dope_w_badread_removal(input1, input2, output_prefix, **kwargs):
    suc, startnum = _dope(input1, input2, output_prefix, **kwargs) 
    if not suc:
        _bad_reads_removal((startnum / 2) - 1e4, 3e4, (input1, input2), output_prefix, kwargs)
    _tryrm(output_prefix + '_short_1.fastq')
    _tryrm(output_prefix + '_short_2.fastq')

def dose_w_badread_removal(input1, output_prefix, **kwargs):
    suc, startnum = _dose(input1, output_prefix, **kwargs)
    if not suc:
        _bad_reads_removal(startnum - 1e4, 3e4, (input1, ), output_prefix, kwargs)
    _tryrm(output_prefix + '_short.fastq')

def _bad_reads_removal(startnum, scope, fs, outpref, kwargs):
    rdss = _extract_suspected_sequences(outpref, startnum, scope, fs)
        
    ex_reads = set()
    _find_fn_reads(rdss, outpref + '_tmp', range(int(scope)), kwargs, ex_reads)
    with open(fs[0]) as i1, open(outpref + ('_tmp_1.fastq' if len(fs) == 2 else '_tmp.fastq'), 'wb') as o1:
        o1.writelines(l for i, l in enumerate(i1) if (i/4) - startnum not in ex_reads)
    if len(fs) == 2:
        with open(fs[1]) as i2, open(outpref + '_tmp_2.fastq', 'wb') as o2:
            o2.writelines(l for i, l in enumerate(i2) if (i/4) - startnum not in ex_reads)
    log_.info('Removed %s reads' % len(ex_reads))
    _tryrm(outpref + '_short.fastq')
    _tryrm(outpref + '_short_1.fastq')
    _tryrm(outpref + '_short_2.fastq')
    
    if len(fs) == 2:
        os.rename(outpref + '_tmp_1.fastq', outpref + '_short_1.fastq')
        os.rename(outpref + '_tmp_2.fastq', outpref + '_short_2.fastq')
        dope_w_badread_removal(outpref + '_short_1.fastq', outpref + '_short_2.fastq', outpref, **kwargs)
    else:
        os.rename(outpref + '_tmp.fastq', outpref + '_short.fastq')
        dose_w_badread_removal(outpref + '_short.fastq', outpref, **kwargs)
    

def _extract_suspected_sequences(outpref, startnum, scope, fs):
    if len(fs) == 2:
        tmp1 = outpref + '_tmp_1.fastq'
        tmp2 = outpref + '_tmp_2.fastq'
    else:
        tmp1 = outpref + '_tmp.fastq'
        
    _shell_command('head -n {hn} {f} | tail -n {tn} > {of}'.\
                  format(hn = int((startnum + scope) * 4), f = fs[0], \
                         tn = int(scope * 4), of = tmp1))
    if len(fs) == 2:
        _shell_command('head -n {hn} {f} | tail -n {tn} > {of}'.\
                  format(hn = int((startnum + scope) * 4), f = fs[1], \
                         tn = int(scope * 4), of = tmp2))
        
    with open(tmp1) as i1:
        rds1 = [r for r in SeqIO.parse(i1, 'fastq')]
    if len(fs) == 2:
        with open(tmp2) as i2:
            rds2 = [r for r in SeqIO.parse(i2, 'fastq')]
    _tryrm(tmp1)
    if len(fs) == 2:
        _tryrm(tmp2)
    return (rds1, rds2) if len(fs) == 2 else (rds1, )
   
def _find_fn_reads(rdss, o_prefix, pos, kwargs, ex_reads):
    rds1_1 = rdss[0][:len(rdss[0])/2]
    rds1_2 = rdss[0][len(rdss[0])/2:]
    pos_1 = pos[:len(rdss[0])/2]
    pos_2 = pos[len(rdss[0])/2:]
    if len(rdss) == 2:
        rds2_1 = rdss[1][:len(rdss[0])/2]
        rds2_2 = rdss[1][len(rdss[0])/2:]
    
    suc1 = _do_one_half((rds1_1, rds2_1) if len(rdss) == 2 else (rds1_1, ), o_prefix, pos_1, ex_reads, kwargs)
    suc2 = _do_one_half((rds1_2, rds2_2) if len(rdss) == 2 else (rds1_2, ), o_prefix, pos_2, ex_reads, kwargs)
        
    if suc1 and suc2:
        for ind in pos:
            ex_reads.add(ind)
        log_.debug('Removed %s reads' % (len(pos)))
    
def _do_one_half(rdss, o_prefix, pos, ex_reads, kwargs):
    if len(rdss[0]) == 0:
        return True
    else:
        ns = _write_reads(rdss, o_prefix)
        if len(ns) == 2:
            suc, _ = _dope(ns[0], ns[1], o_prefix, **kwargs)
        else:
            suc, _ = _dose(ns[0], o_prefix, **kwargs)
        if not suc:
            if len(rdss[0]) == 1:
                ex_reads.add(pos[0])
                log_.debug('Removed one read')
            else:
                _find_fn_reads(rdss, o_prefix, pos, kwargs, ex_reads)
        _tryrm(ns[0])
        if len(ns) == 2:
            _tryrm(ns[1])
        _tryrm(o_prefix + '.map')
        return suc
    
def _write_reads(rdss, o_prefix):
    if len(rdss) == 2:
        n1 = o_prefix + '_1.fastq'
        n2 = o_prefix + '_2.fastq'
    else:
        n1 = o_prefix + '.fastq'
        
    with open(n1, 'wb') as o1:
        SeqIO.write(rdss[0], o1, 'fastq')
    if len(rdss) == 2:
        with open(n2, 'wb') as o2:
            SeqIO.write(rdss[1], o2, 'fastq')
    return (n1, n2) if len(rdss) == 2 else (n1, )
    
def validate_success_ESS(mfile, fqfile):
    mpl = _shell_command('tail -n 1 %s' % mfile)
    fql = _shell_command('tail -n 4 %s | head -n 1' % fqfile)
    return mpl.split()[0] == fql.split()[0][1:], mpl.split()[0], fql.split()[0][1:]  

