from scipy.stats._continuous_distns import norm, betaprime, ncx2
from itertools import combinations, product
from scipy.spatial.distance import pdist, squareform
from scipy.cluster.hierarchy import linkage, fcluster
from scipy.stats.stats import spearmanr
from os.path import exists
from pandas import read_pickle, np, to_pickle, DataFrame, concat,Series
from os.path import basename, join, splitext
import logging
from collections import defaultdict
from glob import glob
import ujson
import bokeh.plotting as bpl 
from os.path import split, realpath
from bokeh.layouts import column
from bokeh.models.sources import ColumnDataSource
from bokeh.models.tools import HoverTool, TapTool
from bokeh.models.ranges import Range1d
from bokeh.models.callbacks import OpenURL
log_ = logging.getLogger('SGVF')

def get_sample_map(delta_fname, x_coverage, average_read_length, rate_param):
    lengthdb = read_pickle(lengthdbpath)
    bin_size = int(rate_param / float(x_coverage))
    with open(delta_fname) as inf:
        delta = ujson.load(inf, precise_float = True)
    bacid_maps = dict()
    for _, mapngs in delta:
        for dest_id, pos1, pos2, used_koef, _ in mapngs:
            if dest_id not in bacid_maps:
                bacid_maps[dest_id] = np.zeros(int(lengthdb[dest_id] / bin_size) + 1)
            ind1 = int((int(pos1) + (int(average_read_length) / 2)) / bin_size)
            if pos2 >= 0:
                used_koef = used_koef / 2.0
                ind2 = int((int(pos2) + (int(average_read_length) / 2)) / bin_size)
                bacid_maps[dest_id][ind2] += used_koef
            bacid_maps[dest_id][ind1] += used_koef
    return {dest_id:cov_map for dest_id, cov_map in bacid_maps.iteritems()\
                           if np.median(cov_map) >= rate_param}

def calculate_by_other(old_deldf_p, old_sgvdf_p, old_frames_df, samp_to_map, \
                       real_del_thresh, dense_perc, weboutputdir, x_coverage, rate_param):
    binsize = int(rate_param / float(x_coverage))
    dichotomize = True
    dichotomize_thresh = 0.5
    old_deldf = read_pickle(old_deldf_p)
    old_sgvdf = read_pickle(old_sgvdf_p)
    
    bac_samps_map = defaultdict(dict)
    for samp, bacid_maps in samp_to_map.iteritems():
        for bacname, bacmap in bacid_maps.iteritems():
            bac_samps_map[bacname][samp] = bacmap
    
    bac_bacdf = dict()
    for bacname, bacdict in bac_samps_map.iteritems():
        bacdf = DataFrame(bacdict).T
        bac_bacdf[bacname] = bacdf
    
    matforweb = dict()        
    delsregions, covregions = [], []
    
    prev_bac = None
    for col in old_deldf.columns:
        cur_bac = col.split(':')[0]
        if not cur_bac in bac_bacdf:
            continue
        if cur_bac != prev_bac:
            if prev_bac is not None:
                matforweb[prev_bac] = [DataFrame(), [], cur_deltn, bac_del_regions]
            t_dip = bac_bacdf[cur_bac].apply(finddip, axis = 1)
            cur_deltn = (~bac_bacdf[cur_bac].le(t_dip, axis = 0)).astype(int)
            prev_bac = cur_bac
            bac_del_regions = []
        if dichotomize:
            delsregions.append((extract_scluster(cur_deltn,[(int(st.split('_')[0]), int(st.split('_')[1])) \
                                                            for st in col.split(':')[1].split(';')])\
                                .mean(1) > dichotomize_thresh).astype(int).to_frame(col))
        else:
            delsregions.append(extract_scluster(cur_deltn,[(int(st.split('_')[0]), int(st.split('_')[1])) \
                                                            for st in col.split(':')[1].split(';')]).\
                               mean(1).to_frame(col))
        bac_del_regions.append(delsregions[-1])
    if prev_bac is not None:
        matforweb[cur_bac] = [DataFrame(), [], cur_deltn, bac_del_regions]
    prev_bac = None
    cache = dict()
    bac_cov_regions = []
    for col in old_sgvdf.columns:
        cur_bac = col.split(':')[0]
        if not cur_bac in bac_bacdf:
            continue
        if cur_bac != prev_bac:
            if prev_bac is not None:
                if prev_bac not in matforweb:
                    matforweb[prev_bac] = [DataFrame(), [], DataFrame(), []]
                dx, _, _ = vsgv_get_tvar_params(1, cur_nodeldf, dense_perc, .02, None, None, cache)
                matforweb[prev_bac][0] = cur_nodeldf.subtract(dx['mean'], axis = 0)\
                                                   .truediv(dx['std'], axis = 0)
                matforweb[prev_bac][1] = bac_cov_regions
            cache = dict()
            if exists(join(old_frames_df, cur_bac + '.df')):
                old_bacdf = read_pickle(join(old_frames_df, cur_bac + '.df'))
                old_dip = old_bacdf.apply(finddip, axis = 1)
                old_deltn = (~old_bacdf.le(old_dip, axis = 0)).astype(int) 
                cur_nodeldf = bac_bacdf[cur_bac].loc[:,(old_deltn ==0).sum() < old_deltn.shape[0]*real_del_thresh]
                prev_bac = cur_bac
                bac_cov_regions = []
            else:
                log_.error('Missing bac frame at ' + join(old_frames_df, cur_bac + '.df'))
                continue
        covregions.append(vsgv_get_region_val([(int(st.split('_')[0]), int(st.split('_')[1])) for st in col.split(':')[1].split(';')], \
                                              cur_nodeldf, dense_perc, .02, None, None, cache).to_frame(col))
        bac_cov_regions.append(covregions[-1])
    if prev_bac is not None:
        if prev_bac not in matforweb:
            matforweb[prev_bac] = [DataFrame(), [], DataFrame(), []]
        dx, _, _ = vsgv_get_tvar_params(1, cur_nodeldf, dense_perc, .02, None, None, cache)
        matforweb[prev_bac][0] = cur_nodeldf.subtract(dx['mean'], axis = 0)\
                                                   .truediv(dx['std'], axis = 0)
        matforweb[prev_bac][1] = bac_cov_regions
    if weboutputdir is not None:
        taxonomy = read_pickle(taxonomypath)
        genepos = read_pickle(genepospath)
        for bacname, (normdf, sgvregions, deldf, delsregions_) in matforweb.iteritems():
            draw_one_region(bacname, binsize, taxonomy, normdf, \
                            deldf, delsregions_, bacdf, sgvregions, weboutputdir, genepos)
    return concat(covregions, axis = 1) if len(covregions) > 0 else DataFrame(), \
        concat(delsregions, axis = 1) if len(delsregions) > 0 else DataFrame()

def collect_from_glob(delta_glob, q, dist_reads, x_coverage, average_read_length, lengthdbp, rate_param, cutoff, outp):
    collect_from_flist(glob(delta_glob), q, dist_reads, x_coverage, average_read_length, lengthdbp, rate_param, cutoff, outp)

def collect_from_flist(deltaflist, q, dist_reads, x_coverage, average_read_length, lengthdbp, rate_param, cutoff, outp):
    samp_to_map = {splitext(basename(f))[0]: \
                   q.method(get_sample_map, (f, dist_reads, x_coverage, \
                                  average_read_length, lengthdbp, \
                                  rate_param, cutoff)) for f in deltaflist}
    samp_to_map = {k:q.waitforresult(v) for k,v in samp_to_map.iteritems()}
    to_pickle(samp_to_map, outp)

def work_on_collection(samp_to_map, min_samp_cutoff, delsdetectthresh, real_del_thresh, dels_cooc_thresh, 
                       vsgv_dissim_thresh, vsgv_clip_quantile, vsgv_fit_interval, vsgv_fit_method, 
                   x_coverage, rate_param, vsgv_dense_perc, browser_path):
    binsize = int(rate_param / float(x_coverage))
    dichotomize = True
    dichotomize_thresh = 0.5
    max_spacing = 10
    taxonomy = read_pickle(taxonomypath)
    genepos = read_pickle(genepospath)
    bac_samps_map = defaultdict(dict)
    for samp, bacid_maps in samp_to_map.iteritems():
        for bacname, bacmap in bacid_maps.iteritems():
                bac_samps_map[bacname][samp] = bacmap
    
    sgvregions_all = []
    delsregions_all = [] 
    for bacname, bacdict in bac_samps_map.iteritems():
        bacdf = DataFrame(bacdict).T
        if bacdf.shape[0] < min_samp_cutoff: continue
        if (bacdf.median() < 1).sum() / float(bacdf.shape[1]) > 0.3:
            continue
        delsregions, deldf = find_deletions(bacdf, bacname, dichotomize, dichotomize_thresh, \
                                            delsdetectthresh, max_spacing, dels_cooc_thresh)
        delsregions_all.extend(delsregions)
        sgvregions, normdf = find_sgvs(bacdf, max_spacing, vsgv_dense_perc, bacname, deldf, \
                               real_del_thresh, vsgv_clip_quantile, vsgv_fit_interval, \
                               vsgv_fit_method, vsgv_dissim_thresh)
        sgvregions_all.extend(sgvregions)
        
        if browser_path is not None:
            draw_one_region(bacname, binsize, taxonomy, normdf, \
                            deldf, delsregions, bacdf, sgvregions, browser_path, genepos)
    return concat(sgvregions_all, axis = 1) if len(sgvregions_all) > 0 else DataFrame(), \
        concat(delsregions_all, axis = 1) if len(delsregions_all) > 0 else DataFrame()
        
def dense_stats(a, perc):
    x = a.copy().values[1:-1]
    x.sort()
    region_size = int(np.ceil(len(x)* perc / 100.0)) 
    st = np.argmin(x[region_size-1:] - x[:-region_size+1])
    x = x[st:st+region_size]
    return x.mean(), x.std()

def vsgv_can_extend(st, en, nodeldf, perc, clip_q, fit_met, fit_int, cache):
    slen = en-st
    dx, qs, cutoff = vsgv_get_tvar_params(slen, nodeldf, perc, clip_q, fit_met, fit_int, cache)
    bt = nodeldf.loc[:, st:en-1].sum(1).subtract(dx['mean'], axis = 0)\
        .truediv(dx['std'], axis = 0).clip(*[q.loc[st:en-1].mean() for q in qs])
    if fit_met == 'betaprime':
        return bt.std() > cutoff
    elif fit_met == 'ncx2':
        return bt.var() > cutoff
    else:
        raise ValueError

def vsgv_get_tvar_params(slen, nodeldf, perc, clip_q, fit_met, fit_int, cache):
    if slen not in cache:
        ndf = nodeldf.rolling(slen, slen, axis = 1).sum().iloc[:, (slen-1)::slen] if slen > 1 else nodeldf
        dx = ndf.apply(dense_stats, args = (perc,), axis = 1).apply(Series).rename(columns={0:'mean', 1:'std'})
        nnormdf =  ndf.subtract(dx['mean'], axis = 0)\
                    .truediv(dx['std'], axis = 0)
        qs = nnormdf.quantile(clip_q), nnormdf.quantile(1-clip_q)
        if fit_met == 'betaprime':
            ncv = nnormdf.clip(*qs, axis = 1).std()
            cutoff = betaprime(*betaprime.fit(ncv)).interval(fit_int)[1]
        elif fit_met == 'ncx2':
            ncv = nnormdf.clip(*qs, axis = 1).var()
            cutoff = ncx2(*ncx2.fit(ncv)).interval(fit_int)[1]
        else:
            ncv, cutoff = None, None
        cache[slen] = dx, qs, cutoff
    return cache[slen]

def vsgv_get_region_val(scluster, nodeldf, perc, clip_q, fit_met, fit_int, cache):
    usclust = extract_scluster(nodeldf, scluster)
    succ = False
    s = usclust.shape[1]
    while not succ:
        dx, _, _ = vsgv_get_tvar_params(s, nodeldf, perc, clip_q, fit_met, fit_int, cache)
        succ = True
    return usclust.sum(1).subtract(dx['mean'], axis = 0)\
        .truediv(dx['std'], axis = 0)
        
def find_sgvs(bacdf, max_spacing, dense_perc, bacname, deldf, real_del_thresh, clip_quantile, \
              fit_interval, fit_method, dissim_thresh):
    assert fit_method in {'betaprime', 'ncx2'}
    nodeldf = bacdf.loc[:,(deldf ==0).sum() < deldf.shape[0]*real_del_thresh]
    cache = dict()
    dx, qs, cutoff = vsgv_get_tvar_params(1, nodeldf, dense_perc, clip_quantile, fit_method, fit_interval, cache)
    normdf = nodeldf.subtract(dx['mean'], axis = 0)\
        .truediv(dx['std'], axis = 0)
    if normdf.isnull().sum().sum() > 0:
        normdf = normdf.replace(np.inf, np.nan).replace(-np.inf, np.nan).dropna(how = 'all') ##needed?
    if fit_method == 'betaprime':
        dfcv = normdf.clip(*qs, axis = 1).std()
    elif fit_method == 'ncx2':
        dfcv = normdf.clip(*qs, axis = 1).var()
    else:
        raise ValueError
    
    seeds = (dfcv > cutoff).astype(int)#.diff()
    def ex_func(s0,s1,snew0,snew1):
        sX = nodeldf.loc[:,s0:s1-1].T.values
        sY = nodeldf.loc[:,snew0:snew1-1].T.values
        return  (np.average([_spearman_dissim(sx,sy) for sx, sy in product(sX, sY)]) <= dissim_thresh) \
                and vsgv_can_extend(min([s0,s1,snew0,snew1]), max(s0,s1,snew0,snew1),  
                                    nodeldf, dense_perc, clip_quantile, fit_method, fit_interval, cache)
    stretches = get_extended_stretches(seeds, bacdf.shape[1], max_spacing, ex_func)
    sclusters = cluster_stretches(stretches, nodeldf, _spearman_dissim, dissim_thresh) 
    covregions = []
    for s in sclusters:
        covregions.append(vsgv_get_region_val(s, nodeldf, \
                                              dense_perc, clip_quantile, fit_method, fit_interval, cache)\
                          .to_frame("{}:{}".format(bacname, scluster_to_txt(s))))
    return covregions, normdf

def finddip(x):
    min_n = 1e9
    prev_b = 1
    for n,b in zip(*np.histogram(x.values, bins = 11, range = (0,11))):
        if b >= 10:
            return 1
        if n > min_n:
            return prev_b
        min_n = n
        prev_b = b

def find_deletions(df_in, bacname, dichotomize, dichotomize_thresh, detectthresh, 
                   max_spacing, cooc_thresh):
    cutoffvar = detectthresh*(1-detectthresh)
    df = df_in.copy()
    dip = df.apply(finddip, axis = 1)
    deltn = (~df.le(dip, axis = 0)).astype(int)
    def getpercdel(s0, s1):
        return (~deltn.iloc[:,s0:s1].any(1)).sum() / float(len(deltn))
    def getpercret(s0,s1):
        return deltn.iloc[:, s0:s1].all(1).sum() / float(len(deltn))
    def ext_cret(s0,s1,snew0,snew1):
        sX = deltn.iloc[:,s0:s1].T.values
        sY = deltn.iloc[:,snew0:snew1].T.values
        return (np.average([_cooc_dissim(sx,sy) for sx, sy in product(sX, sY)]) <= cooc_thresh) \
                and getpercdel(s0, s1) > detectthresh and getpercret(s0, s1) > detectthresh
    seeds = (deltn.var(0) > cutoffvar).astype(int)
    delsregions = []
    stretches = get_extended_stretches(seeds, df.shape[1], max_spacing, ext_cret)
    sclusters = cluster_stretches(stretches, deltn, _cooc_dissim, cooc_thresh) 
    for s in sorted(sclusters, key = lambda x:x[0][0]): 
        if dichotomize:
            delsregions.append((extract_scluster(deltn,s).mean(1) > dichotomize_thresh).\
                                    astype(int).to_frame("{}:{}".format(bacname, scluster_to_txt(s))))
        else:
            delsregions.append(extract_scluster(deltn,s).mean(1).to_frame("{}:{}"\
                                                                        .format(bacname, scluster_to_txt(s))))
    return delsregions, deltn
    
def scluster_to_txt(scluster):
    return ';'.join(['{}_{}'.format(*s) for s in scluster])

def extract_scluster(df, scluster):
    return concat([df.loc[:,s[0]:s[1]-1] for s in scluster], axis=1)
    
def get_extended_stretches(seeds, dfs1, max_spacing, extend_func):
    stretches = [(s, s+1) for s in seeds[seeds == 1].index]
    if len(stretches) == 0:
        return stretches
    while True:
        newstretches = []
        for i,s in enumerate(stretches):
            while (i == 0 or s[0] > newstretches[i-1][1]) and s[0] > 0 and extend_func(s[0], s[1], s[0]-1, s[0]):
                s = (s[0]-1, s[1])
            while (i == len(stretches) - 1 or s[1] < stretches[i+1][0]) and s[1] < dfs1 -1 \
                    and extend_func(s[0], s[1], s[1]+1, s[1]+2):
                s = (s[0], s[1]+1)
            newstretches.append(s)
        newstretches2 = [newstretches[0]]
        for s2 in newstretches[1:]:
            if s2[0] - newstretches2[-1][1] <= max_spacing \
                    and extend_func(newstretches2[-1][0],s2[0],s2[0],s2[1]) \
                    and extend_func(newstretches2[-1][0],newstretches2[-1][1],newstretches2[-1][1],s2[1]):
                newstretches2[-1] = (newstretches2[-1][0], s2[1])
            else:
                newstretches2.append(s2)
        if newstretches2 == stretches:
            break
        else:
            stretches = newstretches2
    return stretches

def _stretch_average(distance, stretches, indices):
    distmat = DataFrame(squareform(distance), index = indices, columns = indices)
    return np.array([np.average(np.ravel(distmat.loc[si[0]:si[1]-1,sj[0]:sj[1]-1])) for si,sj in combinations(stretches,2)])
    
def cluster_stretches(stretches, df, unite_func, cutoff_dist):
    if len(stretches) == 0:
        return []
    if len(stretches) == 1:
        return [(stretches[0],)]
    linkage_method='average' # only average is implemented here
    stretchdf = concat([df.loc[:,s[0]:s[1]-1] for s in stretches], axis = 1)
    distance = pdist(stretchdf.T, unite_func)
    distance = _stretch_average(distance, stretches, stretchdf.columns)
    Z = linkage(distance, method = linkage_method)
    clusterdf = DataFrame({'cluster':fcluster(Z, cutoff_dist, criterion='distance')}, index = stretches)
    return clusterdf.reset_index().sort_values(['cluster','index']).groupby('cluster')['index'].apply(tuple).values
            
def _cooc_dissim(v,u):
    return 1-((v==u).sum()/float(len(u)))

def _spearman_dissim(v,u):
    return 1-((spearmanr(v,u)[0]+1)/2) 

def draw_one_region(bacname, binsize, taxonomy, normdf,
                    deldf, delregions, bacdf, sgvregions, outputdir, geneposs):
    tax = taxonomy.ix[int(bacname.split('.')[0])]['organism']
    if str(tax) == 'nan': tax = bacname
    p1 = bpl.figure(plot_width=1200, plot_height = 400, 
                    title = "{}: {} people, {} variance regions. ".format(tax, bacdf.shape[0], len(sgvregions)))
    p3 = bpl.figure(plot_width=1200, plot_height = 250, x_range = p1.x_range, 
                title = "{} deletion regions".format(len(delregions)))
    # Variable regions
    p1.grid.grid_line_alpha=0.3
    p1.yaxis.axis_label = 'Standardized coverage'
    dfmin = normdf.min().min()
    dfmax = normdf.max().max()
    if len(sgvregions) > 0:
        for clust in concat(sgvregions, axis = 1):
            for reg in clust.split(':')[1].split(';'):
                s = [int(si) for si in reg.split('_')] 
                p1.vbar(s[0] + (s[1]-s[0])/2. -.5, s[1]-s[0]-.5, dfmax, dfmin, 
                                color = 'DarkBlue', alpha = 0.2)
    drawdf = DataFrame({i:normdf[i] if i in normdf.columns else np.nan for i in bacdf.columns})
    p1.line(range(drawdf.shape[1]), drawdf.quantile(0.01), color='Black', alpha = .7)
    p1.line(range(drawdf.shape[1]), drawdf.quantile(0.25), color='Black', alpha = 1)
    p1.line(range(drawdf.shape[1]), drawdf.median(), color='Black', alpha = 1, line_width = 2.5)
    p1.line(range(drawdf.shape[1]), drawdf.quantile(0.75), color='Black', alpha = 1)
    p1.line(range(drawdf.shape[1]), drawdf.quantile(0.99), color='Black', alpha = .7)
    locgeneposs = geneposs.ix[bacname]
    locgeneposs = locgeneposs[~locgeneposs.isnull().any(1)]
    noise = norm.rvs(0,0.75,size = len(locgeneposs))
    source = ColumnDataSource(data = dict(
                x=locgeneposs[['start_pos','end_pos']].mean(1).truediv(binsize).values,
                xs=zip(locgeneposs['start_pos'].truediv(binsize).values, locgeneposs['end_pos'].truediv(binsize).values),
                y=locgeneposs['strand'].replace({'+':2, '-':-2}),
                ys = zip(locgeneposs['strand'].replace({'+':2, '-':-2}).add(noise),locgeneposs['strand'].replace({'+':2, '-':-2}).add(noise)),
                start=locgeneposs.start_pos.apply(int).apply(str).values,
                end = locgeneposs.end_pos.apply(int).apply(str).values,
                strand = locgeneposs.strand.values,
                gene_name = locgeneposs.index.values,
                product = locgeneposs['product'].values,
                gene_type = locgeneposs['gene_type'].values,
                color = locgeneposs['gene_type'].apply(lambda x:'DarkRed' if x == 'CDS' else 'DodgerBlue' if x == 'rRNA' else 'Yellow' if x == 'tRNA' else 'Gray').values))
    hover = HoverTool(tooltips=[
                                ("Gene name", "@gene_name"),
                                ("Position", "@start - @end (@strand)"),
                                ("Gene type", "@gene_type"),
                                ("Product", "@product"),
                                ],
                      mode = 'vline')
    p2 = bpl.figure(plot_width=1200, plot_height = 150, x_range = p1.x_range, 
                    title = None, tools = [hover, 'tap'])
    p2.multi_line('xs','ys', line_width = 4, color = 'color', source=source)
    url = "http://www.google.com/search?q=@gene_name"
    taptool = p2.select(type=TapTool)
    taptool.callback = OpenURL(url=url)
    p2.yaxis.axis_label = 'Genes'
    p2.y_range = Range1d(-5, 5)
    # Deletions
    if p3 is not None:
        p3.xaxis.axis_label = 'Bin position'
        p3.yaxis.axis_label = 'Ratio deleted'
        mednormdf = bacdf.truediv(bacdf.median(1), axis = 0)
        mednormdf[mednormdf > 10] = 10
        if len(delregions) > 0:
            for clust in concat(delregions, axis=1):
                for reg in clust.split(':')[1].split(';'):
                    d = [int(di) for di in reg.split('_')] 
                    p3.vbar(d[0] + (d[1]-d[0])/2. -.5, d[1]-d[0]-.5, 1.5, 0, 
                                    color = 'ForestGreen', alpha = 0.5)
        a = 1- deldf.sum(0).truediv(deldf.count(0))
        p3.line(range(len(a)), a.values, color='FireBrick', alpha = 1, line_width = 2.5)
        # Output
        bpl.output_file(join(outputdir, bacname + '.html'), 
                title=tax)
        if p2 is None:
            bpl.save(column(p1, p3))
        else:
            bpl.save(column(p1, p2, p3))
    else:
        if p2 is None:
            bpl.save(column(p1, ))
        else:
            bpl.save(column(p1, p2))

lengthdbpath = join(split(realpath(__file__))[0], '../DataFiles/representatives.contigs.drepped.dlen')
taxonomypath = join(split(realpath(__file__))[0], '../DataFiles/representatives.genomes.taxonomy.df')
genepospath = join(split(realpath(__file__))[0], '../DataFiles/representatives.genes.drepped.annotations.df')