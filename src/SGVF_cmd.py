import logging
import argparse
from os.path import splitext, basename, split, realpath, join
import ujson
from SGVFinder import work_on_collection, calculate_by_other
from glob import glob

BETAPRIME = 'betaprime'
CHISQ = 'ncx2'

def _addargs(parser):
    parser.add_argument('input_glob', help = 'A glob string that includes all the files processed with the PerFile_cmd')
    parser.add_argument('output_dsgv', help = 'Output path for the deletion-sgv dataframe. By default a pickled pandas dataframe')
    parser.add_argument('output_vsgv', help = 'Output path for the variable-sgv dataframe. By default a pickled pandas dataframe')
    parser.add_argument('--x_coverage', help = 'The desired coverage across the genome in units of 100bp reads. This parameter is used to determine bin size: bin_size = rate_param/x_coverage (Default = 0.1)', type=float, default = 0.01)
    parser.add_argument('--rate_param', help = 'The lower limit for the median number of reads per genomic bin. Genomes with coverage lower than rate_param will be discarded from the analysis (Default = 10)', type = int, default = 10)
    parser.add_argument('--byorig', help = 'Calculate SGVs according to the regions published by Zeevi et al, XXX, 201Y.', action = 'store_true')
    parser.add_argument('--min_samp_cutoff', help = 'Minimum number of samples in which a microbe exists with sufficient coverage to be considered in the analysis (Default=75)', type=int, default = 75)
    parser.add_argument('--dels_detect_thresh', help = 'Determines the minimum and maximum ratio of samples for which a bin is considered a deletion-SGV. (Default=0.25, setting the ratio at 0.25-0.75)', type=float, default = 0.25)
    parser.add_argument('--real_del_thresh', help = 'Threshold above which a bin is considered deleted for all individuals (Default=0.95)', type = float, default = 0.95)
    parser.add_argument('--vsgv_dissim_thresh', help = 'Maximal correlation dissimilarity for concatenation and clustering of variable SGV bins. Correlation dissimilarity is defined as calculated as 1-((rho(u,v)+1)/2), where rho is the Spearman correlation and u, v are the bin vectors being compared (Default=0.125)', type = float, default = 0.125)
    parser.add_argument('--dels_cooc_thresh', help = 'Maximal cooccurrence dissimilarity for concatenation and clustering of deletion SGV bins. Coocurrence dissimilarity is defined as the proportion of samples which are in disagreement on the deletion-state of the two bins being compared (wherein one bin is deleted and one is retained for the same sample) out of all samples that harbor the microbe (Default=0.25)', type = float, default = 0.25) 
    parser.add_argument('--vsgv_clip_quantile', help = 'Determines clipping performed on the distribution of bin values prior to fitting a distribution for the detection of variable-SGVs (Default=0.02 corresponding to clipping outside the 2nd to 98th percentiles)', type = float, default = 0.02)
    parser.add_argument('--vsgv_fit_interval', help = 'Significance cutoff for the fitted distribution above which a bin is considered variable (Default=0.95)', type = float, default = 0.95)
    parser.add_argument('--vsgv_fit_method', help = 'Determines the distribution being fit on bin values (either a Beta-prime or Chi-square distribution; Default=betaprime)', type = str, default = BETAPRIME, choices = [BETAPRIME, CHISQ])
    parser.add_argument('--vsgv_dense_perc', help = 'The percent of the data that is considered when standardizing the bin values of a microbe in a sample. The algorithm chooses the densest part of the data. If a percentage p is selected, the algorithm calculates a subset x of the vector of bins such that max(x)-min(x) is minimal and |x| = p*length(bins). The mean and standard deviation of this vector are calcuated and used to standardize the bin vector (Default=85)', type = float, default = 85)
    parser.add_argument('--browser_path', help = 'Optional; a path for the html output of SGV-Browser (Default=None, resulting in no output)', type=str, default = None)
    parser.add_argument('--csv_output', action = 'store_true', help = 'Will output a csv instead of pandas dataframe')
def _load_ujsn(fname):
    with open(fname) as inf:
        return ujson.load(inf, precise_float = True)

if __name__ == '__main__':
    logging.basicConfig(level=logging.INFO, format = '%(asctime)s-%(levelname)s: %(message)s')
    parser = argparse.ArgumentParser()
    _addargs(parser)
    args = parser.parse_args()
    samp_to_map = {splitext(basename(f))[0]: _load_ujsn(f) for f in glob(args.input_glob)}
    if args.byorig:
        vsgv, dsgv = calculate_by_other(join(split(realpath(__file__))[0], '../DataFiles/orig_dsgv.df'),
                           join(split(realpath(__file__))[0], '../DataFiles/orig_vsgv.df'),
                           join(split(realpath(__file__))[0], '../DataFiles/orig_frames'), 
                           samp_to_map, args.real_del_thresh, args.vsgv_dense_perc, args.browser_path, args.x_coverage, args.rate_param)
    else:
        vsgv, dsgv = work_on_collection(samp_to_map, args.min_samp_cutoff, args.dels_detect_thresh, args.real_del_thresh, \
                       args.dels_cooc_thresh, args.vsgv_dissim_thresh, args.vsgv_clip_quantile, args.vsgv_fit_interval, \
                       args.vsgv_fit_method, args.x_coverage, args.rate_param, args.vsgv_dense_perc, args.browser_path)
    if args.csv_output:
        vsgv.to_csv(args.output_vsgv)
        dsgv.to_csv(args.output_dsgv)
    else:
        vsgv.to_pickle(args.output_vsgv)
        dsgv.to_pickle(args.output_dsgv)
