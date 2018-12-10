from SGVFinder import get_sample_map
import logging
import argparse
import ujson

__author__ = 'Eran Segal (eran.segal@weizmann.ac.il), David Zeevi (dave.zeevi@gmail.com), Tal Korem (tal.korem@gmail.com)'
__version__ = '0.1.0'
__date__ = 'Oct 1st 2018'

def _addargs(parser):
    parser.add_argument('deltaf', help = 'The delta file to process')
    parser.add_argument('outf', help = 'Path to processed file')
    parser.add_argument('average_read_length', help = 'Average read length in the original sample')
    parser.add_argument('--x_coverage', help = 'The desired coverage across the genome in units of 100bp reads. This parameter is used to determine bin size: bin_size = rate_param/x_coverage (Default = 0.01)', type=float, default = 0.01)
    parser.add_argument('--rate_param', help = 'The lower limit for the median number of reads per genomic bin. Genomes with coverage lower than rate_param will be discarded from the analysis (default = 10)', type = int, default = 10)
    
if __name__ == '__main__':
    logging.basicConfig(level=logging.INFO, format = '%(asctime)s-%(levelname)s: %(message)s')
    parser = argparse.ArgumentParser()
    _addargs(parser)
    args = parser.parse_args()
    mp = get_sample_map(args.deltaf, args.x_coverage, args.average_read_length, args.rate_param)
    with open(args.outf, 'wb') as of:
        ujson.dump(mp, of, double_precision=100)
