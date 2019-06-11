from ICRA import CAMI, GENOMES, single_file
import logging
import argparse

__author__ = 'Eran Segal (eran.segal@weizmann.ac.il), David Zeevi (dave.zeevi@gmail.com), Tal Korem (tal.korem@gmail.com)'
__version__ = '0.1.0'
__date__ = 'Oct 1st 2018'

    
def _addargs(parser):
    parser.add_argument('outfol', help = 'Where to write results')
    parser.add_argument('fq', help = 'Name of the fastq file, without the extension. For PE file, _1.fastq or _2.fastq will be appended to what you supply here. For SE, .fastq will')
    parser.add_argument('--pe', help = 'Whether paired end or single end', action = 'store_true')
    parser.add_argument('--max_mismatch', help = 'How many mismatch are considered acceptable', default = 8, type = int)
    
    parser.add_argument('--ignore_lengths', help = 'Should genome lengths be ignored when calculating abundances', default = True, action = 'store_false' )
    
    parser.add_argument('--epsilon', help = 'The stop criteria. Epsilon is the euclidean distance between the internal vectors of element abundances under which ICRA stops.', default = 1e-6, type = float)
    
    parser.add_argument('--max_iterations', help = 'An upper limit to the number of iterations ICRA will run', default = 100, type = int)
    
    parser.add_argument('--min_bins', help = 'The minimum number of bins per genomic element.', default = 10, type = int)
    
    parser.add_argument('--max_bins', help = 'The maximum number of bins per genomic element', default = 100, type = int)
    
    parser.add_argument('--min_reads', help = 'Minimal number of reads mapped to a genomic element for it to be considered present in the sample.', default = 100, type = int)
    
    parser.add_argument('--dense_region_coverage', help = 'The percentage of the genome examined for coverage purposes.', default = 60, type = int)
    
    parser.add_argument('--length_minimum', help = 'Minimal genome length considered,' , default = 1e5, type = float)
    
    parser.add_argument('--length_maximum', help = 'Maximal genome length considered.', default = 2e7, type = float)
    
    parser.add_argument('--usage', help = 'Whether ICRA is mapping to the default db (genomes) or to the db used in the paper for CAMI (cami)', default = GENOMES, type = str, choices = [CAMI, GENOMES])
    
    parser.add_argument('--use_theta', help = 'Theta is an extra algorithmic components that causes ICRA to converge faster', action = 'store_true')
    
    parser.add_argument('--debug', help = 'Logs debug output', action = 'store_true')
                  
if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    _addargs(parser)
    args = parser.parse_args()
    logging.basicConfig(level=(logging.DEBUG if args.debug else logging.INFO), format = '%(asctime)s-%(levelname)s: %(message)s')
    single_file(args.fq + ('_1.fastq' if args.pe else '.fastq'), (args.fq + '_2.fastq') if args.pe else None, \
                 args.outfol, args.max_mismatch, not args.ignore_lengths, args.epsilon, args.max_iterations, \
                 args.min_bins, args.max_bins, args.min_reads, args.dense_region_coverage, args.length_minimum, \
                 args.length_maximum, args.usage, args.use_theta)