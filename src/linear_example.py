#===============================================================================
# This file provides an example of how to use ICRA and SGVFinder to take a folder full of fastq files
# and to generate dataframes containing deletion- and variable-SGVs. 
# NOTE: Memory requirements are at least 20GB per file, and even then - it takes time, depending on the
# sample. So it is highly not recommended to use this file as is, but rather to edit it to work with your
# HPC environment. 
# NOTE2: See the README for requirements etc.
#===============================================================================
from glob import glob
from os.path import join, splitext, basename
from ICRA import single_file
from SGVFinder import get_sample_map, work_on_collection
import os
INPUT_FOLDER = ''
OUTPUT_FOLDER = ''
os.mkdir(join(OUTPUT_FOLDER, 'browser'))
fastqs_1 = glob(join(INPUT_FOLDER, '*_1.fastq'))
fastqs = [f for f in glob(join(INPUT_FOLDER, '*.fastq')) if f not in fastqs_1 and f.replace('_1.fastq', '_2.fastq') not in fastqs_1]

for f in fastqs_1: #Parallelize on your cluster
    single_file(f, f.replace('_1.fastq', '_2.fastq'), OUTPUT_FOLDER, 8, True, 1e-6, 100, 10, 100, 100, 60, 1e5, 2e7, 'genomes', False)
for f in fastqs:#Parallelize on your cluster
    single_file(f, None, OUTPUT_FOLDER, 8, True, 1e-6, 100, 10, 100, 100, 60, 1e5, 2e7, 'genomes', False)
#Parallelize on your cluster
samp_to_map = {splitext(basename(f))[0]: get_sample_map(f, .01, 100, 10) \
               for f in glob(join(OUTPUT_FOLDER, '*.jsdel'))}

vsgv, dsgv = work_on_collection(samp_to_map, 10, 0.25, 0.95, 0.25, 0.125, 0.02, 0.95, 'betaprime', 0.01, 10, 85, join(OUTPUT_FOLDER, 'browser'))
vsgv.to_pickle(join(OUTPUT_FOLDER, 'vsgv.df'))
dsgv.to_pickle(join(OUTPUT_FOLDER, 'dsgv.df'))
