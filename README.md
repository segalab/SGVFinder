# ICRA and SGVFinder

This code corrects read assignments by coverage based redistribution
of ambiguously mapped reads. It then uses these correded assignments
to detect structural variants that are either variable across a cohort
or deleted across 25-75% of it. 
This code was used for the paper "Structural variation in the gut 
microbiome associates with host health", TBP. 


## Requirements

1.  The GEM Mapper (Marco-Sola et al. Nat. Methods 2014), available 
    from https://sourceforge.net/projects/gemlibrary/. Please note 
    that this version will not work with GEM3. You need to have the 
    binaries in your PATH variables, so that a simple call to 
    "gem-mapper" would be successful.
2.  This code was written and tested on python 2.7.8, and requires the following packages:
    - numpy (tested with 1.14.2)
    - biopython (tested with 1.68)
    - ujson (tested with 1.35)
    - pandas (tested with 0.23.4)
    - scipy (tested with 1.1.0)
    - bokeh (tested with 0.12.6)

    If you encounter issues, please try to run in an environment with
    these packages.
3. It additionally requires c++ 11 and cython installed.
    
## Install

1. Download the files from https://weizmann.box.com/v/SGVF-DataFiles 
   to the folder containing the code files. 
2. cat the files to a single archive, using: 
```
   cat DataFiles.tar.gz.xaa DataFiles.tar.gz.xab DataFiles.tar.gz.xac DataFiles.tar.gz.xad > DataFiles.tar.gz
```
3. Extract DataFiles.tar.gz in the root folder of the project (where this file is located). 
4. From the ```cy_ext``` subfolder, run ```python setup.py build_ext```.

(DataFiles.tar.gz also available via Zenodo at https://zenodo.org/record/3237975#.XPc0uJxRWhc)

## Usage

There are two main algorithms here - ICRA and SGVFinder.

### ICRA
ICRA has just a single method needed to operate it - ```single_file```. You 
can use it directly from python (recommended), or run it using the 
command-line wrapper ```ICRA_cmd.py```. This method takes in a (/pair of) 
fastq files and outputs a jsdel file. This file is a json file saved
with python's ujson package. It's a dictionary whose keys are the fastq
read ids, and the values are mapping lists. Each such mapping list is
a list of tuples, where the items in the tuple are: the destination id
in the database, the position of the first read, the position of the 
second read (-1 if SE), the probablity ICRA gives to this mapping, 
and the mapping quality.
You should run that method on each and every sample in your cohort.

### SGVFinder
SGVFinder has two stages, and hence two methods:

```get_sample_map``` - generates coverage maps ber bacteria per sample. You 
can use it directly from python, or run it using the command-line 
wrapper ```SGVF_PerFile_cmd.py```. You should run this method on the jsdel file
of each and every sample in your cohort.

```work_on_collection``` - generates the SV dataframes. You can use it
directly from python or run it using the command-line wrapper ```SGVF_cmd.py```.
You should only run this method once. It takes as input a dictionary
whose keys are the sample names and whose values are the sample_maps 
generated using ```get_sample_map```. This is generated automatically from a
glob string with the command-line wrapper.

**NOTE:** SGVFinder WILL NOT work on a single sample. If you have a small 
cohort we recommend changing the ```min_samp_cutoff``` or running with ```--byorig```.


**See the linear_example.py for a non-parallelized simple implementation.**

