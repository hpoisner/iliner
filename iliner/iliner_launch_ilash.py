#!/bin/usr/python3
from datetime import datetime
import logging
import numpy as np
import os
import subprocess
import sys
import time


class RunIlash(object):
    """RunIlash will prepare the parameter file for ilash and execute ilash.
    The output will be the ilash results.
    """

    def make_param_file(mapfile, pedfile, ilash_path):
        """make_param_file() makes the parameter file for running ilash and launches ilash.
        
        Args:
            mapfile(:obj:`pathway`): pathway to map file 
            pedfile(:obj:`pathway`): pathway to pedigree file
        """
        output_str = pedfile.strip('.ped')
        output_str = output_str + '_out.ilash'
        parameter_file = output_str + '_parameters.txt'
        ## Creating the paramerter file
        paramfile = open(parameter_file, "w")
        paramfile.write('map ' + mapfile + '\n')
        paramfile.write('ped ' + pedfile + '\n')
        paramfile.write('perm_count 12' + '\n')
        paramfile.write('shingle_size 20' + '\n')
        paramfile.write('shingle_overlap 0' + '\n')
        paramfile.write('bucket_count 4' + '\n')
        paramfile.write('max_thread 20' + '\n')
        paramfile.write('match_threshold 0.99' + '\n')
        paramfile.write('interest_threshold 0.70' + '\n')
        paramfile.write('max_error 0' + '\n')
        paramfile.write('min_length 2.9' + '\n')
        paramfile.write('auto_slice 1' + '\n')
        paramfile.write('cm_overlap 1.4' + '\n')
        paramfile.write('output ' + output_str + '\n')
        paramfile.close()
        
        subprocess.call([ilash_path, parameter_file])
        
"""
logfilename = str(sys.argv[1]).strip(".sample")
logfilename = logfilename + ".log"
logging.basicConfig(filename=logfilename, level=logging.INFO)
## Load files and uses class from input_code
start_lf = datetime.now()
sample_file = LoadFiles.load_sample(sys.argv[1])
haps_file = LoadFiles.load_haps(sys.argv[2])
stop_lf = datetime.now()
time_diff_lf = stop_lf - start_lf
print("Loading your files takes {}".format(time_diff_lf))

## Parsing section
start_ps = datetime.now()
FileParser.check_same_sample(sample_file, haps_file)
FileParser.check_alleles(haps_file)

nucfunc = np.vectorize(FileParser.nucleotide_test)
log_array_ref = nucfunc(haps_file[:,3].astype(str))
log_array_alt = nucfunc(haps_file[:,4].astype(str))

FileParser.validate_nucleotides(log_array_ref)
FileParser.validate_nucleotides(log_array_alt)
FileParser.validate_binary(haps_file)
stop_ps = datetime.now()
time_diff_ps = stop_ps - start_ps
print("Parsing your files takes {}".format(time_diff_ps))

## File Convert Section
start_fc = datetime.now()
pedfile_start = FileConvert.start_pedfile(sample_file)
pedlen = int(len(pedfile_start))

hapmap = open(sys.argv[3])

outfile = open(sys.argv[4], "w")
FileConvert.make_gpos_mapfile(haps_file, hapmap, outfile)

hapslist = FileConvert.convert_haps(haps_file)
final_list = FileConvert.make_pedfile(hapslist, pedlen, pedfile_start)

np.savetxt(sys.argv[5], final_list, fmt='%s', delimiter=' ')

stop_fc = datetime.now()
time_diff_fc = stop_fc - start_fc
print("Making the map file and plink ped file takes {}".format(time_diff_fc))


RunIlash.make_param_file(sys.argv[4], sys.argv[5], sys.argv[6])
"""
