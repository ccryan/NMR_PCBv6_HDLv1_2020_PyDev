'''
Created on Mar 30, 2018

@author: Cheng Chen

Converts binary and save data file in matlab format

'''

#!/usr/bin/python

# settings

en_fig = 1  # enable figure
en_remote_dbg = 1  # enable remote debugging. Enable debug server first!
direct_read = 0  # perform direct read from SDRAM. use with caution above!
meas_time = 1  # measure time
process_data = 0  # process data within the SoC

mtch_fltr_sta_idx = 0  # 0 is default or something referenced to SpE, e.g. SpE/4; the start index for match filtering is to neglect ringdown part from calculation
perform_rotation = 1  # perform rotation to the data -> required for reference data for t1 measurement
proc_indv_data = 0  # process individual raw data, otherwise it'll load a sum file generated by C
binary_OR_ascii = 1  # put 1 if the data file uses binary representation, otherwise it is in ascii format
ignore_echoes = 0  # ignore initial echoes #

use_latest_folder = 0

import time
if ( meas_time ):
    start_time = time.time()
    
import math
import csv
import numpy as np
import os
import os
from nmr_std_function.data_parser import parse_simple_info
from nmr_std_function.data_parser import parse_csv_float2col
import matplotlib.pyplot as plt
from scipy import signal
import pydevd
from nmr_std_function import data_parser

data_parent_folder = "V://NMR_DATA"
if ( use_latest_folder ):
    meas_folder = parse_simple_info( data_folder, 'current_folder.txt' )
    data_folder = ( data_parent_folder + '/' + meas_folder + '/' )
else: 
    data_folder = "V:/NMR_DATA/2020_07_31_17_01_34_cpmg/"
    
if ( meas_time ):
    elapsed_time = time.time() - start_time
    print( "load library time: %.3f" % ( elapsed_time ) )
    start_time = time.time()

    # variables from NMR settings
    ( param_list, value_list ) = data_parser.parse_info( 
        data_folder, 'acqu.par' )  # read file
    SpE = int( data_parser.find_value( 
        'nrPnts', param_list, value_list ) )
    NoE = int( data_parser.find_value( 
        'nrEchoes', param_list, value_list ) )
    en_ph_cycle_proc = data_parser.find_value( 
        'usePhaseCycle', param_list, value_list )
    # tE = data_parser.find_value('echoTimeRun', param_list, value_list)
    # Sf = data_parser.find_value(
    #    'adcFreq', param_list, value_list) * 1e6
    # Df = data_parser.find_value(
    #    'b1Freq', param_list, value_list) * 1e6
    # total_scan = int(data_parser.find_value(
    #    'nrIterations', param_list, value_list))
    fpga_dconv = data_parser.find_value( 
        'fpgaDconv', param_list, value_list )
    dconv_fact = data_parser.find_value( 
        'dconvFact', param_list, value_list )

    # compensate for dconv_fact if fpga dconv is used
    if fpga_dconv:
        SpE = int( SpE / dconv_fact )
        Sf = Sf / dconv_fact

    # ignore echoes
    if ignore_echoes:
        NoE = NoE - ignore_echoes
        

if ( proc_indv_data ):
                # read all datas and average it
    data = np.zeros( NoE * SpE )
    for m in range( 1, total_scan + 1 ):
        file_path = ( data_folder + file_name_prefix +
                     '{0:03d}'.format( m ) )
        # read the data from the file and store it in numpy array
        # format
        one_scan = np.array( data_parser.read_data( file_path ) )
        one_scan = ( one_scan - np.mean( one_scan ) ) / \
            total_scan  # remove DC component
        if ( en_ph_cycle_proc ):
            if ( m % 2 ):  # phase cycling every other scan
                data = data - one_scan
            else:
                data = data + one_scan
        else:
            data = data + one_scan
else:
    # read sum data only
    file_path = ( data_folder + 'asum' )
    data = np.zeros( NoE * SpE )

    if binary_OR_ascii:
        data = data_parser.read_hex_float( file_path )  # use binary representation
    else:
        data = np.array( data_parser.read_data( file_path ) )  # use ascii representation
    
with open( data_folder + '/asum.txt', 'a' ) as Table:
    for a in data:
        Table.write( '{:f}\n' .format( a ) )

print('Format Transfer Complete')