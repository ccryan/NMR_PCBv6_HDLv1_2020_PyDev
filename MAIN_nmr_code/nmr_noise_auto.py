'''
Created on Mar 30, 2018

@author: David Ariando

Cheng 07/2020 Option to change matching network and preamp values automatically
'''

#!/usr/bin/python

import os
import time

from nmr_std_function.data_parser import parse_simple_info
from nmr_std_function.nmr_functions import compute_iterate
from nmr_std_function.nmr_functions import compute_stats
from nmr_std_function.nmr_class import tunable_nmr_system_2018
from nmr_std_function.data_parser import parse_csv_float2col
from nmr_std_function import data_parser
import matplotlib.pyplot as plt
from scipy import signal
import pydevd
from datetime import datetime
import shutil

# variables
data_folder = "/root/NMR_DATA"
en_fig = 1
en_remote_dbg = 0

para_folder = "/root/nmr_pcb20_hdl10_2018/MAIN_nmr_code/para"
load_para = 1
target_freq =  2.03

# define the name of the directory to be created
same_folder = 1
now = datetime.now()
datename = now.strftime( "%Y_%m_%d_%H_%M_%S" )
dst_path = data_folder + '/' + datename + '_noiseFolder'

# nmr object declaration
nmrObj = tunable_nmr_system_2018( data_folder, en_remote_dbg )

# measurement settings
samp_freq = 25  # sampling frequency
samples = 40000  # number of points
min_freq = 0.001
max_freq = 4#12.5

# system setup
nmrObj.initNmrSystem()  # necessary to set the GPIO initial setting
nmrObj.assertControlSignal( nmrObj.PSU_15V_TX_P_EN_msk | nmrObj.PSU_15V_TX_N_EN_msk | nmrObj.PSU_5V_TX_N_EN_msk |
                           nmrObj.PSU_5V_ADC_EN_msk | nmrObj.PSU_5V_ANA_P_EN_msk |
                           nmrObj.PSU_5V_ANA_N_EN_msk )
nmrObj.deassertControlSignal( 
    nmrObj.PSU_15V_TX_P_EN_msk | nmrObj.PSU_15V_TX_N_EN_msk )

freq_comp = target_freq+0.10
freqS21_comp = target_freq
if (load_para):
    # parameter from 
    ( FreqList, s11List, CparList, CserList ) = data_parser.parse_csv_float4col_s11( 
        para_folder, '/genS11Table_final_input_10k.txt' )  # read file
    Cpar = int(CparList[[i for i, elem in enumerate( FreqList ) if abs( elem - freq_comp) < 0.01][0]])
    Cser = int(CserList[[i for i, elem in enumerate( FreqList ) if abs( elem - freq_comp) < 0.01][0]])
    
    #Cpar = 0
    #Cser = 0
    
    ( FreqList_S21, PeakVoltage, VvaracList, VbiasList ) = data_parser.parse_csv_float4col_s11( 
        para_folder, '/genS21Table_input_10k.txt' )  # read file
    Vbias = VbiasList[[i for i, elem in enumerate( FreqList_S21 ) if abs( elem - freqS21_comp) < 0.01][0]]
    Vvarac = VvaracList[[i for i, elem in enumerate( FreqList_S21 ) if abs( elem - freqS21_comp) < 0.01][0]]
    
else:
    Cpar = 563
    Cser = 327
    Vbias = -2.0
    Vvarac = 2.8

print(Cpar, Cser)
print(Vbias, Vvarac)

nmrObj.setPreampTuning(Vbias, Vvarac) #-2.5,  2.8)#-2.5,  2.6)# -2.7, 0.3 )  # try -2.7, -1.8 if fail
nmrObj.setMatchingNetwork(Cpar,    Cser)  # 4.25 MHz AFE
nmrObj.assertControlSignal( 
    nmrObj.RX_FL_msk | nmrObj.RX_FH_msk | nmrObj.RX_SEL1_msk | nmrObj.RX2_L_msk | nmrObj.RX2_H_msk | nmrObj.RX1_1L_msk | nmrObj.RX1_1H_msk | nmrObj.PAMP_IN_SEL2_msk )
nmrObj.deassertControlSignal( nmrObj.RX_FH_msk | nmrObj.RX2_L_msk | nmrObj.RX_FH_msk )

i = 0
now = datetime.now()
datename = now.strftime( "%Y_%m_%d_%H_%M_%S" )
dst_path = data_folder + '/' + datename + '_noise'
nmrObj.noise( samp_freq, samples )
while True:

    # time.sleep(0.5)

    nmrObj.noise( samp_freq, samples )

    # process the data
    meas_folder = parse_simple_info( data_folder, 'current_folder.txt' )
    
    if same_folder:
        src_file = ( data_folder + '/' + meas_folder[0] + '/asum')
        dst_file = ( dst_path + '/asum_{}'.format(i))
        shutil.copy2(src_file, dst_file)
    
        i = i + 1
    
    compute_stats( min_freq, max_freq, data_folder, meas_folder[0], 'noise_plot.png', en_fig )

nmrObj.setMatchingNetwork( 0, 0 )
nmrObj.setPreampTuning( 0, 0 )
