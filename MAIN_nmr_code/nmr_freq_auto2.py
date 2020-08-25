'''
Created on Mar 30, 2018

@author: David Ariando

Edits: Cheng Chen, 07/2020, add automatic tuning
08/2020, run in 10KHz spacing
Just one frequency
'''

#!/usr/bin/python

import os
import time

from nmr_std_function.data_parser import parse_simple_info
from nmr_std_function.nmr_functions import compute_iterate
from nmr_std_function.nmr_class import tunable_nmr_system_2018
from nmr_std_function.data_parser import parse_csv_float2col
import matplotlib.pyplot as plt
from scipy import signal
import pydevd
import numpy as np
from datetime import datetime
import shutil
from nmr_std_function import data_parser

# variables
data_folder = "/root/NMR_DATA"
en_scan_fig = 0
en_fig = 0
en_remote_dbg = 0
fig_num = 100
direct_read = 0  # perform direct read from SDRAM. use with caution above!

tx_sd_msk = 1  # 1 to shutdown 8tx opamp during reception, or 0 to keep it powered up during reception
en_dconv = 0  # enable downconversion in the fpga
dconv_fact = 4  # downconversion factor. minimum of 4.

process_data = 0
para_folder = "/root/nmr_pcb20_hdl10_2018/MAIN_nmr_code/para"
load_para = 1 # load parameter for matching network and preamp

# instantiate nmr object
nmrObj = tunable_nmr_system_2018( data_folder, en_remote_dbg )

# cpmg settings
#cpmg_freq = 1.99
pulse1_dtcl = 0.5  # useless with current code
pulse2_dtcl = 0.5  # useless with current code
echo_spacing_us = 350
scan_spacing_us = 22000
samples_per_echo = 1024  # number of points
echoes_per_scan = 40  # number of echos
init_adc_delay_compensation = 20  # acquisition shift microseconds
number_of_iteration = 10 # number of averaging
ph_cycl_en = 1
pulse180_t1_int = 0
delay180_t1_int = 0

pulse1_us = 55  # pulse pi/2 length
pulse2_us = 75 #pulse_us_sw[i]  # pulse pi length

# sweep settings
Freq_step = 1 # number of steps

cpmg_freq_sta = 2.18 # in MHz
cpmg_freq_sto = 2.18  # in MHz
cpmg_freq_sw = np.linspace(cpmg_freq_sta, cpmg_freq_sto, Freq_step)

cpmg_freq_sw[0]=2.04
cpmg_freq_sta = number_of_iteration # in MHz
cpmg_freq_sto = number_of_iteration  # in MHz
cpmg_NI_sw = np.linspace(cpmg_freq_sta, cpmg_freq_sto, Freq_step)

cpmg_NI_sw[0]=8000

# system setup
nmrObj.initNmrSystem()  # necessary to set the GPIO initial setting
# nmrObj.turnOnPower()
nmrObj.assertControlSignal( nmrObj.PSU_15V_TX_P_EN_msk | nmrObj.PSU_15V_TX_N_EN_msk | nmrObj.PSU_5V_TX_N_EN_msk |
                           nmrObj.PSU_5V_ADC_EN_msk | nmrObj.PSU_5V_ANA_P_EN_msk |
                           nmrObj.PSU_5V_ANA_N_EN_msk )



now = datetime.now()
datename = now.strftime( "%Y_%m_%d_%H_%M_%S" )

# define the name of the directory to be created
dst_path = data_folder + '/' + datename + '_FreqSweep'

try:
    os.mkdir(dst_path)
except OSError:
    print ("Creation of the directory %s failed" % dst_path)
else:
    print ("Successfully created the directory %s " % dst_path)

a_integ_table = np.zeros( Freq_step )
for i in range( 0, Freq_step ):
    print( '----------------------------------' )
    print( 'frequency = ' + str( cpmg_freq_sw[i] ) + ' MHz' )

    
    #cpmg_freq=cpmg_freq_sw[i]
    
     # compensate for setup 
    freq_comp = cpmg_freq_sw[i] #+ 0.03
    freqS21_comp = cpmg_freq_sw[i]
    number_of_iteration=cpmg_NI_sw[i]
    if (load_para):
        # parameter from 
        ( FreqList, s11List, CparList, CserList ) = data_parser.parse_csv_float4col_s11( 
            para_folder, '/genS11Table_final_input_10k.txt' )  # read file
        Cpar = int(CparList[[i for i, elem in enumerate( FreqList ) if abs( elem - freq_comp) < 0.01][0]])
        Cser = int(CserList[[i for i, elem in enumerate( FreqList ) if abs( elem - freq_comp) < 0.01][0]])
        
        ( FreqList_S21, PeakVoltage, VvaracList, VbiasList ) = data_parser.parse_csv_float4col_s11( 
            para_folder, '/genS21Table_input_10k.txt' )  # read file
        Vbias = VbiasList[[i for i, elem in enumerate( FreqList_S21 ) if abs( elem - freqS21_comp) < 0.01][0]]
        Vvarac = VvaracList[[i for i, elem in enumerate( FreqList_S21 ) if abs( elem - freqS21_comp) < 0.01][0]]
    else:
        Cpar = 563
        Cser = 327
        Vbias = -2.0
        Vvarac = 2.8
    
    nmrObj.setPreampTuning(Vbias, Vvarac)
    nmrObj.setMatchingNetwork(Cpar, Cser)

    nmrObj.assertControlSignal( 
    nmrObj.RX1_1H_msk | nmrObj.RX1_1L_msk | nmrObj.RX2_L_msk | nmrObj.RX2_H_msk | nmrObj.RX_SEL1_msk | nmrObj.RX_FL_msk | nmrObj.RX_FH_msk | nmrObj.PAMP_IN_SEL2_msk )
     
    nmrObj.cpmgSequence( cpmg_freq_sw[i], pulse1_us, pulse2_us, pulse1_dtcl, pulse2_dtcl, echo_spacing_us, scan_spacing_us, samples_per_echo,
                        echoes_per_scan, init_adc_delay_compensation, number_of_iteration, ph_cycl_en, pulse180_t1_int, delay180_t1_int,
                         tx_sd_msk, en_dconv, dconv_fact  )
    datain = []  # set datain to 0 because the data will be read from file instead
    meas_folder = parse_simple_info( data_folder, 'current_folder.txt' )

    if  en_dconv:
        src_file = ( data_folder + '/' + meas_folder[0] + '/dconv')
        dst_file = ( dst_path + '/dconv_{}'.format(i))
        shutil.copy2(src_file, dst_file)
    else:
        src_file = ( data_folder + '/' + meas_folder[0] + '/asum')
        dst_file = ( dst_path + '/asum_{}'.format(i))
        shutil.copy2(src_file, dst_file)      
    
    if process_data:
        ( a, a_integ, a0, snr, T2, noise, res, theta, data_filt, echo_avg, Df, t_echospace ) = compute_iterate( 
        nmrObj, data_folder, meas_folder[0], 0, 0, 0, direct_read, datain, en_scan_fig)

        a_integ_table[i] = a_integ
        
    if en_fig:
        plt.ion()
        fig = plt.figure( fig_num )
        fig.clf()
        ax = fig.add_subplot( 1, 1, 1 )
        line1, = ax.plot( cpmg_freq_sw[0:i + 1], a_integ_table[0:i + 1], 'b-o' )
        # ax.set_ylim(-50, 0)
        # ax.set_xlabel('Frequency [MHz]')
        # ax.set_ylabel('S11 [dB]')
        # ax.set_title("Reflection Measurement (S11) Parameter")
        ax.grid()
        fig.canvas.draw()
        # fig.canvas.flush_events()

# turn off system
nmrObj.deassertControlSignal( 
    nmrObj.RX1_1H_msk | nmrObj.RX1_1L_msk | nmrObj.RX2_L_msk | nmrObj.RX2_H_msk | nmrObj.RX_SEL1_msk | nmrObj.RX_FL_msk | nmrObj.RX_FH_msk | nmrObj.PAMP_IN_SEL2_msk )

nmrObj.setMatchingNetwork( 0, 0 )
nmrObj.setPreampTuning( 0, 0 )
nmrObj.deassertControlSignal( nmrObj.PSU_15V_TX_P_EN_msk | nmrObj.PSU_15V_TX_N_EN_msk | nmrObj.PSU_5V_TX_N_EN_msk |
                             nmrObj.PSU_5V_ADC_EN_msk | nmrObj.PSU_5V_ANA_P_EN_msk | nmrObj.PSU_5V_ANA_N_EN_msk )

shutil.copy2(data_folder + '/' + meas_folder[0] + '/acqu.par', dst_path)

for kk in range( 0, Freq_step ):
    data_parser.write_text_append(dst_path, 'acqu.par', 'Freq = {}'.format(cpmg_freq_sw[kk]))

if en_fig:
    fig.savefig( dst_path + '/' + datename + '_pulsesw.pdf' )


pass
pass
