'''
Created on Mar 30, 2018

@author: David Ariando

in progress:
*** direct_read, a setting to read data directly from SDRAM using Python is implemented using the same CPMG_iterate C-function (with the loop handled by Python instead of C).
It is a slower operation than using C to write to SDCard directly. It might be due to every iteration will
generate acqu.par separately in separate folder. Or MAYBE [less likely] writing to SDRAM via DMA is slower than writing directly
to a text file. Or MAYBE because every iteration runs CPMG_iterate with 1-iteration, the C program is called
number_of_iteration times. Also the computation of SNR/echo/scan is wrong in compute_iterate due to having
only 1 as a number_of_iteration in acqu.par. It seems to be no point of using Python if C still needs to be called.
So many problems yet not having better speed. It might be better to write everything in Python instead.

Cheng 07/2020 Option to change matching network and preamp values automatically
'''

#!/usr/bin/python

# settings
data_folder = "/root/NMR_DATA"  # the nmr data folder
en_fig = 1  # enable figure
en_remote_dbg = 0  # enable remote debugging. Enable debug server first!
direct_read = 0  # perform direct read from SDRAM. use with caution above!
meas_time = 1  # measure time
process_data = 0  # process data within the SoC

para_folder = "/root/nmr_pcb20_hdl10_2018/MAIN_nmr_code/para"
load_para = 1

import time
if ( meas_time ):
    start_time = time.time()

import os
from nmr_std_function.data_parser import parse_simple_info
from nmr_std_function.nmr_functions import compute_iterate
from nmr_std_function.nmr_functions import calcP90
from nmr_std_function.nmr_class import tunable_nmr_system_2018
from nmr_std_function.data_parser import parse_csv_float2col

from nmr_std_function import data_parser
import matplotlib.pyplot as plt
from scipy import signal
import pydevd

if ( meas_time ):
    elapsed_time = time.time() - start_time
    print( "load library time: %.3f" % ( elapsed_time ) )
    start_time = time.time()

# cpmg settings
cpmg_freq = 2.0 #4.286 + ( 9 - 99 + 20 + 19 - 32.6 - 3.45 + 22.8 - 11 + 7.65 ) * 1e-3
pulse1_us = 30  # 75 for Cheng's coil. pulse pi/2 length.
pulse2_us = 37 # pulse pi length
pulse1_dtcl = 0.5  # useless with current code
pulse2_dtcl = 0.5  # useless with current code
echo_spacing_us = 300  # 200
scan_spacing_us = 20000
samples_per_echo = 1024  # 3072
echoes_per_scan = 30  # 20
# put to 10 for broadband board and 6 for tunable board
init_adc_delay_compensation = 20 #6  # acquisition shift microseconds.
number_of_iteration = 100  # number of averaging
ph_cycl_en = 1
pulse180_t1_int = 0
delay180_t1_int = 0
tx_sd_msk = 1  # 1 to shutdown 8tx opamp during reception, or 0 to keep it powered up during reception
en_dconv = 1  # enable downconversion in the fpga
dconv_fact = 4  # downconversion factor. minimum of 4.

# coil param and measured voltage across the coil
Vpp = 312  # 190
rs = 1.2
L = 2.438e-6
coilLength = 36e-3
numTurns = 37
coilFactor = 0.675  # measured_eff_p90/calc'ed_p90. Equal to 1 for calc'ed_p90
# magnet param
B0 = 0.044  # T
gamma = 42.57  # MHz/T
print( "freq estimate: %3.3f MHz" % ( gamma * B0 ) )
P90, Pwatt = calcP90( Vpp, rs, L, cpmg_freq * 1e6,
                     numTurns, coilLength, coilFactor )
print( "P90 len estimate: %3.3f us, power estimate: %3.3f Watts" %
      ( P90 * 1e6, Pwatt ) )

# instantiate nmr object
nmrObj = tunable_nmr_system_2018( data_folder, en_remote_dbg )

# system setup
nmrObj.initNmrSystem()  # necessary to set the GPIO initial setting. Also fix the
nmrObj.assertControlSignal( nmrObj.PSU_15V_TX_P_EN_msk | nmrObj.PSU_15V_TX_N_EN_msk | nmrObj.PSU_5V_TX_N_EN_msk |
                           nmrObj.PSU_5V_ADC_EN_msk | nmrObj.PSU_5V_ANA_P_EN_msk |
                           nmrObj.PSU_5V_ANA_N_EN_msk )
# nmrObj.deassertControlSignal(
#    nmrObj.PSU_15V_TX_P_EN_msk | nmrObj.PSU_15V_TX_N_EN_msk)

freq_comp = cpmg_freq #+ 0.03
freqS21_comp = cpmg_freq

if (load_para):
    # parameter from 
    ( FreqList, s11List, CparList, CserList ) = data_parser.parse_csv_float4col_s11( 
        para_folder, '/genS11Table_final_input_10k.txt' )  # read file
    Cpar = int(CparList[[i for i, elem in enumerate( FreqList ) if abs( elem - freq_comp) < 0.01][0]])
    Cser = int(CserList[[i for i, elem in enumerate( FreqList ) if abs( elem - freq_comp) < 0.01][0]])
    
    ( FreqList_S21, PeakVoltage, VvaracList, VbiasList ) = data_parser.parse_csv_float4col_s11( 
        para_folder, '/genS21Table_input_10k.txt' )  # read file
    Vbias = VbiasList[[i for i, elem in enumerate( FreqList_S21 ) if abs( elem - freq_comp) < 0.01][0]]
    Vvarac = VvaracList[[i for i, elem in enumerate( FreqList_S21 ) if abs( elem - freq_comp) < 0.01][0]]
    
else:
    Cpar = 563
    Cser = 327
    Vbias = -2.0
    Vvarac = 2.8    

#print(Cpar, Cser, Vbias, Vvarac)
nmrObj.setPreampTuning(Vbias, Vvarac)  # try -2.7, -1.8 if fail
nmrObj.setMatchingNetwork(Cpar,    Cser)  # 4.25 MHz AFE
nmrObj.setMatchingNetwork(Cpar,    Cser)

if ( nmrObj.PCBVer == 'v4.0_and_below' ):
    nmrObj.assertControlSignal( nmrObj.AMP_HP_LT1210_EN_msk |
                               nmrObj.PAMP_IN_SEL_RX_msk | nmrObj.RX_IN_SEL_1_msk )
elif ( nmrObj.PCBVer == 'v5.0' ):
    nmrObj.assertControlSignal( 
        nmrObj.RX1_1L_msk | nmrObj.RX1_1H_msk | nmrObj.RX2_L_msk | nmrObj.RX2_H_msk | nmrObj.RX_SEL1_msk | nmrObj.RX_FL_msk | nmrObj.RX_FH_msk | nmrObj.PAMP_IN_SEL2_msk )
nmrObj.deassertControlSignal( nmrObj.RX1_1H_msk | nmrObj.RX_FH_msk )

if ( meas_time ):
    elapsed_time = time.time() - start_time
    print( "set parameter time: %.3f" % ( elapsed_time ) )
    start_time = time.time()

if ( direct_read ):
    datain = nmrObj.cpmgSequenceDirectRead( cpmg_freq, pulse1_us, pulse2_us, pulse1_dtcl, pulse2_dtcl, echo_spacing_us, scan_spacing_us, samples_per_echo,
                                           echoes_per_scan, init_adc_delay_compensation, number_of_iteration, ph_cycl_en,
                                           pulse180_t1_int, delay180_t1_int , tx_sd_msk )
else:
    nmrObj.cpmgSequence( cpmg_freq, pulse1_us, pulse2_us, pulse1_dtcl, pulse2_dtcl, echo_spacing_us, scan_spacing_us, samples_per_echo,
                        echoes_per_scan, init_adc_delay_compensation, number_of_iteration,
                        ph_cycl_en, pulse180_t1_int, delay180_t1_int , tx_sd_msk, en_dconv, dconv_fact )
    datain = []  # set datain to 0 because the data will be read from file instead

if ( meas_time ):
    elapsed_time = time.time() - start_time
    print( "cpmgSequence acquisition time: %.3f" % ( elapsed_time ) )
    start_time = time.time()

if ( nmrObj.PCBVer == 'v4.0_and_below' ):
    nmrObj.deassertControlSignal( nmrObj.AMP_HP_LT1210_EN_msk |
                                 nmrObj.PAMP_IN_SEL_RX_msk | nmrObj.RX_IN_SEL_2_msk )
elif ( nmrObj.PCBVer == 'v5.0' ):
    nmrObj.deassertControlSignal( 
        nmrObj.RX1_1H_msk | nmrObj.RX1_1L_msk | nmrObj.RX2_L_msk | nmrObj.RX2_H_msk | nmrObj.RX_SEL1_msk | nmrObj.RX_FL_msk | nmrObj.RX_FH_msk | nmrObj.PAMP_IN_SEL2_msk )

nmrObj.setMatchingNetwork( 0, 0 )
nmrObj.setPreampTuning( 0, 0 )
nmrObj.deassertControlSignal( nmrObj.PSU_15V_TX_P_EN_msk | nmrObj.PSU_15V_TX_N_EN_msk | nmrObj.PSU_5V_TX_N_EN_msk |
                             nmrObj.PSU_5V_ADC_EN_msk | nmrObj.PSU_5V_ANA_P_EN_msk | nmrObj.PSU_5V_ANA_N_EN_msk )

if ( process_data ):
    meas_folder = parse_simple_info( data_folder, 'current_folder.txt' )
    ( a, a_integ, a0, snr, T2, noise, res, theta, data_filt, echo_avg, Df, t_echospace ) = compute_iterate( 
      nmrObj, data_folder, meas_folder[0], 0, 0, 0, direct_read, datain, en_fig )

if ( meas_time ):
    elapsed_time = time.time() - start_time
    print( "data processing time: %.3f" % ( elapsed_time ) )
    start_time = time.time()
