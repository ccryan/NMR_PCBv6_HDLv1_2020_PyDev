'''
Created on Oct 30, 2018

@author: David Ariando
'''

import os
import time

from nmr_std_function.nmr_functions import compute_iterate
from nmr_std_function.nmr_functions import compute_wobble
from nmr_std_function.data_parser import parse_simple_info
from nmr_std_function import data_parser
from nmr_std_function.nmr_class import tunable_nmr_system_2018

# variables
data_parent_folder = "/root/NMR_DATA"

en_remote_dbg = 0
fig_num = 1
en_fig = 1
meas_time = 1

para_folder = "/root/nmr_pcb20_hdl10_2018/MAIN_nmr_code/para"
load_para = 1

# measurement properties
sta_freq = 1.5
sto_freq = 2.5
spac_freq = 0.01
samp_freq = 25
target_freq = 2.2

# instantiate nmr object
nmrObj = tunable_nmr_system_2018( data_parent_folder, en_remote_dbg )

work_dir = os.getcwd()
os.chdir( data_parent_folder )

# system setup
nmrObj.initNmrSystem()  # necessary to set the GPIO initial setting

nmrObj.deassertAll()

nmrObj.assertControlSignal( nmrObj.PSU_5V_TX_N_EN_msk |
                           nmrObj.PSU_5V_ADC_EN_msk | nmrObj.PSU_5V_ANA_P_EN_msk |
                           nmrObj.PSU_5V_ANA_N_EN_msk )

nmrObj.setPreampTuning(-1.80, 3.9)#-2.5,  2.6)  # try -2.7, -1.8 if fail


def runExpt( cparVal, cserVal, S11mV_ref, useRef ):
    # useRef: use the pregenerated S11mV_ref as a reference to compute reflection. If this option is 0, then the compute_wobble will instead generated S11 in mV format instead of dB format

    if meas_time:
        start_time = time.time()

    # enable power and signal path
    nmrObj.assertControlSignal( nmrObj.RX1_2L_msk | nmrObj.RX_SEL2_msk | nmrObj.RX_FL_msk )
    nmrObj.assertControlSignal( nmrObj.PSU_15V_TX_P_EN_msk | nmrObj.PSU_15V_TX_N_EN_msk | nmrObj.PSU_5V_TX_N_EN_msk |
                               nmrObj.PSU_5V_ADC_EN_msk | nmrObj.PSU_5V_ANA_P_EN_msk | nmrObj.PSU_5V_ANA_N_EN_msk )

    # change matching network values (twice because sometimes it doesnt' work once due to transient
    nmrObj.setMatchingNetwork( cparVal, cserVal )
    nmrObj.setMatchingNetwork( cparVal, cserVal )

    # do measurement
    nmrObj.wobble( sta_freq, sto_freq, spac_freq, samp_freq )

    # disable all to save power
    nmrObj.deassertAll()

    if meas_time :
        elapsed_time = time.time() - start_time
        start_time = time.time()  # reset the start time
        print( "### time elapsed for running wobble exec: %.3f" % ( elapsed_time ) )

    # compute the generated data
    meas_folder = parse_simple_info( data_parent_folder, 'current_folder.txt' )
    S11dB, S11_fmin, S11_fmax, S11_bw, minS11, minS11_freq = compute_wobble( nmrObj, data_parent_folder, meas_folder[0], -10, S11mV_ref, useRef, en_fig, fig_num )
    print( '\t\tfmin={:0.3f} fmax={:0.3f} bw={:0.3f} minS11={:0.2f} minS11_freq={:0.3f} cparVal={:d} cserVal={:d}'.format( S11_fmin, S11_fmax, S11_bw, minS11, minS11_freq, cparVal, cserVal ) )

    if meas_time :
        elapsed_time = time.time() - start_time
        print( "### time elapsed for compute_wobble: %.3f" % ( elapsed_time ) )

    return S11dB, minS11_freq

if ((target_freq < sta_freq) | (target_freq > sto_freq)):
    print("Warning: Target frequency out of range")

# find reference
print( 'Generate reference.' )
S11mV_ref, minS11Freq_ref = runExpt( 0, 0, 0, 0 )  # background is computed with no capacitor connected -> max reflection

if (load_para):
    # parameter from 
    ( FreqList, s11List, CparList, CserList ) = data_parser.parse_csv_float4col_s11( 
        para_folder, '/genS11Table_final_input.txt' )  # read file
    Cpar = int(CparList[[i for i, elem in enumerate( FreqList ) if abs( elem - target_freq) < 0.05][0]])
    Cser = int(CserList[[i for i, elem in enumerate( FreqList ) if abs( elem - target_freq) < 0.05][0]])
    
else:
    Cpar = 598
    Cser = 352
    
while True:
    runExpt(Cpar,    Cser,  S11mV_ref, 1 )

