'''
Created on May 04, 2020

This module characterizes the preamp gain and show the gain over frequency

@author: David Ariando

Edits: Cheng Chen, 07/2020, add automatic tuning
'''

import os
import time
from nmr_std_function.nmr_functions import compute_iterate
from nmr_std_function.nmr_functions import compute_gain
from nmr_std_function.data_parser import parse_simple_info
from nmr_std_function.nmr_class import tunable_nmr_system_2018
from nmr_std_function import data_parser

# variables
data_parent_folder = "/root/NMR_DATA"
en_remote_dbg = 0
fig_num = 1
en_fig = 1

para_folder = "/root/nmr_pcb20_hdl10_2018/MAIN_nmr_code/para"
load_para = 1

# measurement properties
sta_freq = 1
sto_freq = 3
spac_freq = 0.01
samp_freq = 25

target_freq = 1.97

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

if ((target_freq < sta_freq) | (target_freq > sto_freq)):
    print("Warning: Target frequency out of range")

if (load_para):
    # parameter from 
    ( FreqList_S21, PeakVoltage, VvaracList, VbiasList ) = data_parser.parse_csv_float4col_s11( 
        para_folder, '/genS21Table_input_10k.txt' )  # read file
    Vbias = VbiasList[[i for i, elem in enumerate( FreqList_S21 ) if abs( elem - target_freq) < 0.01][0]]
    Vvarac = VvaracList[[i for i, elem in enumerate( FreqList_S21 ) if abs( elem - target_freq) < 0.01][0]]
    
else:
    Vbias = -2.0
    Vvarac = 2.8 

# (vbias vvarac)
nmrObj.setPreampTuning(Vbias, Vvarac)  #-2.5,  2.8)  # try -2.7, -1.8 if fail 
nmrObj.setMatchingNetwork( 0, 0 )
nmrObj.setMatchingNetwork( 0, 0 )

while True:
    nmrObj.assertControlSignal( nmrObj.RX1_1L_msk | nmrObj.RX2_L_msk | nmrObj.RX2_H_msk | nmrObj.RX_SEL1_msk |
                               nmrObj.RX_FL_msk )
    nmrObj.assertControlSignal( nmrObj.PSU_15V_TX_P_EN_msk | nmrObj.PSU_15V_TX_N_EN_msk | nmrObj.PSU_5V_TX_N_EN_msk |
                               nmrObj.PSU_5V_ADC_EN_msk | nmrObj.PSU_5V_ANA_P_EN_msk |
                               nmrObj.PSU_5V_ANA_N_EN_msk )
    time.sleep( 0.1 )

    nmrObj.pamp_char ( sta_freq, sto_freq, spac_freq, samp_freq )

    nmrObj.deassertAll()

    meas_folder = parse_simple_info( data_parent_folder, 'current_folder.txt' )
    # meas_folder[0] = "nmr_wobb_2018_06_25_12_44_48"

    maxS21, maxS21_freq, _ = compute_gain( nmrObj, data_parent_folder, meas_folder[0], en_fig, fig_num )
    print( 'maxS21={0:0.2f} maxS21_freq={1:0.2f}'.format( 
         maxS21, maxS21_freq ) )
    
    dst_pat = ( data_parent_folder + '/' +  meas_folder[0] + '/' )
    data_parser.write_text_append(dst_pat, 'acqu.par', 'target_freq = {}'.format(target_freq))
    data_parser.write_text_append(dst_pat, 'acqu.par', 'Vbias = {}'.format(Vbias))
    data_parser.write_text_append(dst_pat, 'acqu.par', 'Vvarac = {}'.format(Vvarac))
