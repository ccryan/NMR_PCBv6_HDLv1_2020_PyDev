'''
Created on Mar 30, 2018

@author: David Ariando
'''

#!/usr/bin/python

import os
import time

from nmr_std_function.data_parser import parse_simple_info
from nmr_std_function.nmr_functions import compute_iterate
from nmr_std_function.nmr_class import tunable_nmr_system_2018

# variables
data_folder = "/root/NMR_DATA"
en_fig = 1

# system setup
nmrObj = tunable_nmr_system_2018(data_folder)
nmrObj.turnOnRemoteDebug()
nmrObj.initNmrSystem()
nmrObj.turnOnPower()
nmrObj.setPreampTuning()
nmrObj.setMatchingNetwork()
nmrObj.setSignalPath()

# cpmg settings
cpmg_freq = 4.268  # 4.253 original
pulse1_us = 5  # pulse pi/2 length
pulse2_us = pulse1_us * 1.6  # pulse pi length
pulse1_dtcl = 0.5  # useless with current code
pulse2_dtcl = 0.5  # useless with current code
echo_spacing_us = 150
scan_spacing_us = 400000
samples_per_echo = 1024  # number of points
echoes_per_scan = 64  # number of echos
init_adc_delay_compensation = 10  # acquisition shift microseconds
number_of_iteration = 4  # number of averaging
ph_cycl_en = 1
pulse180_t1_int = 0
delay180_t1_int = 0

nmrObj.cpmgSequence(cpmg_freq, pulse1_us, pulse2_us, pulse1_dtcl, pulse2_dtcl, echo_spacing_us, scan_spacing_us, samples_per_echo,
                    echoes_per_scan, init_adc_delay_compensation, number_of_iteration, ph_cycl_en, pulse180_t1_int, delay180_t1_int)

nmrObj.turnOffSystem()

meas_folder = parse_simple_info(data_folder, 'current_folder.txt')
(a, a0, snr, T2, noise, res, theta, data_filt, echo_avg, Df) = compute_iterate(
    data_folder, meas_folder[0], en_fig)
