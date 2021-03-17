# -*- coding: utf-8 -*-
"""
Created on Wed Mar 17 10:43:57 2021
Script to analyze oscilloscope trace.
@author: matias
"""

#%% Importing libraries
import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import find_peaks
from scipy import signal
import os

#%% Function definitions
def calc_snr(f,p,f0):
    # function to calculate signal to noie Ratio
    qq_base = np.where((f>f0-15e6)&(f<f0-5e6))
    qq_peak = np.where((f>f0-10e6)&(f<f0+10e6))
    SNR = np.max(p[qq_peak[0]]) - np.mean(p[qq_base[0]])
    
    return SNR

#%% Get files in folder path
path = 'gator_amp_2/'
files = os.listdir(path)    
files = list(filter(lambda f: f.endswith('.csv'), files))

time = np.array([])
data = []

for f in files:     
    my_data = np.genfromtxt(path+f, delimiter=',',skip_header=6)
    data.append(my_data[:,4])
    if len(time)==0:
        time = my_data[:,3]

data = np.array(data)

#%% Varying the delay pulse.

delay_array = np.linspace(-0.4e-9, 0.4e-9,51)
P_delay = []

for ii in delay_array:
    pwm = signal.square(2*np.pi*250e6*(time+ii),duty=0.25)+1
    sig_gate = data*pwm
    X_gate = np.fft.fft(sig_gate)
    P_med = 10*np.log10(np.mean(np.abs(X_gate)**2,axis=0)/50/1e-3)
    P_delay.append(P_med)

#%% Varying the duty cycle.


#%% Gated approach
    
#%% A few attempts with FFT
n = time.size
timestep = 40e-12
frequ = np.linspace(0, 1/timestep, n)

X_data = np.fft.fft(data)


plt.figure(),plt.semilogy(frequ,np.mean(X_data**2,axis=0))
