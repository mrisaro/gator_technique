# -*- coding: utf-8 -*-
"""
Created on Fri Mar 19 11:36:17 2021
Script to analyze data GATOR
@author: matias
"""
#%% Importing libraries

from pylab import *
import os
from scipy import signal
from scipy.optimize import curve_fit
from numpy import *
from scipy.signal import *
from scipy.fft import fft
from pandas import *
# from mpl_toolkits.axes_grid.inset_locator import (inset_axes, InsetPosition,mark_inset)
# from matplotlib.ticker import FuncFormatter

#%% Function defining
def normfft(s):
    si=(2/len(s) *abs(fft(s)))**2/2 # multiply by 2 to account for negative freqs, divide by two to account for effective value
    out=si[0:len(signal)//2]
    frout= fftfreq(len(s), ts)[0:len(s)//2]
    return frout,out

def calc_snr(f,p,f0):
    # function to calculate signal to noie Ratio
    qq_base = np.where((f>10e6)&(f<50e6))
    qq_peak = np.where((f>f0-10e6)&(f<f0+10e6))
    SNR = np.max(p[qq_peak[0]]) - np.mean(p[qq_base[0]])
    
    return SNR
#%% Loading one data frame
fr= loadtxt('gator_data_raw/210316_140348_Ch1.csv',usecols=(3,4),delimiter=',',unpack=True)
frep=250e6
ts=40e-12
fs=1/ts

time=fr[0]
signal=fr[1]
gate=square(2*pi*frep*(time-0.02e-9),duty=0.03)+1
prod=gate*signal
fxs,fftsig=normfft(signal)
fxp,fftprod=normfft(prod)

plt.figure('PSDsP')
plt.loglog(fxp,fftsig/max(fftsig[((fxp>90e6)&(fxp<100e6))]),color='gray')
plt.loglog(fxp,fftprod/max(fftprod[((fxp>90e6)&(fxp<100e6))]),color='red')
plt.grid(True,'both')

#%% Doing the loop for all the other measurement.
fft_sig = []
fft_gat = []

P_sig = []
P_gat = []

path = 'gator_data_raw/'
files = os.listdir(path)    
files = list(filter(lambda f: f.endswith('.csv'), files))

for f in files:
    fr= loadtxt(path+f,usecols=(3,4),delimiter=',',unpack=True)     
    time=fr[0]
    signal=fr[1]
    gate=square(2*pi*frep*(time-0.02e-9),duty=0.2)+1
    prod=gate*signal
    fxs,fftsig=normfft(signal)
    fxp,fftprod=normfft(prod)
    
    fft_sig.append(fftsig)
    fft_gat.append(fftprod)    

fft_sig = np.array(fft_sig)        
fft_gat = np.array(fft_gat)        

#%% Mean values 
freq = fxs
Pmed_sig = 10*np.log10(np.mean(np.abs(fft_sig)**2,axis=0)/50/1e-3)
Pmed_gat = 10*np.log10(np.mean(np.abs(fft_gat)**2,axis=0)/50/1e-3)

#%% Calculate peak and noise level
qq = np.where((fxs>80e6)&(fxs<100e6))
max_sig = np.max(Pmed_sig[qq])
max_gat = np.max(Pmed_gat[qq])

#%% Save file
#A = [freq,Pmed_sig-max_sig,Pmed_gat-max_gat]
#A = np.array(A)
#np.savetxt("gator_test.csv",A.T,delimiter=",",
#           header="fr[Hz],Pow Sig (dBm),Pow Gate (dBm)")

#%% Make a few plots
fig = plt.figure(87,figsize=(10,6))
ax = fig.add_subplot(111)
ax.plot(freq,Pmed_sig-max_sig,label='Ref signal')
ax.plot(freq,Pmed_gat-max_gat,label='Gated signal')
ax.set_xlim(0,350e6)
ax.set_ylim(-88,2)
ax.set_xlabel(r'Freq', fontsize=14)
ax.set_ylabel(r'Power (dBm)', fontsize=14)    
ax.grid(linestyle='--')
fig.tight_layout()

gain = calc_snr(freq,Pmed_gat,95e6)-calc_snr(freq,Pmed_sig,95e6)

#%% Load all the traces
time = np.array([])
data = []

for f in files:     
    my_data = loadtxt(path+f,usecols=(3,4),delimiter=',',unpack=True)     
    data.append(my_data[1])
    if len(time)==0:
        time = my_data[0]

data = np.array(data)

#%% Varying the delay
fxs,fftsig=normfft(data)
delay_array = np.linspace(-0.5e-9,0.5e-9,101)

for ii in delay_array:
    pwm = square(2*pi*frep*(time-ii),duty=0.1)+1
    prod=pwm*data
    fxp,fftprod=normfft(prod)
    
    
    
#%%    
    