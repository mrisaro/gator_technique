####################################################
######  PSD: DSP internal vs DAC phase comparison ##
####################################################
from pylab import *

from scipy import signal
from scipy.optimize import curve_fit
from numpy import *
from scipy.signal import *
from scipy.fft import fft
from pandas import *
# from mpl_toolkits.axes_grid.inset_locator import (inset_axes, InsetPosition,mark_inset)
# from matplotlib.ticker import FuncFormatter

def normfft(s):
    si=(2/len(s) *abs(fft(s)))**2/2 # multiply by 2 to account for negative freqs, divide by two to account for effective value
    out=si[0:len(signal)//2]
    frout= fftfreq(len(s), ts)[0:len(s)//2]
    return frout,out

fr= loadtxt('gator_data_raw/210316_140348_Ch1.csv',usecols=(3,4),delimiter=',',unpack=True)
frep=250e6
ts=40e-12
fs=1/ts

time=fr[0]
signal=fr[1]
gate=square(2*pi*frep*(time-0.02e-9),duty=0.03)+1
prod=gate*signal

# frep=250e6
# fbeat= 30e6
# fs=10*2.5e9
# B=2e9
# ts=1/fs
# time=arange(0,10000,1)*ts
# v0=1
# alpha=1
# SNRg=180
# SNRs=140


# gate=(sin(2*pi*frep*time)>0.9999)*v0
# signal=gate*alpha*sin(2*pi*fbeat*time)

# noisegate=(random.rand(len(gate))-0.5)*2*sqrt(2)*sqrt(10**(-SNRg/10)*v0**2/2*fs/2)
# noisesig=(random.rand(len(signal))-0.5)*2*sqrt(2)*sqrt(10**(-SNRs/10)*(v0*alpha/2)**2/2*fs/2)
#     #for adding noise, keep this order, otherwise the signal gets also the noise of the gate
# signal=signal+noisesig
# gate=gate+noisegate
# prod=gate*signal

# fxg,fftgate=normfft(gate)
fxs,fftsig=normfft(signal)
fxp,fftprod=normfft(prod)
    
fig=plt.figure('timeseries3')
ax1=fig.add_subplot(3,1,1)
ax2=fig.add_subplot(3,1,2,sharex=ax1)
ax3=fig.add_subplot(3,1,3,sharex=ax1)
ax1.plot(time*1e9,gate,color='red')
ax2.plot(time*1e9,signal,color='gray')
ax3.plot(time*1e9,prod,color='red')

plt.figure('PSDs_alone')
plt.loglog(fxp,fftsig/max(fftsig[((fxp>90e6)&(fxp<100e6))]),color='gray')

plt.figure('PSDsP')
plt.loglog(fxp,fftsig/max(fftsig[((fxp>90e6)&(fxp<100e6))]),color='gray')
plt.loglog(fxp,fftprod/max(fftprod[((fxp>90e6)&(fxp<100e6))]),color='red')
plt.grid(True,'both')

   


    ########filtered version #########
c='green'
B=2e9
gate=square(2*pi*frep*(time+0.15e-9),duty=0.03)+1

b, a = butter(3, B/(fs/2), btype='low', analog = False)
gate=lfilter(b, a, gate)

prod=abs(gate)*signal

fxp,fftprod=normfft(prod)
    
fig=plt.figure('timeseries3')

ax1.plot(time*1e9,gate,color=c)
ax3.plot(time*1e9,prod,color=c)
    
plt.figure('PSDs_alone')
plt.loglog(fxp,fftprod/max(fftprod[((fxp>90e6)&(fxp<100e6))]),color=c)
 
#     ########square-wave version ########
c='orange'
# gsi=(sin(2*pi*frep*time)>0.9999)*v0
# b, a = butter(3, B/(fs/2), btype='low', analog = False)
# gsi=lfilter(b, a, gsi)
# signal=gsi*alpha*sin(2*pi*fbeat*time)
gate=square(2*pi*frep*(time+0.1e-9),duty=0.1)+1

# signal=signal+noisesig
# gate=gate+noisegate
prod=gate*signal

fxg,fftgate=normfft(gate)
fxs,fftsig=normfft(signal)
fxp,fftprod=normfft(prod)
    
fig=plt.figure('timeseries3')
# ax1=fig.add_subplot(3,1,1)
# ax2=fig.add_subplot(3,1,2,sharex=ax1)
# ax3=fig.add_subplot(3,1,3,sharex=ax1)
ax1.plot(time*1e9,gate,color=c)
# ax2.plot(time*1e9,signal,color=c)
ax3.plot(time*1e9,prod,color=c)
    
plt.figure('PSDsP')
plt.loglog(fxp,fftprod/max(fftprod[((fxp>90e6)&(fxp<100e6))]),color=c)
plt.grid(True,'both')

#     ########square-wave version, narrow ########
# c='green'
# gsi=(sin(2*pi*frep*time)>0.9999)*v0
# b, a = butter(3, B/(fs/2), btype='low', analog = False)
# gsi=lfilter(b, a, gsi)
# signal=gsi*alpha*sin(2*pi*fbeat*time)
# gate=square(2*pi*frep*(time-1.1e-9),duty=0.05)*v0+v0

# signal=signal+noisesig
# gate=gate+noisegate
# prod=gate*signal

# fxg,fftgate=normfft(gate)
# fxs,fftsig=normfft(signal)
# fxp,fftprod=normfft(prod)
    
# fig=plt.figure('timeseries3')
# ax1=fig.add_subplot(3,1,1)
# ax2=fig.add_subplot(3,1,2,sharex=ax1)
# ax3=fig.add_subplot(3,1,3,sharex=ax1)
# ax1.plot(time,gate,color=c)
# ax2.plot(time,signal,color=c)
# ax3.plot(time,prod,color=c)
    
# plt.figure('PSDsP')
# plt.loglog(fxp,fftprod/max(fftprod[((fxp>28e6)&(fxp<32e6))]),color=c)
# plt.grid(True,'both')

# c='blue'
# gsi=(sin(2*pi*frep*time)>0.9999)*v0
# b, a = butter(3, B/(fs/2), btype='low', analog = False)
# gsi=lfilter(b, a, gsi)
# signal=gsi*alpha*sin(2*pi*fbeat*time)
# gate=square(2*pi*frep*(time-1.15e-9),duty=0.005)*v0+v0

# signal=signal+noisesig
# gate=gate+noisegate
# prod=gate*signal

# fxg,fftgate=normfft(gate)
# fxs,fftsig=normfft(signal)
# fxp,fftprod=normfft(prod)
    
# fig=plt.figure('timeseries3')
# ax1=fig.add_subplot(3,1,1)
# ax2=fig.add_subplot(3,1,2,sharex=ax1)
# ax3=fig.add_subplot(3,1,3,sharex=ax1)
# ax1.plot(time,gate,color=c)
# ax2.plot(time,signal,color=c)
# ax3.plot(time,prod,color=c)
    
# plt.figure('PSDsP')
# plt.loglog(fxp,fftprod/max(fftprod[((fxp>28e6)&(fxp<32e6))]),color=c)
# plt.grid(True,'both')