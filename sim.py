####################################################
######  PSD: DSP internal vs DAC phase comparison ##
####################################################
from pylab import *

from scipy import signal
from scipy.optimize import curve_fit
from numpy import *
from scipy.signal import *
from scipy.fft import fft
# from mpl_toolkits.axes_grid.inset_locator import (inset_axes, InsetPosition,mark_inset)
# from matplotlib.ticker import FuncFormatter

def normfft(s):
    si=(2/len(s) *abs(fft(s)))**2/2
    out=si[0:len(signal)//2]
    frout= fftfreq(len(s), ts)[0:len(s)//2]
    return frout,out

frep=250e6
fbeat= 30e6
fs=10*2.5e9
ts=1/fs
time=arange(0,10000,1)*ts
v0=1
alpha=1
SNRg=180
SNRs=120

gate=(sin(2*pi*frep*time)>0.9999)
gate = (square(2*pi*frep*(time-1.99e-9),duty=0.1)+1)*0.5
# gate=v0*sin(2*pi*frep*time)
signal=gate*alpha*sin(2*pi*fbeat*time)

fig=plt.figure('timeseries3')
ax1=fig.add_subplot(3,1,1)
ax2=fig.add_subplot(3,1,2,sharex=ax1)
ax3=fig.add_subplot(3,1,3,sharex=ax1)

ax1.plot(time,gate,'gray')
ax2.plot(time,signal,'gray')

# ax1.set_xlim(0,400e-9)

noisegate=(random.rand(len(gate))-0.5)*2*sqrt(2)*sqrt(10**(-SNRg/10)*v0**2/2*fs/2)
noisesig=(random.rand(len(signal))-0.5)*2*sqrt(2)*sqrt(10**(-SNRs/10)*(v0*alpha/2)**2/2*fs/2)

# fg,psdgate=welch(gate,fs)


fxg,fftgate=normfft(gate)
fxs,fftsig=normfft(signal)

plt.figure('PSDs3')
# plt.plot(fg,psdgate/time[-1],'gray')
# plt.plot(fxg,fftgate,'black')
plt.yscale('log')


#for adding noise, keep this order, otherwise the signal gets also the noise of the gate
signal=signal+noisesig
gate=gate+noisegate


prod=gate*signal

fig=plt.figure('timeseries3')
ax1.plot(time,gate,'red')
ax2.plot(time,signal,'blue')
ax3.plot(time,prod,'green')

fxg,fftgate=normfft(gate)
fxs,fftsig=normfft(signal)
fxp,fftprod=normfft(prod)

expectedgate=10**(-SNRg/10)*(v0**2)/2*1/time[-1]
expectedsig=10**(-SNRs/10)*(v0*alpha/4)**2/2*1/time[-1]

plt.figure('PSDs3')
plt.ylim(1e-15,1)
plt.xlim(0,250e6)
plt.plot(fxg,fftgate,'red')
plt.hlines(expectedgate,0,fs/2,'red')

plt.hlines(v0**2/2,0,fs/2,'gray')
plt.plot(fxs,fftsig,'blue')

plt.plot(fxp,fftprod,'green')
plt.hlines(max(expectedsig,expectedgate)/(2*fs/2/frep),0,fs/2,'green') #improvement = 2B/frep; B=fs/2; there are 2 beats / frep

plt.hlines(expectedsig,0,fs/2,'blue')

# plt.yscale('log')



# ax4=fig.add_subplot(2,2,4,sharex=ax2,sharey=ax3)

# ax1.tick_params(direction='in',labelbottom=False)#,labelsize=26)
# ax2.tick_params(direction='in',labelbottom=False,labelleft=False)
# ax3.tick_params(direction='in')#,labelsize=26)
# ax4.tick_params(direction='in',labelleft=False)#,labelsize=26)
# ax1.text(-0.76,1.25,'a')#,fontsize=26)
# ax2.text(-0.18,1.25,'b')#,fontsize=26)
# ax3.text(-0.76,1.25,'c')#,fontsize=26)
# ax4.text(-0.18,1.25,'d')#,fontsize=26)
# # for A in  [ax1,ax3]:
#     # A.yaxis.set_major_locator(FixedLocator(levels))
#     # A.yaxis.set_major_formatter(FixedFormatter(labels))
# for A in [ax1,ax2,ax3,ax4]:
#     A.grid(axis='y')
#     A.set_ylim(-0.2,1.2)
# axseq=[ax1,ax1,ax2,ax2,ax3,ax4]
# col=['red','black','red','black','red','red']
# for F,axe,c in zip(fi,axseq,col):
    
#     f = load(F)
#     t=f['arr_0'][:,0]
#     y=f['arr_0'][:,1]
#     y=(y+1)/2
#     if F== 'trace_20201014_fringes_counts_11_4Hz_counter.npz' or F=='trace_20201014_fringes_counts_11_4Hz_PD.npz' or F=='trace_20201014_fringes_counts_attenuated_0_4Hz_counter.npz':
#         y = y[t<5]
#         t=t[t<5]
        
#     axe.plot(t,y,color=c)


# f=load('trace_20201014_fringes_counts_0_free_counter.npz')
# t=f['arr_0'][:,0]
# y=f['arr_0'][:,1]
# y=(y+1)/2
# plt.figure(figsize=(4,2))
# plt.plot(t,y,'red')
# f=load('trace_20201014_fringes_counts_0_free_PD.npz')
# t=f['arr_0'][:,0]
# y=f['arr_0'][:,1]
# y=(y+1)/2
# plt.plot(t,y,'black')
# plt.xlabel('Time / s')#,fontsize=26,labelpad=20)
# plt.ylabel('Normalised interference')#,fontsize=26, labelpad=100)
