# -*- coding: utf-8 -*-
"""
Created on Fri Apr 22 09:30:36 2022

@author: Christian Escobar
Code adapted from: 
https://pythonnumericalmethods.berkeley.edu/notebooks/chapter24.04-FFT-in-Python.html
This notebook contains an excerpt from the Python Programming and Numerical 
Methods - A Guide for Engineers and Scientists, the content is also available
 at Berkeley Python Numerical Methods.

The copyright of the book belongs to Elsevier. We also have this interactive
 book online for a better learning experience. The code is released under the 
 MIT license. If you find this content useful, please consider supporting the
 work on Elsevier or Amazon!

"""

import matplotlib.pyplot as plt
import numpy as np
from numpy.fft import fft
from scipy.signal import chirp
from scipy.signal import butter, lfilter, freqz

## TIME DOMAIN

# Create the figure. 
plt.figure(figsize = (8, 15)) # figsize = (width, height)

# Sampling rate
fs = 2000
# Sampling interval
ts = 1/fs
t = np.arange(0,1,ts) # Create the time array from 0 to 0.05 seconds

# Signal 1
freq = 15 # Frequency in Hz
x1 =np.cos(2*np.pi*freq*t) # cos(2*pi*f*t)
plt.subplot(4, 1, 1)
plt.plot(t, x1, 'b')
plt.ylabel('Amplitude')
plt.title('x1 = cos('+str(freq)+'*2$\pi$t)')
plt.xlim(0, .25) # Change xlim for display purposes

# Signal 2
freq = 450
x2= np.cos(2*np.pi*freq*t) # cos(2*pi*f*t)
plt.subplot(4, 1, 2)
plt.plot(t, x2, 'b')
plt.ylabel('Amplitude')
plt.title('x2 = cos('+str(freq)+'*2$\pi$t)')
plt.xlim(0, .25) # Change xlim for display purposes

freq = 300  
x3 = chirp(t, f0=80, f1=350, t1=1, method='linear')
#x3= np.cos(2*np.pi*freq*t) # cos(2*pi*f*t)
plt.subplot(4, 1, 3)
plt.plot(t, x3, 'b')
plt.ylabel('Amplitude')
plt.title('x3 = Chirp: 80 - 350 Hz')
plt.xlim(0, .25) # Change xlim for display purposes

x4= x1+x2+x3
plt.subplot(4, 1, 4)
plt.plot(t, x4, 'b')
plt.ylabel('Amplitude')
plt.title('x4 = x1 + x2 + x3')
plt.xlim(0, .25) # Change xlim for display purposes

# plt.show()

## FREQUENCY DOMAIN
#Signal 1
y1 = fft(x1)
N = len(y1)
n = np.arange(N)
T = N/fs
freq = n/T 
freq=np.arange(0,(fs/2),(fs/N))
P2 = np.abs(y1/N);
P1 = P2[1:int(N/2+1)];
P1[2:-1] = 2*P1[2:-1];

# Create the figure. 
plt.figure(figsize = (8, 15)) # figsize = (width, height)
plt.subplot(4, 1, 1)
# plt.subplot(4, 1, 1)
plt.plot(freq, P1, 'b')
plt.ylabel('FFT Amplitude |X(freq)|')
plt.xlim(0, 1000) # Change xlim for display purposes
plt.grid()
plt.title('fft of: x1')


#Signal 2
y2 = fft(x2)
N = len(y2)
n = np.arange(N)
T = N/fs
freq = n/T 
freq=np.arange(0,(fs/2),(fs/N))
P2 = np.abs(y2/N);
P1 = P2[1:int(N/2+1)];
P1[2:-1] = 2*P1[2:-1];

# Create the figure. 
# plt.figure(figsize = (8, 15)) # figsize = (width, height)
plt.subplot(4, 1, 2)
plt.subplot(4, 1, 2)
plt.plot(freq, P1, 'b')
plt.ylabel('FFT Amplitude |X(freq)|')
plt.xlim(0, 1000) # Change xlim for display purposes
plt.grid()
plt.title('fft of: x2')


#Signal 3
y3 = fft(x3)
N = len(y3)
n = np.arange(N)
T = N/fs
freq = n/T 
freq=np.arange(0,(fs/2),(fs/N))
P2 = np.abs(y3/N);
P1 = P2[1:int(N/2+1)];
P1[2:-1] = 2*P1[2:-1];

# Create the figure. 
# plt.figure(figsize = (8, 15)) # figsize = (width, height)
plt.subplot(4, 1, 3)
plt.subplot(4, 1, 3)
plt.plot(freq, P1, 'b')
plt.ylabel('FFT Amplitude |X(freq)|')
plt.xlim(0, 1000) # Change xlim for display purposes
plt.grid()
plt.title('fft of: x3')


#Signal 4
y4 = fft(x4)
N = len(y4)
n = np.arange(N)
T = N/fs
freq = n/T 
freq=np.arange(0,(fs/2),(fs/N))
P2 = np.abs(y4/N);
P1 = P2[1:int(N/2+1)];
P1[2:-1] = 2*P1[2:-1];

# Create the figure. 
plt.subplot(4, 1, 4)
plt.plot(freq, P1, 'b')
plt.ylabel('FFT Amplitude |X(freq)|')
plt.xlim(0, 1000) # Change xlim for display purposes
plt.grid()
plt.title('fft of: x4')


def butter_bandstop(lowcut,highcut, fs, order=5):
    nyq = 0.5 * fs
    low = lowcut / nyq
    high = highcut / nyq
    b, a = butter(order, [low, high], btype='bandstop')#, analog=False)
    return b, a

def butter_bandstop_filter(data, lowcut,highcut, fs, order=4):
    b, a = butter_bandstop(lowcut,highcut, fs, order=order)
    y = lfilter(b, a, data)
    return y


lowcut=15
highcut= 400
order=5
filtered_sine = butter_bandstop_filter(x4,lowcut,highcut,fs,order)

plt.figure(figsize = (8, 15)) # figsize = (width, height)
plt.subplot(3, 1, 3)
plt.plot(t, filtered_sine, 'b')
plt.plot(t, x4, 'k')

plt.xlim(0, .1) # Change xlim for display purposes

b, a = butter_bandstop(lowcut,highcut, fs, order)

# Plotting the frequency response.
w, h = freqz(b, a, worN=1000)

plt.subplot(3, 1, 1)
plt.plot(0.5*fs*w/np.pi, np.abs(h), 'b')
plt.xlim(0, 0.5*fs)
plt.title("Lowpass Filter Frequency Response")
plt.xlabel('Frequency [Hz]')
plt.grid()
plt.show()
