#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 11 00:37:17 2022

@author: avramkachura
"""

import numpy as np
import matplotlib.pyplot as plt
import time  
from numpy import *
import math
from impedance.models.circuits import CustomCircuit
import h5py
import sys
from scipy.fft import fft, ifft
from numpy import sum,isrealobj,sqrt
from numpy.random import standard_normal

def ERROR_MODEL_SNR_v1(f, Zre, Zim, SNR):
    n = len(f)
    Z = Zre+Zim*1j
    
    Zout = np.zeros(n, dtype = "complex_")
    Zoutn = np.zeros(n, dtype = "complex_")
    
    for i in range(0, n):
        fi = f[i]
        
        #print('Testing f =', fi)
        #system buffer
        if fi > 244:
            fs = 250000 #sampling rate (samples/s) (Can't set higher)
            Ts = 1/fs #sampling time s/samples
            
            x1 = 0 
            x2 = Ts*(1024-1) + x1 #ending time  = Ts*(samples-1) + x1
            t1 = np.linspace(x1,x2,1024)
            N = len(t1) #should always be 1024
            
        elif fi < 244 and fi > 1:
            #to test this need to try frequency that dont give out factors to 250000 i.e. 100, 150,...
            #print x2 and 1/fi*cycles to see how close it is
            #if x2 > 1/fi*cycles fs is less than fsdes: capture more than indended cycles (nbest > nideal)
            #if x2 < 1/fi*cycles fs is more than fsdes: capture less than indended cycles (nbest > nideal)
            
            #nideal = 250000/fs
            #nbest = round(250000/fs)
            #if nbest > nideal: capture more than indended cycles (fsbest < fideal)
            #if nbest < nideal: capture less than indended cycles (fsbest > fideal)
            cycles = 3
            fs = fi/cycles * 1024 #desired fs = fsideal
            nd = round(250000/fs) #closest decimation factor = rounding to closest factor of fs = fsbest
            fs = 250000/nd #actual fs
            Ts = 1/fs
            
            x1 = 0
            x2 = Ts*(1024-1) + x1
            t1 = np.linspace(x1,x2,1024)
            N = len(t1)
            
        elif fi < 1:
            fs = fi * 1024 #desired fs
            nd = round(250000/fs) #closest decimation factor = rounding to closest factor of fs
            fs = 250000/nd
            Ts = 1/fs
            
            x1 = 0
            x2 = Ts*(1024-1) + x1
            t1 = np.linspace(x1,x2,1024)
            N = len(t1)
        
        else: #fi == 244Hz
            fs = 250000 #sampling rate (samples/s) (Can't set higher)
            Ts = 1/fs
            
            x1 = 0
            x2 = Ts*(1024-1) + x1
            t1 = np.linspace(x1,x2,1024)
            N = len(t1)
        
        #Current I(t)
        x = sin(2*pi*fi*t1)
        
        #Voltage V(t) = I(t)*Z(fi)
        y = x*Z[i]
        
        #Impedance Z(fi) = V(fi)/I(fi), no noise
        Zout[i] = fi_FFT(y, fi, fs, N)/fi_FFT(x, fi, fs, N)
        
        #Noise
        #Current In(t)
        xn = awgn(x, SNR)
        
        #Voltage Vn(t)
        yn = awgn(y, SNR)
        
        #Impedance Z(fi) = V(fi)/I(fi), with noise
        Zoutn[i] = fi_FFT(yn, fi, fs, N)/fi_FFT(xn, fi, fs, N)
        
    return Zoutn

#FFT function
def fi_FFT(x, fi, fs, N): #signal, signal freq, sampling freq, number of samples
    y = fft(x)
    P2 = y/N
    P1 = P2[0:N//2+1] #// makes it an int
    P1[1:-2] = 2*P1[1:-2] #-2 is second last element
    fbins = fs*np.arange(0,N/2+1)/N
    n = np.argmin(np.abs(fbins-fi))
    return P1[n]

def awgn(s,SNRdB,L=1):
    """
    AWGN channel
    Add AWGN noise to input signal. The function adds AWGN noise vector to signal 's' to generate a resulting signal vector 'r' of specified SNR in dB. It also
    returns the noise vector 'n' that is added to the signal 's' and the power spectral density N0 of noise added
    Parameters:
        s : input/transmitted signal vector
        SNRdB : desired signal to noise ratio (expressed in dB) for the received signal
        L : oversampling factor (applicable for waveform simulation) default L = 1.
    Returns:
        r : received signal vector (r=s+n)
"""
    gamma = 10**(SNRdB/10) #SNR to linear scale
    if s.ndim==1:# if s is single dimensional vector
        P=L*sum(abs(s)**2)/len(s) #Actual power in the vector
    else: # multi-dimensional signals like MFSK
        P=L*sum(sum(abs(s)**2))/len(s) # if s is a matrix [MxN]
    N0=P/gamma # Find the noise spectral density
    if isrealobj(s):# check if input is real/complex object type
        n = sqrt(N0/2)*standard_normal(s.shape) # computed noise
    else:
        n = sqrt(N0/2)*(standard_normal(s.shape)+1j*standard_normal(s.shape))
    r = s + n # received signal
    return r    
