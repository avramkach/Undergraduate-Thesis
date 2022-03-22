#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 11 00:37:17 2022

@author: avramkachura

Function Inputs: (ordered with corresponding impedances)
    -f (Hz)
    -Zre (+) (same units as Zim, ohm)
    -Zim (-) needs -Zim as input
(Error inputs:)
    -SNR (dB) 
    -td (s) (constant time delay between I(t) and V(t))
    
Outputs:
    Zout if SNR is inf (a+bj)
    Zoutn for SNR (a+bj)
    
#Current: 1.5mA each bit represents 1.5mA
5/1024 (5V/1024 bits how much), analog amplification factor

#Voltage: LSB 0.763 uV after enob

#use 5 or 6 noise bits -> fixed amplitude randomly sampled noise
#nyquist with mean and stdev
#until stdev stops changing 
#check results

#effective number of bit enob = 16

#

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

def ERROR_MODEL_v4(f, Zre, Zim, SNR, td):
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
            
            #Current samples
            x1I = 0 
            x2I = Ts*(1024-1) + x1I #ending time  = Ts*(samples-1) + x1
            tI = np.linspace(x1I,x2I,1024)
            N = len(tI) #should always be 1024
            
            #Voltage samples
            x1V = x1I + td
            x2V = Ts*(1024-1) + x1V
            tV = np.linspace(x1V,x2V,1024)
            
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
            
            #Current samples
            x1I = 0
            x2I = Ts*(1024-1) + x1I
            tI = np.linspace(x1I,x2I,1024)
            N = len(tI)
            
            #Voltage samples
            x1V = x1I + td
            x2V = Ts*(1024-1) + x1V
            tV = np.linspace(x1V,x2V,1024)
            
        elif fi < 1:
            fs = fi * 1024 #desired fs
            nd = round(250000/fs) #closest decimation factor = rounding to closest factor of fs
            fs = 250000/nd
            Ts = 1/fs
            
            #Current samples
            x1I = 0
            x2I = Ts*(1024-1) + x1I
            tI = np.linspace(x1I,x2I,1024)
            N = len(tI)
            
            #Voltage samples
            x1V = x1I + td
            x2V = Ts*(1024-1) + x1V
            tV = np.linspace(x1V,x2V,1024)
        
        else: #fi == 244Hz
            fs = 250000 #sampling rate (samples/s) (Can't set higher)
            Ts = 1/fs
            
            #Current samples
            x1I = 0
            x2I = Ts*(1024-1) + x1I
            tI = np.linspace(x1I,x2I,1024)
            N = len(tI)
            
            #Voltage samples
            x1V = x1I + td
            x2V = Ts*(1024-1) + x1V
            tV = np.linspace(x1V,x2V,1024)
        
        #Current I(t)
        Ishift = 0
        x = sin_wave(fi, tI, Ishift)
        
        #Voltage V(t)
        #V|<V = Io/root(2)<I * |Z|<Z = Io/root(2)|Z|<I+Z
        #V(t) = Io|Z|sin(wt + <I + <Z)
        #V(t) = |Z|sin(wt + <Z) #assuming Io = 1, <I = 0
        Vshift = 0 #might use this as a variable to test error -> good simple test
        #adding +shift will cause ECM not to work because it filters out Zim>0 -> need to make robust
        y = sin_wave(fi, tV, np.angle(Z[i]) + Ishift + Vshift, np.abs(Z[i]))
        
        #Impedance Z(fi) = V(fi)/I(fi), no noise
        Zout[i] = fi_FFT(y, fi, fs, N)/fi_FFT(x, fi, fs, N)
    
        #return cycles and SNR without noise think of the build needs to be inside loop
        if SNR == np.inf:
            #outputs Zout
            Zoutn[i] = Zout[i]
        else:
            #Noise with ENOB -> 
            
            #Current In(t)
            xn = awgn(x, SNR)
            
            #use measured here
            
            #Voltage Vn(t)
            yn = awgn(y, SNR)
            
            #Impedance Zn(fi) = Vn(fi)/In(fi), with noise
            Zoutn[i] = fi_FFT(yn, fi, fs, N)/fi_FFT(xn, fi, fs, N)
            
    return Zoutn
    
#sinusoidal generator
def sin_wave(f, t1, phi, A = 1):
    out = A*np.sin(2*np.pi*f*t1 + phi)
    return out

#FFT function
def fi_FFT(x, fi, fs, N): #signal, signal freq, sampling freq, number of samples
    y = fft(x)
    P2 = y/N
    P1 = P2[0:N//2+1] #// makes it an int
    P1[1:-2] = 2*P1[1:-2] #-2 is second last element
    fbins = fs*np.arange(0,N/2+1)/N
    n = np.argmin(np.abs(fbins-fi))
    return P1[n]

#Measured function: 
def measured(s, A, SNRdB):
    #computes ENOB and LSB based using equation in ADS1675
    #Same SNRdB for voltage and current measurements
    
    #Inputs: signal, amplitude (A), SNRdB
    #Outputs: signal measurements
    
    #A = ampltude of signal
    
    #RMS signal = 2*Amplitude (peak)/(root(2)) 
    RMSsig = 2*A/np.sqrt(2)
    
    ##RMS noise = RMS signal/((10^(SNR(dB)/20))**0.5)
    RMSnoise = RMSsig/((10**(SNRdB/20))**0.5)
    
    #ENOB = ln(6V/(RMS noise))/(ln2)
    ENOB = np.log(6/RMSnoise)/np.log(2)

    #LSB = Vref/(2^(ENOB)) = 3V/(2^(ENOB)) 
    LSB = 3/(2**ENOB)
    
    #initialized Smeasured to return 
    Smeas = np.zeros(len(s))
    
    for i in range(0,len(s)):
        sample = s[i]
        
        #take all sample that give values less than ENOB
        if np.floor(np.abs(sample/LSB)) == 0:
            Smeas[i] = 0
        #case if it goes above or equal to reference 3V/-3V
        elif np.abs(sample) >= 3:
            Smeas[i] = (3 - LSB)*np.sign(sample)
        #case if its within 3V/-3V range
        else:
            Smeas[i] = ((np.floor(np.abs(sample)/LSB))*LSB)*np.sign(sample)
    
    #include gain....
    
    return Smeas
    
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
