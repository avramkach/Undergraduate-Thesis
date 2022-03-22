#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 11 13:55:42 2022

@author: avramkachura

V1 to evaluate function_v1 with outdated plots
"""

import numpy as np
import matplotlib.pyplot as plt
import time  
from numpy import *
import math
from impedance.models.circuits import CustomCircuit
import sys
import h5py
import matplotlib as mpl
mpl.rcParams['figure.dpi'] = 200 #set figure to higher default dpi

#Import function in different directory
sys.path.insert(1, '/Users/avramkachura/Desktop/Fourth Year/ESC499/ECM_Extraction/ECM_EXT_v2.py') 

from ECM_EXT_v2 import *

sys.path.insert(1, '/Users/avramkachura/Desktop/Fourth Year/ESC499/Python Error Model/ERROR_MODEL_SNR_function_v1.py')

from ERROR_MODEL_SNR_function_v1 import *

'''sys.path.insert(1, '/Users/avramkachura/Desktop/Fourth Year/ESC499/Python Error Model/ERROR_MODEL_SNR_function_v2.py')

from ERROR_MODEL_SNR_function_v2 import *'''


def rmse(a, b):
    "Root Mean Square Error = root(|sum(a - b)|/n)"
    
    '''if isinstance(a, np.ndarray) and isinstance(b, np.ndarray):
        n = len(a)
    else:
        n = 1'''
        
    #Make a andb np arrays and flatten if arrays to make it work with array inputs
    a = np.array(a).flatten()
    b = np.array(b).flatten()
    n = len(a)
    return np.linalg.norm(a - b) / np.sqrt(n)

#circuit model, set param defaults to 0
def Zgen(fi, R0 = 0, R1 = 0, C1 = 0, sig = 0, L = 0): 
    '''L = 0
    sig = 0
    R0 = 0.5
    R1 = 1
    C1 = 1'''
    
    #Cell5_0 from EIS Analyzer fit
    R0 = 0.26424
    C1 = 0.42183
    R1 = 0.21616
    sig = 0.11111

    Zc = 1/(2*pi*fi*C1*1j)
    ZL = 1j*2*pi*fi*L;
    Zw = sig/(np.sqrt(2*pi*fi))*(1-1j) #Semifinite warburg
    
    #R + R//C circuit
    #Z = R0 + 1/(1/R1+1/Zc)
    
    #ECM EXT paper 2nd order randles model   
    #Z = R0 + 1/(1/(R1+Zw)+1/Zc) + ZL + 1/(1/0.0135778+2*pi*fi*0.0344249)
    
    #1st order randles model
    Z = R0 + 1/(1/(R1+Zw)+1/Zc)
    
    #n = len(fi)
    #Z = np.ones(n)
    return R0, R1, C1, Z

#test frequencies
f = np.array([1952,1708,1464,1220,976,732,488,
              244,122,61,30.50000,15.25000,7.62500,
              3.81250,1.90625,0.95300,0.47600,
              0.40000,0.30000,0.20000,0.10000])

#Evaluate RC elements for now
R0, R1, C1, Z = Zgen(f)
Zre = Z.real
Zim = Z.imag

#Need to add robustness to taking peak at arround 10%>
SNRarray = np.array([np.inf])
#SNRarray = np.linspace(1,100,100)
sn = len(SNRarray)

#Error Arrays % = |(estimate - actual)/actual| * 100%
R0ERR = np.zeros(sn)
R1ERR = np.zeros(sn)
C1ERR = np.zeros(sn)

countf = 3 #figure counter in loop

for i in range(0,sn):
    print('---------' + 'Testing SNR = ' + str(SNRarray[i]) + '---------')
    
    #Noise Model
    Zoutn = ERROR_MODEL_SNR_v1(f, Zre, Zim, SNRarray[i]) #v2 matches
    
    #ECM extract estimates for R0e
    R0e, R1e, C1e, ZoutECM = ECM_EXT_v2(f, Zoutn.real, Zoutn.imag)

    R0ERR[i] = np.abs((R0e-R0)/R0)*100
    R1ERR[i] = np.abs((R1e-R1)/R1)*100
    C1ERR[i] = np.abs((C1e-C1)/C1)*100

    plt.figure(1 + i*countf)
    plt.plot(Zre, -Zim, label = 'Actual Data', marker = 'o', color  = 'purple', linestyle = 'None')
    plt.plot(Zoutn.real, -Zoutn.imag, label = 'System Buffer with Noise', marker = 'x', color = 'orange', linestyle = 'None')
    plt.plot(ZoutECM.real, -ZoutECM.imag, label = 'ECM_EXT', color  = 'green', linestyle = '-')
    
    lmin = (diff(np.sign(diff(-Zoutn.imag))) > 0).nonzero()[0] + 1 # local min indices
    lmax = (diff(np.sign(diff(-Zoutn.imag))) < 0).nonzero()[0] + 1 # local max indices
    plt.plot(Zoutn.real[lmin],-Zoutn.imag[lmin], 'kv', label = 'local min')
    plt.plot(Zoutn.real[lmax],-Zoutn.imag[lmax], 'r^', label = 'local max')
    
    plt.title('Nyquist Plot: SNR = ' + str(SNRarray[i]))
    plt.xlabel('Zreal (\u03A9)')
    plt.ylabel('-Zim (\u03A9)')
    plt.grid(b=True, which='major', color='#999999', linestyle='-')
    plt.legend()
    
    plt.figure(2 + i*countf)
    plt.title('Bode Magnitude Plot: SNR = ' + str(SNRarray[i]))
    plt.plot(f, np.abs(Z), label = 'Actual Data', marker = 'o', color  = 'purple', linestyle = '-')
    plt.plot(f, np.abs(Zoutn), label = 'System Buffer with Noise', marker = 'x', color = 'orange', linestyle = '-')
    plt.plot(f, np.abs(ZoutECM), label = 'ECM_EXT', marker = 'D', color  = 'green', linestyle = '-')
    
    plt.xscale('log')
    plt.grid(b=True, which='major', color='#999999', linestyle='-')
    plt.xlabel('f (Hz)')
    plt.ylabel('|Z|')
    plt.legend()
    
    plt.figure(3 + i*countf)
    plt.title('Bode Magnitude Plot: SNR = ' + str(SNRarray[i]))
    plt.plot(f, np.angle(Z), label = 'Actual Data', marker = 'o', color  = 'purple', linestyle = '-')
    plt.plot(f, np.angle(Zoutn), label = 'System Buffer with Noise', marker = 'x', color = 'orange', linestyle = '-')
    plt.plot(f, np.angle(ZoutECM), label = 'ECM_EXT', marker = 'D', color  = 'green', linestyle = '-')
    plt.xscale('log')
    plt.grid(b=True, which='major', color='#999999', linestyle='-')
    plt.xlabel('f (Hz)')
    plt.ylabel('\u03C6 (rad)')
    plt.legend()
    
'''plt.figure(4 + (sn-1)*countf)
plt.title('ECM_EXT Parameter Error for Varying SNR with 1st-order Randles Model')
plt.plot(SNRarray, R0ERR, label = 'R0')
plt.plot(SNRarray, R1ERR, label = 'R1')
plt.plot(SNRarray, C1ERR, label = 'C1')
plt.xlabel('SNR')
plt.ylabel('% Error')
plt.legend()'''