#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 11 13:55:42 2022

@author: avramkachura
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

sys.path.insert(1, '/Users/avramkachura/Desktop/Fourth Year/ESC499/Python Error Model/ERROR_MODEL_SNR_function_v2.py')

from ERROR_MODEL_SNR_function_v2 import *


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
def Zgen(f, R0 = 0, R1 = 0, C1 = 0, sig = 0, L = 0): 
    '''L = 0
    sig = 0
    R0 = 0.25
    R1 = 1
    C1 = 0.25'''
    
    #ECM EXT paper extracted from simulated    
    #L = 0.17549*10**(-6)
    '''R0 = 0.5509587
    R1 = 0.1193365
    C1 = 1.14637*
    sig = 0.0346208'''
    
    #Cell5_0 from EIS Analyzer fit
    R0 = 0.26424
    C1 = 0.42183
    R1 = 0.21616
    sig = 0.11111

    if C1 != 0:
        Zc = 1/(2*pi*f*C1*1j)
        
    ZL = 1j*2*pi*f*L;
    
    if sig != 0:
        Zw = sig/(np.sqrt(2*pi*f))*(1-1j) #Semifinite warburg
    
    #R + R//C circuit
    #Z = R0 + 1/(1/R1+1/Zc)
    
    #ECM EXT paper 2nd order randles model   
    #Z = R0 + 1/(1/(R1+Zw)+1/Zc) + ZL + 1/(1/0.0135778+2*pi*fi*0.0344249)
    
    #1st order randles model Cell5_0 EIS Analyzer fit
    Z = R0 + 1/(1/(R1+Zw)+1/Zc)
    
    '''R0 = 0.1
    n = len(f)
    Z = np.ones(n)'''
    
    return R0, R1, C1, Z

#add frequency labels for noise plotting with frequency cutoff
#https://queirozf.com/entries/add-labels-and-text-to-matplotlib-plots-annotation-examples
def addlabels(f, Zre, Zim):
    fcutoff = 10
    inval = np.where(f > fcutoff)
    f = np.delete(f, inval[0])
    ys = np.delete(-Zim, inval[0])
    xs = np.delete(Zre, inval[0])
    n = len(f)
    for i in range(0, n):
        
        label = "{:.2f}".format(f[i])
        
        if i % 2 == 0:
            plt.annotate(label, # this is the text
                         (xs[i],-ys[i]), # these are the coordinates to position the label
                         textcoords="offset points", # how to position the text
                         xytext=(0,5), # distance from text to points (x,y)
                         ha='center') # horizontal alignment can be left, right or center
        else:
            plt.annotate(label, # this is the text
                         (xs[i],-ys[i]), # these are the coordinates to position the label
                         textcoords="offset points", # how to position the text
                         xytext=(0,-12), # distance from text to points (x,y)
                         ha='center') # horizontal alignment can be left, right or center


#test frequencies
'''f = np.array([1952,1708,1464,1220,976,732,488,
              244, 122, 61,30.50000,15.25000,7.62500,
              3.81250,1.90625,0.95300,0.47600,
              0.40000,0.30000,0.20000,0.10000]) #dont see distortion with 200 and 100Hz, too small'''

f = np.array([1952,1708,1464,1220,976,732,488,
              244, 122, 61,30.50000,15.25000,7.62500,
              3.81250,1.90625, 1.745, 0.95300,0.47600,
              0.40000,0.30000,0.20000,0.10000]) #1.745 should be peak

#Evaluate RC elements for now
R0, R1, C1, Z = Zgen(f)
Zre = Z.real
Zim = Z.imag

#Need to add robustness to taking peak at arround 10%>, input is SNRdB
#SNRarray = np.array([np.inf]) #SNR = inf disables SNR
SNRarray = np.linspace(100,1,100) #start at higher SNR (low noise) to lower SNR (high noise)

#SNRarray = np.array([0.5, 0.1])

sn = len(SNRarray)

#Error Arrays % = |(estimate - actual)/actual| * 100%
R0ERR = np.zeros(sn)
R1ERR = np.zeros(sn)
C1ERR = np.zeros(sn)

countf = 3 #figure counter in loop

for i in range(0,sn):
    print('---------' + 'Testing SNR = ' + str(SNRarray[i]) + '---------')
    
    #Noise Model
    Zoutn = ERROR_MODEL_SNR_v2(f, Zre, Zim, SNRarray[i]) 
    
    #ECM extract estimates for R0e
    R0e, R1e, C1e, ZoutECM = ECM_EXT_v2(f, Zoutn.real, Zoutn.imag)

    R0ERR[i] = np.abs((R0e-R0)/R0)*100
    
    if R1 != 0 and C1 != 0:
        R1ERR[i] = np.abs((R1e-R1)/R1)*100
        C1ERR[i] = np.abs((C1e-C1)/C1)*100

    if sn <= 3: #length check to auto do it on SNR
        plt.figure(1 + i*countf)
        plt.plot(Zre, -Zim, label = 'Actual Data', marker = 'o', color  = 'purple', linestyle = 'None')
        plt.plot(Zoutn.real, -Zoutn.imag, label = 'System Buffer with SNR = ' + str(SNRarray[i]), marker = 'x', color = 'orange', linestyle = '-.')
        
        addlabels(f, Zoutn.real, -Zoutn.imag)
        
        plt.plot(ZoutECM.real, -ZoutECM.imag, label = 'ECM_EXT', color  = 'green', linestyle = '-')
        
        #lmin = (diff(np.sign(diff(-Zoutn.imag))) > 0).nonzero()[0] + 1 # local min indices
        lmin = (diff(np.sign(diff(-Zoutn.imag)/diff(Zoutn.real))) > 0).nonzero()[0] + 1
        #lmax = (diff(np.sign(diff(-Zoutn.imag))) < 0).nonzero()[0] + 1 # local max indices
        lmax = (diff(np.sign(diff(-Zoutn.imag)/diff(Zoutn.real))) < 0).nonzero()[0] + 1  #gets local mins? No
        
        plt.plot(Zoutn.real[lmin],-Zoutn.imag[lmin], 'kv', label = 'local min')
        plt.plot(Zoutn.real[lmax],-Zoutn.imag[lmax], 'r^', label = 'local max')
        
        plt.title('Nyquist Plot')
        plt.xlabel('Zreal (\u03A9)')
        plt.ylabel('-Zim (\u03A9)')
        plt.grid(b=True, which='both', linestyle='-')
        plt.legend()
        
        plt.figure(2 + i*countf)
        plt.title('Bode Magnitude Plot')
        plt.plot(f, np.abs(Z), label = 'Actual Data', marker = 'o', color  = 'purple', linestyle = '-')
        plt.plot(f, np.abs(Zoutn), label = 'System Buffer with SNR = ' + str(SNRarray[i]), marker = 'x', color = 'orange', linestyle = '-')
        plt.plot(f, np.abs(ZoutECM), label = 'ECM_EXT', marker = 'D', color  = 'green', linestyle = '-')
        
        plt.xscale('log')
        plt.grid(b=True, which='both', color='#DDDDDD', linestyle='-')
        plt.xlabel('f (Hz)')
        plt.ylabel('|Z| (\u03A9)')
        plt.legend()
        
        plt.figure(3 + i*countf)
        plt.title('Bode Magnitude Plot')
        plt.plot(f, np.angle(Z), label = 'Actual Data', marker = 'o', color  = 'purple', linestyle = '-')
        plt.plot(f, np.angle(Zoutn), label = 'System Buffer with SNR = ' + str(SNRarray[i]), marker = 'x', color = 'orange', linestyle = '-')
        plt.plot(f, np.angle(ZoutECM), label = 'ECM_EXT', marker = 'D', color  = 'green', linestyle = '-')
        plt.xscale('log')
        plt.grid(b=True, which='both', color='#DDDDDD', linestyle='-')
        plt.xlabel('f (Hz)')
        plt.ylabel('\u03C6 (rad)')
        plt.legend()
    
if sn > 1:
    plt.figure(4 + (sn-1)*countf)
    plt.title('ECM_EXT Parameter Error for Varying SNR with 1st-order Randles Model')
    
    #print Zoutn Zout residuals %, chisquar, rmse
    
    
    plt.plot(SNRarray, R0ERR, label = 'R0')
    plt.plot(SNRarray, R1ERR, label = 'R1')
    plt.plot(SNRarray, C1ERR, label = 'C1')
    plt.xlabel('SNR (dB)')
    #plt.xscale('log')
    #should be log scale... since its db or not can just consider it as is
    plt.ylabel('% Error')
    plt.legend()