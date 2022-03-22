#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 11 13:55:42 2022

@author: avramkachura

To test tdelay for now
"""

import numpy as np
import matplotlib.pyplot as plt
import time  
from numpy import *
import math
from impedance.models.circuits import CustomCircuit
import sys
import h5py
from sklearn.metrics import r2_score #r^2 value
from tabulate import tabulate #printing table

import matplotlib as mpl
mpl.rcParams['figure.dpi'] = 200 #set figure to higher default dpi

#Import function in different directory needs to be folder path not file
sys.path.insert(1, '/Users/avramkachura/Desktop/Fourth Year/ESC499/ECM_Extraction/') 

from ECM_EXT_v2 import *

from ECM_EXT_v3 import *

from ECM_EXT_C_v1 import *

sys.path.insert(1, '/Users/avramkachura/Desktop/Fourth Year/ESC499/Python Error Model/')

from ERROR_MODEL_function_v4 import *


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

def RSS(a):
    "Residual Sum of Squares = sum(a^2), a = residual = yi - f(xi)"
    return np.sum(a**2)

#circuit model, set param defaults to 0
def Zgen(f, R0 = 0, R1 = 0, C1 = 0, sig = 0, L = 0): 
    '''L = 0
    sig = 0
    R0 = 0.25
    R1 = 1
    C1 = 0.25'''
    
    #ECM EXT paper extracted from simulated    
    L = 0.17549*10**(-6)
    R0 = 0.5509587
    R1 = 0.1193365
    C1 = 1.14637
    sig = 0.0346208
    
    #Cell5_0 from EIS Analyzer fit
    '''R0 = 0.26424
    C1 = 0.42183
    R1 = 0.21616
    sig = 0.11111'''

    if C1 != 0:
        Zc = 1/(2*pi*f*C1*1j)
        
    ZL = 1j*2*pi*f*L;
    
    if sig != 0:
        Zw = sig/(np.sqrt(2*pi*f))*(1-1j) #Semifinite warburg
    
    #R + R//C circuit
    #Z = R0 + 1/(1/R1+1/Zc)
    
    #ECM EXT paper 2nd order randles model   
    Z = R0 + 1/(1/(R1+Zw)+1/Zc) + ZL + 1/(1/0.0135778+2*pi*f*0.0344249)
    
    #1st order randles model Cell5_0 EIS Analyzer fit
    #Z = R0 + 1/(1/(R1+Zw)+1/Zc)
    
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
f = np.array([1952,1708,1464,1220,976,732,488,
              244, 122, 61,30.50000,15.25000,7.62500,
              3.81250,1.90625,0.95300,0.47600,
              0.40000,0.30000,0.20000,0.10000]) #dont see distortion with 200 and 100Hz, too small

'''f = np.array([1952,1708,1464,1220,976,732,488,
              244, 122, 61,30.50000,15.25000,7.62500,
              3.81250,1.90625, 1.745, 0.95300,0.47600,
              0.40000,0.30000,0.20000,0.10000]) #1.745 should be peak'''

#Evaluate RC elements for now
R0, R1, C1, Z = Zgen(f)
Zre = Z.real
Zim = Z.imag


#sn = len(SNRarray)
#0.00001 = 10us = 1e-5
#tdarray = np.array([0.000001, 0.00001, 0.0001, 0.0005, 0.001, 0.005]) # tests shown zhe, change to us and print at top of plot
tdarray = np.array([0.00000]) #10us

#tdarray = np.array([0.00005]) #50us

#tdarray = np.array([0.001]) #at this frequency no inductive element robust to this change

#tdarray = np.array([0.01]) #10ms really breaks end at 1e-2
#tdarray = np.array([0.00000001]) #0.01us fimax = 1952, T = 0.0051229 /100 = 5us of period
#tdarray = np.array([0.0000001]) #0.1us start at 1e-7

#tdarray = np.logspace(-8,-2,100)

#tdarray = np.linspace(1e-7, 1e-2, 100)

#tdarray = np.linspace(1e-8, 1e-2, 100) #0.01us to 10ms

#tdarray = np.append([0], np.linspace(1e-8, 1e-2, 100))
#tdarray = np.append([0], np.logspace(-8,-2,100))

#breaks down at 15us and 50us delay looking at the imaginary impedance statistics
#tdarray = np.array([0.000015, 0.000050])

#tdarray = np.array([0, 0.000001, 0.00001, 0.0001, 0.0005, 0.001, 0.005])

#tdarray = np.array([0, 0.000001, 0.00001, 0.00005, 0.0001, 0.0005, 0.001])

#tdarray = np.linspace(1e-7, 5e-4, 100)

#tdarray = np.array([0.0004]) #400us spike

#tdarray = np.array([0])

tn = len(tdarray)

#Disable SNR by setting to inf same length as tn
SNRarray = np.ones(tn)*np.inf
#SNRarray = np.ones(tn)*1

#Param Error Arrays % = |(estimate - actual)/actual| * 100%
R0ERR = np.zeros(tn)
R1ERR = np.zeros(tn)
C1ERR = np.zeros(tn)

#Residual sum of squares RSS = sum(residual^2)
ZreRSS_M = np.zeros(tn)
ZimRSS_M = np.zeros(tn)

#RMSE Arrays M
ZreRMSE_M = np.zeros(tn)
ZimRMSE_M = np.zeros(tn)
MagRMSE_M = np.zeros(tn)
PhiRMSE_M = np.zeros(tn)


#R^2 arrays of coefficient of determination, degree of any linear correlation between y_true and y_pred, r2_score(y_true, y_pred)
ZreR2_M = np.zeros(tn)
ZimR2_M = np.zeros(tn)


'''#Chi-squared measure of the difference between the observed and actual values, implement later
ZreChi_M = np.zeros(tn)
ZimChi_M = np.zeros(tn)'''



countf = 5 #figure counter in loop + 2 for residual

for i in range(0,tn):
    print('-----' + 'Testing Time Delay (us) = ' + str(tdarray[i]*10**6) + '-----')
    print('-----' + 'Testing Noise (dB) = ' + str(SNRarray[i]) + '-----')
    
    #us conversion for printing
    tdarrayu = tdarray[i]*10**6
    
    #Noise Model
    Zoutn = ERROR_MODEL_v4(f, Zre, Zim, SNRarray[i], tdarray[i]) 
    
    #ECM extract estimates for R0e
    R0e, R1e, C1e, ZoutECM = ECM_EXT_C_v1(f, Zoutn.real, Zoutn.imag)

    #Calculate Errors
    R0ERR[i] = np.abs((R0e-R0)/R0)*100
    if R1 != 0 and C1 != 0: #R1 = 0 iff C1 = 0 where Zim peak isnt at 0
        R1ERR[i] = np.abs((R1e-R1)/R1)*100
        C1ERR[i] = np.abs((C1e-C1)/C1)*100
        
    #Calculate Residuals of Measured Z (Residual=actual−predicted = Actual - Measured)
    reResid_M = Zre - Zoutn.real
    imResid_M = Zim - Zoutn.imag
    
    #Calculated RSS of Measured Z
    ZreRSS_M[i] = RSS(reResid_M)
    ZimRSS_M[i] = RSS(imResid_M)
    
    #Calculae RMSE of Measured Z
    ZreRMSE_M[i] = rmse(Zre, Zoutn.real)
    ZimRMSE_M[i] = rmse(Zim, Zoutn.imag)
    MagRMSE_M[i] = rmse(np.abs(Z), np.abs(Zoutn))
    PhiRMSE_M[i] = rmse(np.angle(Z), np.angle(Zoutn))
    
    #Calculae R^2 of Measured Z
    ZreR2_M[i] = r2_score(Zre, Zoutn.real)
    ZimR2_M[i] = r2_score(Zim, Zoutn.imag)

    #Actual vs Measured 
    if tn <= 2: #length check to not show many plots
        #Nyq plot
        plt.figure(1 + i*countf)
        plt.plot(Zre, -Zim, label = 'Actual', marker = 'o', color  = 'purple', linestyle = 'None')
        plt.plot(Zoutn.real, -Zoutn.imag, label = 'Measured with tdelay (us) = ' + str(tdarrayu), marker = 'x', color = 'orange', linestyle = '-.')
        #plt.plot(ZoutECM.real, -ZoutECM.imag, label = 'ECM_EXT_v2', marker = 'x', color  = 'green', linestyle = '-')
        plt.plot(ZoutECM.real, -ZoutECM.imag, label = 'ECM_EXT_C_v1', marker = 'x', color  = 'green', linestyle = '-')
        addlabels(f, Zoutn.real, -Zoutn.imag)
        plt.title('Nyquist Plot')
        plt.xlabel('Zreal (\u03A9)')
        plt.ylabel('-Zim (\u03A9)')
        plt.grid(b=True, which='both', linestyle='-')
        plt.legend()
    
    if tn <= 2 and False: #turn off extra plots
        #Mag plot
        plt.figure(2 + i*countf)
        plt.title('Bode Magnitude Plot')
        plt.plot(f, np.abs(Z), label = 'Actual', marker = 'o', color  = 'purple', linestyle = '-')
        plt.plot(f, np.abs(Zoutn), label = 'Measured with tdelay (us) = ' + str(tdarrayu), marker = 'x', color = 'orange', linestyle = '-') 
        plt.xscale('log')
        plt.grid(b=True, which='both', color='#DDDDDD', linestyle='-')
        plt.xlabel('f (Hz)')
        plt.ylabel('|Z| (\u03A9)')
        plt.legend()

        #Phase plot
        plt.figure(3 + i*countf)
        plt.title('Bode Phase Plot')
        plt.plot(f, np.angle(Z), label = 'Actual', marker = 'o', color  = 'purple', linestyle = '-')
        plt.plot(f, np.angle(Zoutn), label = 'Measured with tdelay (us) = ' + str(tdarrayu), marker = 'x', color = 'orange', linestyle = '-')
        plt.xscale('log')
        plt.grid(b=True, which='both', color='#DDDDDD', linestyle='-')
        plt.xlabel('f (Hz)')
        plt.ylabel('\u03C6 (rad)')
        plt.legend()
        
        #Zre Resid plot
        plt.figure(4 + i*countf)
        plt.title('Zreal Residuals' + ', Measured with tdelay (us) = ' + str(tdarrayu))
        plt.plot(f, reResid_M, marker = 'o', color  = 'blue', linestyle = 'None')
        plt.xscale('log')
        plt.grid(b=True, which='both', color='#DDDDDD', linestyle='-')
        plt.xlabel('f (Hz)')
        plt.ylabel('Zreal (\u03A9)')
        plt.axhline(0, color='black')
        
        #Zim Resid plot
        plt.figure(5 + i*countf)
        plt.title('Zim Residuals' + ', Measured with tdelay (us) = ' + str(tdarrayu))
        plt.plot(f, -imResid_M, marker = 'o', color  = 'red', linestyle = 'None') #Use negative residuals for negative imaginary impedance ploted
        plt.xscale('log')
        plt.grid(b=True, which='both', color='#DDDDDD', linestyle='-')
        plt.xlabel('f (Hz)')
        plt.ylabel('-Zim (\u03A9)')
        plt.axhline(0, color='black')

    '''if tn <= 4: #length check to not show many plots
        plt.figure(1 + i*countf)
        plt.plot(Zre, -Zim, label = 'Actual Data', marker = 'o', color  = 'purple', linestyle = 'None')
        plt.plot(Zoutn.real, -Zoutn.imag, label = 'System Buffer with tdelay (s) = ' + str(tdarray[i]), marker = 'x', color = 'orange', linestyle = '-.')
        
        addlabels(f, Zoutn.real, -Zoutn.imag)
        
        plt.plot(ZoutECM.real, -ZoutECM.imag, label = 'ECM_EXT', color  = 'green', linestyle = '-')
        
        lmin = (diff(np.sign(diff(-Zoutn.imag)/diff(Zoutn.real))) > 0).nonzero()[0] + 1
        lmax = (diff(np.sign(diff(-Zoutn.imag)/diff(Zoutn.real))) < 0).nonzero()[0] + 1
        
        #plt.plot(Zoutn.real[lmin],-Zoutn.imag[lmin], 'kv', label = 'local min')
        #plt.plot(Zoutn.real[lmax],-Zoutn.imag[lmax], 'r^', label = 'local max')
        
        plt.title('Nyquist Plot')
        plt.xlabel('Zreal (\u03A9)')
        plt.ylabel('-Zim (\u03A9)')
        plt.grid(b=True, which='both', linestyle='-')
        plt.legend()
        
        plt.figure(2 + i*countf)
        plt.title('Bode Magnitude Plot')
        plt.plot(f, np.abs(Z), label = 'Actual Data', marker = 'o', color  = 'purple', linestyle = '-')
        plt.plot(f, np.abs(Zoutn), label = 'System Buffer with tdelay (s) = ' + str(tdarray[i]), marker = 'x', color = 'orange', linestyle = '-')
        #plt.plot(f, np.abs(ZoutECM), label = 'ECM_EXT', marker = 'D', color  = 'green', linestyle = '-')
        
        plt.xscale('log')
        plt.grid(b=True, which='both', color='#DDDDDD', linestyle='-')
        plt.xlabel('f (Hz)')
        plt.ylabel('|Z| (\u03A9)')
        plt.legend()
        
        plt.figure(3 + i*countf)
        plt.title('Bode Phase Plot')
        plt.plot(f, np.angle(Z), label = 'Actual Data', marker = 'o', color  = 'purple', linestyle = '-')
        plt.plot(f, np.angle(Zoutn), label = 'System Buffer with tdelay (s) = ' + str(tdarray[i]), marker = 'x', color = 'orange', linestyle = '-')
        #plt.plot(f, np.angle(ZoutECM), label = 'ECM_EXT', marker = 'D', color  = 'green', linestyle = '-')
        plt.xscale('log')
        plt.grid(b=True, which='both', color='#DDDDDD', linestyle='-')
        plt.xlabel('f (Hz)')
        plt.ylabel('\u03C6 (rad)')
        plt.legend()'''
    
if tn > 1:
    plt.figure(countf+1 + (tn-1)*countf)
    plt.title('ECM_EXT Parameter Error for Varying tdelay with 1st-order Randles Model')
    
    #print Zoutn Zout residuals %, chisquar, rmse
    
    tdscale = tdarray*10**6
    
    plt.plot(tdscale, R0ERR, label = 'R0')
    plt.plot(tdscale, R1ERR, label = 'R1')
    plt.plot(tdscale, C1ERR, label = 'C1')
    plt.xlabel('Time Delay (us)')
    #plt.xscale('log')
    #should be log scale... since its db or not can just consider it as is
    plt.ylabel('% Error')
    plt.legend()
    
    
#Print Table for RMSE, r^2, (chi-squared, l2 norm of (residual)
#Models that have worse predictions than this baseline will have a negative R^2
print('\n')
print('Actual vs Measured Impedance')

#large table with R^2 and RSS
#headers = ["Tdelay (us)", "Zre RMSE (\u03A9)", "Zre R\u00b2", 'Zre RSS (\u03A9\u00b2)', "Zim RMSE (\u03A9) ", "Zim R\u00b2", 'Zim RSS (\u03A9\u00b2)', '|Z| RMSE (\u03A9)', '\u03C6 RMSE (rad)']
#Tdata = np.array([tdarray*10**6, ZreRMSE_M, ZreR2_M, ZreRSS_M, ZimRMSE_M, ZimR2_M, ZimRSS_M, MagRMSE_M, PhiRMSE_M]).T

'''headers = ["Tdelay (us)", "Zre RMSE (\u03A9)", "Zre R\u00b2", 'Zre RSS (\u03A9\u00b2)', "Zim RMSE (\u03A9) ", "Zim R\u00b2", 'Zim RSS (\u03A9\u00b2)']
Tdata = np.array([tdarray*10**6, ZreRMSE_M, ZreR2_M, ZreRSS_M, ZimRMSE_M, ZimR2_M, ZimRSS_M,]).T'''

'''headers = ["Tdelay (us)", "Zre RMSE", "Zre R\u00b2", "Zim RMSE", "Zim R\u00b2"]
Tdata = np.array([tdarray*10**6, ZreRMSE_M, ZreR2_M, ZimRMSE_M, ZimR2_M]).T'''

#concise table only RMSEs for impedance generated from noise, most important
headers = ["Tdelay (us)", "Zre RMSE (\u03A9)", "Zim RMSE (\u03A9) ", '|Z| RMSE (\u03A9)', '\u03C6 RMSE (rad)']
Tdata = np.array([tdarray*10**6, ZreRMSE_M, ZimRMSE_M, MagRMSE_M, PhiRMSE_M]).T


# tabulate data
#table = tabulate(m, headers, tablefmt="fancy_grid") #weord but clear table
table = tabulate(Tdata, headers)

# output
print(table)