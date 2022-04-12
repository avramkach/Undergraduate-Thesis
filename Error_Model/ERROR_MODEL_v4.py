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
import os
import h5py
from sklearn.metrics import r2_score #r^2 value
from tabulate import tabulate #printing table


import matplotlib as mpl
mpl.rcParams['figure.dpi'] = 300 #set figure to higher default dpi

#Import function in different directory needs to be folder path not file
sys.path.insert(1, '/Users/avramkachura/Desktop/Fourth Year/ESC499/ECM_Extraction/') 

from ECM_EXT_v2 import *

from ECM_EXT_v3 import *

from ECM_EXT_C_v1 import *

sys.path.insert(1, '/Users/avramkachura/Desktop/Fourth Year/ESC499/Python Error Model/')

from ERROR_MODEL_function_v4 import *

sys.path.insert(1, '/Users/avramkachura/Desktop/Fourth Year/ESC499/ECM_Extraction/NLLS') 

from NLLS_Rand1 import *

#Disabling and Enabling Print statements
# Disable
class HiddenPrints:
    def __enter__(self):
        self._original_stdout = sys.stdout
        sys.stdout = open(os.devnull, 'w')

    def __exit__(self, exc_type, exc_val, exc_tb):
        sys.stdout.close()
        sys.stdout = self._original_stdout

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

#circuit model, set param defaults to 0, output mohm
def Zgen(f, R0 = 0, R1 = 0, C1 = 0, sig = 0, L = 0): 
    '''L = 0
    sig = 0
    R0 = 0.25
    R1 = 1
    C1 = 0.25'''
    
    #ECM EXT paper extracted from simulated    
    '''L = 0.17549*10**(-6)
    R0 = 0.5509587
    R1 = 0.1193365
    C1 = 1.14637
    sig = 0.0346208'''
    
    #Cell5_0 from EIS Analyzer fit
    R0 = 0.26424
    C1 = 0.42183
    R1 = 0.21616
    sig = 0.11111
    
    '''R0 = 0.26424
    C1 = 0.42183/np.pi
    R1 = 0.21616
    sig = 0.11111'''
    
    
    '''R0 = 0.26424*10**(-3) #mOhms
    C1 = 0.42183*10**(3) #KF
    R1 = 0.21616*10**(-3)
    sig = 0.11111*10**(-3)'''

    if C1 != 0:
        Zc = 1/(2*pi*f*C1*1j)
        
    ZL = 1j*2*pi*f*L;
    
    if sig != 0:
        Zw = sig/(np.sqrt(2*pi*f))*(1-1j) #Semifinite warburg
    
    #R + R//C circuit
    #Z = R0 + 1/(1/R1+1/Zc)
    
    #ECM EXT paper 2nd order randles model   
    #Z = R0 + 1/(1/(R1+Zw)+1/Zc) + ZL + 1/(1/0.0135778+2*pi*f*0.0344249)
    
    
    
    
    #1st order randles model Cell5_0 EIS Analyzer fit
    Z = R0 + 1/(1/(R1+Zw)+1/Zc)
    
    '''R0 = 0.1
    n = len(f)
    Z = np.ones(n)'''
    #Z = Zw
    
    return R0, R1, C1, Z

#add frequency labels for noise plotting with frequency cutoff
#https://queirozf.com/entries/add-labels-and-text-to-matplotlib-plots-annotation-examples
def addlabels(f, Zre, Zim):
    fcutoff = 0
    inval = np.where(f > fcutoff)
    f = np.delete(f, inval[0])
    ys = np.delete(-Zim, inval[0])
    xs = np.delete(Zre, inval[0])
    n = len(f)
    for i in range(0, n):
        
        label = "{:.2f}Hz".format(f[i])
        #label = "{:.2f}".format(f[i])
        if i % 2 == 0:
            plt.annotate(label, # this is the text
                         (xs[i],-ys[i]), # these are the coordinates to position the label
                         textcoords="offset points", # how to position the text
                         xytext=(0,10), # distance from text to points (x,y)
                         ha='center') # horizontal alignment can be left, right or center
        else:
            plt.annotate(label, # this is the text
                         (xs[i],-ys[i]), # these are the coordinates to position the label
                         textcoords="offset points", # how to position the text
                         xytext=(0,-15), # distance from text to points (x,y)
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
    
#f = np.logspace(4,-2,30)

#f = np.array([1000,20,15,10,1])

#Evaluate RC elements for now, turn to mOhm
R0, R1, C1, Z = Zgen(f)
Z = Z
Zre = Z.real
Zim = Z.imag


#Define Time Delays
#0.00001 = 10us = 1e-5
#tdarray = np.array([0.000001, 0.00001, 0.0001, 0.0005, 0.001, 0.005]) # tests shown zhe, change to us and print at top of plot
tdarray = np.array([0.00000]) #10us

#tdarray = np.array([0.00005]) #50us

#tdarray = np.array([0.001]) #at this frequency no inductive element robust to this change

#tdarray = np.array([0.01]) #10ms really breaks end at 1e-2
#tdarray = np.array([0.00000001]) #0.01us fimax = 1952, T = 0.0051229 /100 = 5us of period
#tdarray = np.array([0.0000001]) #0.1us start at 1e-7
#tdarray = np.linspace(1e-8, 1e-2, 100) #0.01us to 10ms
#tdarray = np.append([0], np.linspace(1e-8, 1e-2, 100))
#tdarray = np.append([0], np.logspace(-8,-2,100))
#breaks down at 15us and 50us delay looking at the imaginary impedance statistics
tn = len(tdarray)

#Disable SNR by setting to inf same length as tn
#SNRarray = np.ones(tn)*np.inf

#SNRarray = np.array([10, 40])
#SNRarray = np.arange(1,101)
#SNRarray = np.linspace(1,100,10)
SNRarray = np.array([1])
sn = len(SNRarray)

'''#Param Error Arrays % = |(estimate - actual)/actual| * 100%
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
ZimR2_M = np.zeros(tn)'''


#Measured statistics
#Residual sum of squares RSS = sum(residual^2)
ZreRSS_M = np.zeros(sn)
ZimRSS_M = np.zeros(sn)

#RMSE Arrays M
ZreRMSE_M = np.zeros(sn)
ZimRMSE_M = np.zeros(sn)
MagRMSE_M = np.zeros(sn)
PhiRMSE_M = np.zeros(sn)

#R^2 arrays of coefficient of determination, degree of any linear correlation between y_true and y_pred, r2_score(y_true, y_pred)
ZreR2_M = np.zeros(sn)
ZimR2_M = np.zeros(sn)

#stdev arrays
ZreSTDEV_M = np.zeros(sn)
ZimSTDEV_M = np.zeros(sn)


#iterations for error model random noise
iters = 1

#initialize impedance storage measured test matrices
Zre_M = np.zeros((len(f), iters))
Zim_M = np.zeros((len(f), iters)) #each row is a impedance measurement for a specific frequency
#i.e. Zim_M[0] = Zim1,...Zim100 for 1952Hz
#set Zim_M[:,i] = Zoutn.image(populates each column)

#ECM EXT Fit statistics
#ECM_EXT_c_v1 stats
#initialize fit impedance storage matrices
Zoutre_c_v1 = np.zeros((len(f), iters)) #generated impedances for each iteration
Zoutim_c_v1 = np.zeros((len(f), iters))

#Fit parameters
R0out_c_v1 = np.zeros(iters)
R1out_c_v1 = np.zeros(iters)
C1out_c_v1 = np.zeros(iters)

#Fit parameter relative RSD = 100*stdev/mean
R0out_c_v1_RSD = np.zeros(sn)
R1out_c_v1_RSD = np.zeros(sn)
C1out_c_v1_RSD = np.zeros(sn)

#Param Error Arrays % = |(estimate - actual)/actual| * 100%
R0ERR_c_v1 = np.zeros(sn)
R1ERR_c_v1 = np.zeros(sn)
C1ERR_c_v1 = np.zeros(sn)

#ECM EXT impedance rmse
ZreRMSE_out_c_v1 = np.zeros(sn)
ZimRMSE_out_c_v1 = np.zeros(sn)

#ECM EXT average time over each SNR
t_mean_c_v1 = np.zeros(sn)
t_c_v1 = np.zeros(iters)

#NLLS stats
#initialize fit impedance storage matrices
Zoutre_NLLS = np.zeros((len(f), iters)) #generated impedances for each iteration
Zoutim_NLLS = np.zeros((len(f), iters))

#Fit parameters
R0out_NLLS = np.zeros(iters)
R1out_NLLS = np.zeros(iters)
C1out_NLLS = np.zeros(iters)

#Fit parameter relative RSD = 100*stdev/mean
R0out_NLLS_RSD = np.zeros(sn)
R1out_NLLS_RSD = np.zeros(sn)
C1out_NLLS_RSD = np.zeros(sn)

#Param Error Arrays % = |(estimate - actual)/actual| * 100%
R0ERR_NLLS = np.zeros(sn)
R1ERR_NLLS = np.zeros(sn)
C1ERR_NLLS = np.zeros(sn)

#NLLS impedance rmse
ZreRMSE_out_NLLS = np.zeros(sn)
ZimRMSE_out_NLLS = np.zeros(sn)

#NLLS average time over each SNR
t_mean_NLLS = np.zeros(sn)
t_NLLS = np.zeros(iters)
#iter_NLLS = np.zeros(tn*sn) # not used

countf = 3 #figure counter in loop + 2 for residual

for i in range(0,sn):

    print('-----' + 'Testing Time Delay (us) = ' + str(tdarray[0]*10**6) + '-----')
    print('-----' + 'Testing Noise (dB) = ' + str(SNRarray[i]) + '-----')
    
    for j in range(0, iters):
        #us conversion for printing
        tdarrayu = tdarray[0]*10**6
        
        #Noise Model
        Zoutn = ERROR_MODEL_v4(f, Zre, Zim, SNRarray[i], tdarray[0]) 
        
        #ECM extract estimates for R0e
        #HideenPrints block ECM prints
        #with HiddenPrints():
        t0_c_v1 = time.time()
        R0e_c_v1, R1e_c_v1, C1e_c_v1, ZoutECM_c_v1 = ECM_EXT_C_v1(f, Zoutn.real, Zoutn.imag)
        t1_c_v1 = time.time()
        t_c_v1[j] = t1_c_v1 - t0_c_v1
        print(t_c_v1)
        #R0e, R1e, C1e, ZoutECM = ECM_EXT_v3(f, Zoutn.real, Zoutn.imag)
        
        #1st order randles
        t0_NLLS = time.time()
        R0e_NLLS, R1e_NLLS, C1e_NLLS, ZoutECM_NLLS = NLLS_Rand1(f, Zoutn.real, Zoutn.imag)
        t1_NLLS = time.time()
        t_NLLS[j] = t1_NLLS - t0_NLLS
        print(t_NLLS)
        
        #Populate and store ECM params
        #ECM_EXT_c_v1
        R0out_c_v1[j] = R0e_c_v1
        R1out_c_v1[j] = R1e_c_v1
        C1out_c_v1[j] = C1e_c_v1
        #NLLS
        R0out_NLLS[j] = R0e_NLLS
        R1out_NLLS[j] = R1e_NLLS
        C1out_NLLS[j] = C1e_NLLS
            
        #Populate impedance storage measured test matrices
        Zre_M[:,j] = Zoutn.real
        Zim_M[:,j] = Zoutn.imag
        
        #Populate impedance fit matrices
        #ECM_EXT_c_v1
        Zoutre_c_v1[:,j] = ZoutECM_c_v1.real
        Zoutim_c_v1[:,j] = ZoutECM_c_v1.imag
        #NLLS
        Zoutre_NLLS[:,j] = ZoutECM_NLLS.real
        Zoutim_NLLS[:,j] = ZoutECM_NLLS.imag
            
        #Actual vs Measured 
        if sn <= 2 and iters <= 5: #length check to not show many plots
            #Nyq plot
            plt.figure(1 + i*countf + j)
            plt.plot(Zre, -Zim, label = 'Actual', marker = 'o', color  = 'purple', linestyle = 'None')
            plt.plot(Zoutn.real, -Zoutn.imag, label = 'Measured with SNR (dB) = ' + str(SNRarray[i]), marker = 'x', color = 'orange', linestyle = '-.')
            
            #for powerpoint
            #plt.plot(Zre, -Zim, marker = 'o', color  = 'purple', linestyle = 'None')
            #plt.plot(Zre, -Zim, label = 'Manually Fitted ECM', marker = 'o', markersize = 2, color  = 'green', linestyle = '-')
            #plt.plot(Zoutn.real, -Zoutn.imag, label = 'Measured with SNR (dB) = ' + str(SNRarray[i]), marker = 'x', color = 'orange', linestyle = '-.')
            
            #addlabels(f, Zre, -Zim)
            #plt.plot(ZoutECM_c_v1.real, -ZoutECM_c_v1.imag, label = 'ECM_EXT_c_v1', marker = 'x', color  = 'green', linestyle = '-')
            #plt.plot(ZoutECM.real, -ZoutECM.imag, label = 'ECM_EXT_v3', marker = 'x', color  = 'green', linestyle = '-')
            #plt.plot(ZoutECM_C_v1.real, -ZoutECM_C_v1.imag, label = 'ECM_EXT_C_v1', marker = '.', color  = 'black', linestyle = '-')
            plt.plot(ZoutECM_c_v1.real, -ZoutECM_c_v1.imag, label = 'ECM_EXT_c_v1', marker = '.', color  = 'black', linestyle = '-')
            plt.plot(ZoutECM_NLLS.real, -ZoutECM_NLLS.imag, label = 'NLLS', marker = '.', color  = 'grey', linestyle = '-')
            
            if False: #len(f) < 30:
                addlabels(f, Zoutn.real, -Zoutn.imag)
            plt.title('Nyquist Plot')
            plt.xlabel('Zreal (m\u03A9)')
            plt.ylabel('-Zim (m\u03A9)')
            plt.grid(b=True, which='both', linestyle='-')
            plt.legend()
            
            '''plt.figure(1 + i*countf + j+30)
            plt.plot(Zre, -Zim, label = 'Actual', marker = 'o', color  = 'purple', linestyle = 'None')
            plt.plot(Zoutn.real, -Zoutn.imag, label = 'Measured with SNR (dB) = ' + str(SNRarray[i]), marker = 'x', color = 'orange', linestyle = '-.')
            
            plt.plot(ZoutECM_v3.real, -ZoutECM_v3.imag, label = 'ECM_EXT_v3', marker = 'x', color  = 'green', linestyle = '-')
            #plt.plot(ZoutECM_C_v1.real, -ZoutECM_C_v1.imag, label = 'ECM_EXT_C_v1', marker = '.', color  = 'black', linestyle = '-')
            
            addlabels(f, Zoutn.real, -Zoutn.imag)
            plt.title('Nyquist Plot')
            plt.xlabel('Zreal (m\u03A9)')
            plt.ylabel('-Zim (m\u03A9)')
            plt.grid(b=True, which='both', linestyle='-')
            plt.legend()'''
    
    #get average times
    t_mean_c_v1[i] = np.mean(t_c_v1)
    t_mean_NLLS[i] = np.mean(t_NLLS)
    
    #Get impedance mean and variance from rows
    #Measured
    Zre_Mmean = np.mean(Zre_M, axis = 1)
    Zim_Mmean = np.mean(Zim_M, axis = 1)
    Zre_Mstdev = np.std(Zre_M, axis = 1)
    Zim_Mstdev = np.std(Zim_M, axis = 1)
    
    #ECM Fit 
    #ECM_EXT_v1_c
    Zoutre_c_v1_mean = np.mean(Zoutre_c_v1, axis = 1)
    Zoutim_c_v1_mean = np.mean(Zoutim_c_v1, axis = 1)
    Zoutre_c_v1_stdev = np.std(Zoutre_c_v1, axis = 1)
    Zoutim_c_v1_stdev = np.std(Zoutim_c_v1, axis = 1)
    R0out_c_v1_mean = np.mean(R0out_c_v1)
    R1out_c_v1_mean = np.mean(R1out_c_v1)
    C1out_c_v1_mean = np.mean(C1out_c_v1)
    #NLLS
    Zoutre_NLLS_mean = np.mean(Zoutre_NLLS, axis = 1)
    Zoutim_NLLS_mean = np.mean(Zoutim_NLLS, axis = 1)
    Zoutre_NLLS_stdev = np.std(Zoutre_NLLS, axis = 1)
    Zoutim_NLLS_stdev = np.std(Zoutim_NLLS, axis = 1)
    R0out_NLLS_mean = np.mean(R0out_NLLS)
    R1out_NLLS_mean = np.mean(R1out_NLLS)
    C1out_NLLS_mean = np.mean(C1out_NLLS)
    
    
    '''R0ERR[i] = np.abs((R0e-R0)/R0)*100
    if R1 != 0 and C1 != 0: #R1 = 0 iff C1 = 0 where Zim peak isnt at 0
        R1ERR[i] = np.abs((R1e-R1)/R1)*100
        C1ERR[i] = np.abs((C1e-C1)/C1)*100'''
    
    R0ERR_c_v1[i] = np.abs((R0out_c_v1_mean-R0)/R0)*100
    R0out_c_v1_RSD[i] = np.std(R0out_c_v1)/R0out_c_v1_mean*100
    
    R0ERR_NLLS[i] = np.abs((R0out_NLLS_mean-R0)/R0)*100
    R0out_NLLS_RSD[i] = np.std(R0out_NLLS)/R0out_NLLS_mean*100
    
    if R1 != 0 or C1 != 0:
        R1ERR_c_v1[i] = np.abs((R1out_c_v1_mean-R1)/R1)*100
        C1ERR_c_v1[i] = np.abs((C1out_c_v1_mean-C1)/C1)*100
        R1out_c_v1_RSD[i] = np.std(R1out_c_v1)/R1out_c_v1_mean*100
        C1out_c_v1_RSD[i] = np.std(C1out_c_v1)/C1out_c_v1_mean*100
        
        R1ERR_NLLS[i] = np.abs((R1out_NLLS_mean-R1)/R1)*100
        C1ERR_NLLS[i] = np.abs((C1out_NLLS_mean-C1)/C1)*100
        R1out_NLLS_RSD[i] = np.std(R1out_NLLS)/R1out_NLLS_mean*100
        C1out_NLLS_RSD[i] = np.std(C1out_NLLS)/C1out_NLLS_mean*100
    
    ZreRMSE_out_c_v1[i] = rmse(Zre, Zoutre_c_v1_mean)
    ZimRMSE_out_c_v1[i] = rmse(Zim, Zoutim_c_v1_mean)
    
    ZreRMSE_out_NLLS[i] = rmse(Zre, Zoutre_NLLS_mean)
    ZimRMSE_out_NLLS[i] = rmse(Zim, Zoutim_NLLS_mean)
    
    #Calculate and store measured stats
    #Calculate Residuals of Measured Z (Residual=actual−predicted = Actual - Measured)
    reResid_M = Zre - Zre_Mmean
    imResid_M = Zim - Zim_Mmean
    
    #Calculated RSS of Measured Z
    ZreRSS_M[i] = RSS(reResid_M)
    ZimRSS_M[i] = RSS(imResid_M)
    
    #Calculae RMSE of Measured Z
    ZreRMSE_M[i] = rmse(Zre, Zre_Mmean)
    ZimRMSE_M[i] = rmse(Zim, Zim_Mmean)
    MagRMSE_M[i] = rmse(np.abs(Z), np.abs(Zre_Mmean + Zim_Mmean*1j))
    PhiRMSE_M[i] = rmse(np.angle(Z), np.angle(Zre_Mmean + Zim_Mmean*1j))
    
    #Calculae R^2 of Measured Z
    ZreR2_M[i] = r2_score(Zre, Zre_Mmean)
    ZimR2_M[i] = r2_score(Zim, Zim_Mmean)
    
    #Calculate stdev of Measured Z by summing up variances and sqrt
    ZreSTDEV_M[i] = np.sum(Zre_Mstdev**2)**0.5
    ZimSTDEV_M[i] = np.sum(Zim_Mstdev**2)**0.5
    
    #Main nyquist plot for 
    if sn <= 4: #length check to not show many plots
        #Nyq plot
        plt.figure(1 + i*countf + iters)
        plt.plot(Zre, -Zim, label = 'Actual', marker = 'o', color  = 'purple', linestyle = 'None')
        #plt.plot(Zre_Mmean, -Zim_Mmean, label = 'Measured Mean with SNR (dB) = ' + str(SNRarray[i]), marker = 'x', color = 'orange', linestyle = '-.')
        #plt.plot(ZoutECM.real, -ZoutECM.imag, label = 'ECM_EXT_C_v1', marker = 'x', color  = 'green', linestyle = '-')
        
        plt.errorbar(Zre_Mmean, -Zim_Mmean, yerr = Zim_Mstdev, xerr = Zre_Mstdev, label = 'Measured Mean with SNR (dB) = ' + str(SNRarray[i]), fmt = 'o', markersize = 5, color = 'orange', capsize = 2)
        
        plt.errorbar(Zoutre_c_v1_mean, -Zoutim_c_v1_mean, yerr = Zoutre_c_v1_stdev, xerr = Zoutim_c_v1_stdev, label = 'ECM_EXT_c_v1 Mean with SNR (dB) = ' + str(SNRarray[i]), fmt = '.', markersize = 5, color = 'black', capsize = 2)
        plt.errorbar(Zoutre_NLLS_mean, -Zoutim_NLLS_mean, yerr = Zoutre_NLLS_stdev, xerr = Zoutim_NLLS_stdev, label = 'NLLS Mean with SNR (dB) = ' + str(SNRarray[i]), fmt = '.', markersize = 5, color = 'grey', capsize = 2)
        
        if len(f) < 30:
            addlabels(f, Zre_Mmean, -Zim_Mmean)
        plt.title('Mean Nyquist Plot')
        plt.xlabel('Zreal (m\u03A9)')
        plt.ylabel('-Zim (m\u03A9)')
        plt.grid(b=True, which='both', linestyle='-')
        plt.legend()
    
    #Bode plots, mohm
    if False: #sn <= 2 and iters <= 5:
        #Mag plot
        plt.figure(2 + i*countf + iters)
        plt.title('Bode Magnitude Plot')
        plt.plot(f, np.abs(Z), label = 'Actual', marker = 'o', color  = 'purple', linestyle = '-')
        plt.plot(f, np.abs(Zoutn), label = 'Measured with SNR (dB) = ' + str(SNRarray[i]), marker = 'x', color = 'orange', linestyle = '-') 
        plt.xscale('log')
        plt.grid(b=True, which='both', color='#DDDDDD', linestyle='-')
        plt.xlabel('f (Hz)')
        plt.ylabel('|Z| (m\u03A9)')
        plt.legend()

        #Phase plot
        plt.figure(3 + i*countf + iters)
        plt.title('Bode Phase Plot')
        plt.plot(f, np.angle(Z), label = 'Actual', marker = 'o', color  = 'purple', linestyle = '-')
        plt.plot(f, np.angle(Zoutn), label = 'Measured with SNR (dB) = ' + str(SNRarray[i]), marker = 'x', color = 'orange', linestyle = '-')
        plt.xscale('log')
        plt.grid(b=True, which='both', color='#DDDDDD', linestyle='-')
        plt.xlabel('f (Hz)')
        plt.ylabel('\u03C6 (rad)')
        plt.legend()
        
        #Zre Resid plot
        plt.figure(4 + i*countf + iters)
        plt.title('Zreal Residuals' + ', Measured with tdelay (us) = ' + str(tdarrayu))
        plt.plot(f, reResid_M, marker = 'o', color  = 'blue', linestyle = 'None')
        plt.xscale('log')
        plt.grid(b=True, which='both', color='#DDDDDD', linestyle='-')
        plt.xlabel('f (Hz)')
        plt.ylabel('Zreal (m\u03A9)')
        plt.axhline(0, color='black')
        
        #Zim Resid plot
        plt.figure(5 + i*countf + iters)
        plt.title('Zim Residuals' + ', Measured with tdelay (us) = ' + str(tdarrayu))
        plt.plot(f, -imResid_M, marker = 'o', color  = 'red', linestyle = 'None') #Use negative residuals for negative imaginary impedance ploted
        plt.xscale('log')
        plt.grid(b=True, which='both', color='#DDDDDD', linestyle='-')
        plt.xlabel('f (Hz)')
        plt.ylabel('-Zim (m\u03A9)')
        plt.axhline(0, color='black')

#
'''if tn > 1:
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
    plt.legend()'''

#SNR plotting, C_v1
if sn > 1:
    plt.figure(countf+1 + (sn-1)*countf*iters)
    plt.title('ECM_EXT_c_v1 Parameter Error for Varying SNR with 1st-order Randles Model')
    plt.errorbar(SNRarray, R0ERR_c_v1, yerr = R0out_c_v1_RSD, label = 'R0', linestyle = 'solid', capsize = 2)
    plt.errorbar(SNRarray, R1ERR_c_v1, yerr = R1out_c_v1_RSD, label = 'R1', linestyle = 'solid', capsize = 2)
    plt.errorbar(SNRarray, C1ERR_c_v1, yerr = C1out_c_v1_RSD, label = 'C1', linestyle = 'solid', capsize = 2)
    plt.xlabel('SNR (dB)')
    #plt.xscale('log')
    #should be log scale... since its db or not can just consider it as is
    plt.ylabel('% Error')
    plt.legend()
    
    plt.figure(countf+1 + (sn-1)*countf*iters+1)
    plt.title('ECM_EXT_c_v1 Computation Time for Varying SNR with 1st-order Randles Model')
    plt.plot(SNRarray, t_mean_c_v1, label = 'Mean')
    plt.xlabel('SNR (dB)')
    plt.ylabel('Time (s)')
    plt.legend()
    
#SNR plotting, NLLS
if sn > 1:
    plt.figure(countf+1 + (sn-1)*countf*iters+2)
    plt.title('NLLS Parameter Error for Varying SNR with 1st-order Randles Model')
    plt.errorbar(SNRarray, R0ERR_NLLS, yerr = R0out_NLLS_RSD, label = 'R0', linestyle = 'solid', capsize = 2)
    plt.errorbar(SNRarray, R1ERR_NLLS, yerr = R1out_NLLS_RSD, label = 'R1', linestyle = 'solid', capsize = 2)
    plt.errorbar(SNRarray, C1ERR_NLLS, yerr = C1out_NLLS_RSD, label = 'C1', linestyle = 'solid', capsize = 2)
    plt.xlabel('SNR (dB)')
    plt.ylabel('% Error')
    plt.legend()
    
    plt.figure(countf+1 + (sn-1)*countf*iters+3)
    plt.title('NLLS Computation Time for Varying SNR with 1st-order Randles Model')
    plt.plot(SNRarray, t_mean_NLLS, label = 'Mean')
    plt.xlabel('SNR (dB)')
    plt.ylabel('Time (s)')
    plt.legend()

#_v3
'''if sn > 1:
    plt.figure(countf+1 + (sn-1)*countf*iters + 1)
    plt.title('ECM_EXT_v3 Parameter Error for Varying SNR with 1st-order Randles Model')
    
    plt.errorbar(SNRarray, R0ERR, yerr = R0out_RSD, label = 'R0', linestyle = 'solid', capsize = 2)
    plt.errorbar(SNRarray, R1ERR, yerr = R1out_RSD, label = 'R1', linestyle = 'solid', capsize = 2)
    plt.errorbar(SNRarray, C1ERR, yerr = C1out_RSD, label = 'C1', linestyle = 'solid', capsize = 2)
    
    plt.xlabel('SNR (dB)')
    #plt.xscale('log')
    #should be log scale... since its db or not can just consider it as is
    plt.ylabel('% Error')
    plt.legend()'''
    
    
#Print Table for RMSE, r^2, (chi-squared, l2 norm of (residual)
#Models that have worse predictions than this baseline will have a negative R^2
print('\n')
print('Actual vs Mean Measured Impedance')

#large table with R^2 and RSS
#headers = ["Tdelay (us)", "Zre RMSE (\u03A9)", "Zre R\u00b2", 'Zre RSS (\u03A9\u00b2)', "Zim RMSE (\u03A9) ", "Zim R\u00b2", 'Zim RSS (\u03A9\u00b2)', '|Z| RMSE (\u03A9)', '\u03C6 RMSE (rad)']
#Tdata = np.array([tdarray*10**6, ZreRMSE_M, ZreR2_M, ZreRSS_M, ZimRMSE_M, ZimR2_M, ZimRSS_M, MagRMSE_M, PhiRMSE_M]).T

'''headers = ["Tdelay (us)", "Zre RMSE (\u03A9)", "Zre R\u00b2", 'Zre RSS (\u03A9\u00b2)', "Zim RMSE (\u03A9) ", "Zim R\u00b2", 'Zim RSS (\u03A9\u00b2)']
Tdata = np.array([tdarray*10**6, ZreRMSE_M, ZreR2_M, ZreRSS_M, ZimRMSE_M, ZimR2_M, ZimRSS_M,]).T'''

'''headers = ["Tdelay (us)", "Zre RMSE", "Zre R\u00b2", "Zim RMSE", "Zim R\u00b2"]
Tdata = np.array([tdarray*10**6, ZreRMSE_M, ZreR2_M, ZimRMSE_M, ZimR2_M]).T'''

'''#concise table only RMSEs for impedance generated from noise, most important
headers = ["Tdelay (us)", "Zre RMSE (\u03A9)", "Zim RMSE (\u03A9) ", '|Z| RMSE (\u03A9)', '\u03C6 RMSE (rad)']
Tdata = np.array([tdarray*10**6, ZreRMSE_M, ZimRMSE_M, MagRMSE_M, PhiRMSE_M]).T'''

#Measured
#SNR concise table only RMSEs for impedance generated from noise, most important
headers = ["SNR (dB)", "Zre RMSE (m\u03A9)", "Zre \u03C3 (m\u03A9)", "Zim RMSE (m\u03A9) ", "Zim \u03C3 (m\u03A9)", '|Z| RMSE (m\u03A9)', '\u03C6 RMSE (rad)']
Tdata = np.array([SNRarray, ZreRMSE_M, ZreSTDEV_M, ZimRMSE_M, ZimSTDEV_M, MagRMSE_M, PhiRMSE_M]).T

# tabulate data
#table = tabulate(m, headers, tablefmt="fancy_grid") #weord but clear table
table = tabulate(Tdata, headers)

# output
print(table)



#print('Actual vs ECM_EXT_v3')

print('\n')
print('Actual vs ECM_EXT_c_v1')
headers = ["SNR (dB)", "Zre RMSE (m\u03A9)", "Zim RMSE (m\u03A9)", "R0 Mean Error (%)", "R1 Mean Error (%)", "C1 Mean Error (%)"]
Tdata = np.array([SNRarray, ZreRMSE_out_c_v1, ZimRMSE_out_c_v1, R0ERR_c_v1, R1ERR_c_v1, C1ERR_c_v1]).T
# tabulate data
table = tabulate(Tdata, headers)
# output
print(table)

print('\n')
print('Actual vs NLLS')
headers = ["SNR (dB)", "Zre RMSE (m\u03A9)", "Zim RMSE (m\u03A9)", "R0 Mean Error (%)", "R1 Mean Error (%)", "C1 Mean Error (%)"]
Tdata = np.array([SNRarray, ZreRMSE_out_NLLS, ZimRMSE_out_NLLS, R0ERR_NLLS, R1ERR_NLLS, C1ERR_NLLS]).T
# tabulate data
table = tabulate(Tdata, headers)
# output
print(table)