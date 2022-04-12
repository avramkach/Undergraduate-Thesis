#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr  4 23:42:45 2022

NLLS for 2RC circuits

Inputs: f (w = 2pif), Zre, Zim (preprocessed)

Function Inputs: (ordered with corresponding impedances)
    -f (Hz)
    -Zre (+) (same units as Zim, ohm)
    -Zim (+) (no negative, raw data)

Outputs: R0, R1, C1, ZoutECM (same units inputed as impedance data)

Assumes order: Takes in highest to lowest frequencies with corresponding impedances

Filters: 
    -All f points higher than lowest f Zim>0 
    -Zre<0

import sys
import os
sys.path.insert(1, '/Users/avramkachura/Desktop/Fourth Year/ESC499/ECM_Extraction/NLLS') 
from NLLS_Rand1 import *

@author: avramkachura
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import least_squares


def NLLS_Rand1(f, Zre, Zim):
    print('-----------' + 'NLLS_Rand1' + '-----------')
    
    #Filter out Zre<
    invalre = np.where(Zre < 0)[0]
    Zre = np.delete(Zre, invalre)
    Zim = np.delete(Zim, invalre)
    ft = np.delete(f, invalre) #dont delete from f since need to generate ECM with same size and same spectrum
    #use ft to analyze filtered out data and f for ZoutECM
    
    #Removes 0<Zim based on last frequency
    if np.size(np.where(Zim > 0)[0]) != 0: #if no -Zim>0 measurements  
        filim = np.max(np.where(Zim > 0)[0]) #take last element 
        Zre = np.delete(Zre, np.arange(0,filim+1))
        Zim = np.delete(Zim, np.arange(0,filim+1))
        ft = np.delete(ft, np.arange(0,filim+1)) 
    
    #NLLS Parameters, lm  = Levenberg-Marquardt algorithm
    #no bounds
    
    #R0, R1, C1, sig
    guess = (1, 1, 1, 1)
    
    #set gradient tolerance very low
    gtol = 5e-16 
    
    #use ft for fitting filtered data
    ECM = least_squares(Z_resid, guess, args=(Zre, Zim, ft), gtol = gtol, method = 'lm')

    #Circuit Parameters
    R0e = ECM.x[0]
    R1e = ECM.x[1]
    C1e = ECM.x[2]
    sig = ECM.x[3]
    
    #number of calls to the function:
    count = ECM.nfev
    
    print('Rohm (Rinf, R0): %.5f' % R0e)
    print('Rct (R1): %.5f' % R1e)
    print('Cdl (C1): %.5f' % C1e)
    print("Warburg Param Aw/\u03C3 = %.7f" % sig)
    print('Iterations: ', count)
    print(ECM.message)
    
    #sub in parameters with frequency given
    Zoutre, Zoutim = f_Randles1(ECM.x, f)
    ZoutECM = Zoutre + Zoutim*1j
    
    return R0e, R1e, C1e, ZoutECM

'''ecm real and imaginary function, randles circuit
inputs 
    -x: circuit parameter array: [R0, R1, C1]
    -f: frequency points of actual data

output: Real and Imaginary Impedance'''
def f_Randles1(ECM, f):
    R0 = ECM[0]
    R1 = ECM[1]
    C1 = ECM[2]
    sig = ECM[3]
    
    Zc1 = 1/(2*np.pi*f*C1*1j)
    
    Zw = sig/(np.sqrt(2*np.pi*f))*(1-1j)
    
    Z = R0 + 1/(1/(R1 + Zw) + 1/Zc1)
    
    return Z.real, Z.imag


'''ecm residual function, error(x)
inputs 
    -x: circuit parameter array parameters: [R0, R1, C1]
    -f: frequency points of actual data
outputs
    -stack of residuals of weighted real and imaginary impedance'''
def Z_resid(ECM, Zre, Zim, f):
    Zre_sim, Zim_sim = f_Randles1(ECM, f)
    
    #weight vector (LKK)
    w = 1/(Zre**2 + Zim**2)
    
    Zre_error = w*(Zre - Zre_sim)**2
    Zim_error = w*(Zim - Zim_sim)**2
    
    #1d array
    return np.hstack((Zre_error, Zim_error))