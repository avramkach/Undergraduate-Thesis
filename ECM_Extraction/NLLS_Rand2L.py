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
from NLLS_Rand2 import *

@author: avramkachura
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import least_squares


def NLLS_Rand2L(f, Zre, Zim):
    print('-----------' + 'NLLS_Rand2' + '-----------')
    
    #Filter out Zre<
    invalre = np.where(Zre < 0)[0]
    Zre = np.delete(Zre, invalre)
    Zim = np.delete(Zim, invalre)
    ft = np.delete(f, invalre) #dont delete from f since need to generate ECM with same size and same spectrum
    #use ft to analyze filtered out data and f for ZoutECM
    
    
    #NLLS Parameters, lm  = Levenberg-Marquardt algorithm
    #no bounds
    
    #L, R0, R1, C1, R2, C2 (-> R2, C2 SEI HF, low values), sig
    guess = (0.000000001, 1, 1, 2, 0.1, 0.01, 1)
    
    #set gradient tolerance very low
    gtol = 5e-16 
    
    #use ft for fitting filtered data
    ECM = least_squares(Z_resid, guess, args=(Zre, Zim, ft), gtol = gtol, method = 'lm')

    #Circuit Parameters
    Le = ECM.x[0]
    R0e = ECM.x[1]
    R1e = ECM.x[2]
    C1e = ECM.x[3]
    R2e = ECM.x[4]
    C2e = ECM.x[5]
    sig = ECM.x[6]
    
    #number of calls to the function:
    count = ECM.nfev
    
    print('L: %.10f' % Le)
    print('Rohm (Rinf, R0): %.7f' % R0e)
    print('Rct (R1): %.7f' % R1e)
    print('Cdl (C1): %.7f' % C1e)
    print('Rsei (R2): %.7f' % R2e)
    print('Csei (C2): %.7f' % C2e)
    print("Warburg Param Aw/\u03C3 = %.7f" % sig)
    print('Iterations: ', count)
    print(ECM.message)
    
    #sub in parameters with frequency given
    Zoutre, Zoutim = f_Randles2(ECM.x, f)
    ZoutECM = Zoutre + Zoutim*1j
    
    return R0e, R1e, C1e, ZoutECM

'''ecm real and imaginary function, 2RC circuit
inputs 
    -x: circuit parameter array: [R0, R1, C1]
    -f: frequency points of actual data

output: Real and Imaginary Impedance'''
def f_Randles2(ECM, f):
    L = ECM[0]
    R0 = ECM[1]
    R1 = ECM[2]
    C1 = ECM[3]
    R2 = ECM[4]
    C2 = ECM[5]
    sig = ECM[6]
    
    Zc1 = 1/(2*np.pi*f*C1*1j)
    Zc2 = 1/(2*np.pi*f*C2*1j)
    Zw = sig/(np.sqrt(2*np.pi*f))*(1-1j)
    
    Z = 2*np.pi*f*L*1j + R0 + 1/(1/(R1+Zw)+1/Zc1) + 1/(1/R2+1/Zc2)
    
    return Z.real, Z.imag


'''ecm residual function, error(x)
inputs 
    -x: circuit parameter array parameters: [R0, R1, C1]
    -f: frequency points of actual data
outputs
    -stack of residuals of weighted real and imaginary impedance'''
def Z_resid(ECM, Zre, Zim, f):
    Zre_sim, Zim_sim = f_Randles2(ECM, f)
    
    #weight vector (LKK)
    w = 1/(Zre**2 + Zim**2)
    
    Zre_error = w*(Zre - Zre_sim)**2
    Zim_error = w*(Zim - Zim_sim)**2
    
    #1d array
    return np.hstack((Zre_error, Zim_error))