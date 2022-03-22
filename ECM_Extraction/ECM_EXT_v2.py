#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec 21 13:03:36 2021

ECM extraction algo 

Inputs: f (w = 2pif), Zre, Zim (preprocessed)

Function Inputs: (ordered with corresponding impedances)
    -f (Hz)
    -Zre (+) (same units as Zim, ohm)
    -Zim (+) (no negative, raw data)

Outputs: R0, R1, C1, ZoutECM (same units inputed as impedance data)

Assumes order: Takes in highest to lowest frequencies with corresponding impedances

Filters: 
    -Rounds impedance data to 15 digits, fine for now
    -Zre<0

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


def ECM_EXT_v2(f, Zre, Zim):
    #round to 15 decimals 32 bit -> 10 decimals (digits)
    Zim = np.round(Zim, 15)
    Zre = np.round(Zre, 15)
    
    print('-----------' + 'ECM_EXT_v2' + '-----------')
    
    
    '''print(Zre)
    print(Zim)'''
     
    #Output Impedance from ECM
    ZoutECM = np.zeros(len(f), dtype = "complex_")
    
    #Filter out Zre<
    invalre = np.where(Zre < 0)[0]
    Zre = np.delete(Zre, invalre)
    Zim = np.delete(Zim, invalre)
    ft = np.delete(f, invalre) #dont delete from f since need to generate ECM with same size and same spectrum
    w = 2*np.pi*ft #w is used for ECM EXT
    
    #Find Rohm = Rinf
    Rohm = Rinf(Zre, Zim)
    print('Rohm (Rinf, R0): %.5f' % Rohm)
    
    #Add Ohmic Resistance
    ZoutECM += Rohm
    
    #L inductive element
    L = 0
    XL = np.abs(np.min(-Zim)) 
    if XL >= 0 and XL == Zim[np.argmax(w)]: #robust to frequencies that arent the maximum frequency
        #Find fL which is just the maximum frequency
        #(19) Calculate L incudance
        #L = XL/2pif = XL/w
        L = XL/np.max(w)
    print('L: %.10f' % L)
      
    #Add Inductive Element
    ZoutECM += L*2*np.pi*f*1j
    
    
    #Find peaks being R and C elements
    #Zpeaks = (np.diff(np.sign(np.diff(-Zim))) < 0).nonzero()[0] + 1 # local max indices
    Zpeaks = (np.diff(np.sign(np.diff(-Zim)/np.diff(Zre))) < 0).nonzero()[0] + 1 # local max indices better? Think so
    N = len(Zpeaks)
    print("Number of Z Peaks = %d" % N)
    '''print("Z Peaks = ", Zpeaks)
    print("-Zim Peak Values = ", -Zim[Zpeaks])
    print("Zre Peak Values = ", -Zre[Zpeaks])'''
    
    #Warburg Element
    sig = Warburg(w, Zre, Zim)
    Zw = -sig/(np.sqrt(2*np.pi*f))*(1-1j)
    
    if N == 1:
        #(22, 23) Calculate Rct
        #0.5Rct = 0.5|Xct| = |Zimfctres| = |second[peak(-Zim)]|
        #Rct = 2*np.abs(-Zim[Zpeaks[0]])
        Rct = 2*np.abs(-Zim[Zpeaks[-1]])
        
        #Zpeak = Zpeaks[0] since theres only 1 peak
        print('Rct (R1): %.5f' % Rct)
        #Rctmid = 0.5*Rct+Rohm
        
        #wct is just the frequency of this peak
        #wct = w[Zpeaks[0]]
        wct = w[Zpeaks[-1]]
        
        #Calculate Cdl
        #C = 1/(Rw) (T = 1/w = RC)
        Cdl = 1/(Rct*wct)
        print('Cdl (C1): %.5f' % Cdl)
        #R3 = Rct and C3 = Cdl with R1C1=R2C2 being equivalent to R3C3 
        
        R1e = Rct
        C1e = Cdl
        T1e = R1e*C1e
        
        #Add R1C1
        Zc1 = 1/(C1e*2*np.pi*f*1j)
        
        ZoutECM += 1/(1/(R1e+Zw)+1/Zc1)
        
    elif N >= 2: #NEED MORE CONDITIONING FOR ROBUSTNESS, changed it to >=2
        #(20, 21) Calculate Rsei
        #0.5Rsei = 0.5|Xsei| = |Zimfseires| = |first[peak(-Zim)]|
        #Rsei = 2*np.abs(-Zim[Zpeaks[0]])
        Rsei = 2*np.abs(-Zim[Zpeaks[len(Zpeaks)-2]]) #takes second lowest frequency peak
        
        print('Rsei (R2): %.5f' % Rsei)
    
        #w is just the frequency of at the first peak
        #wsei = w[Zpeaks[0]]
        wsei = w[Zpeaks[len(Zpeaks)-2]] #takes second lowest frequency peak
        
        #Calculate Csei
        Csei = 1/(Rsei*wsei)
        print('Csei (C2): %.5f' % Csei)
        
        #(22, 23) Calculate Rct
        #0.5Rct = 0.5|Xct| = |Zimfctres| = |second[peak(-Zim)]|
        #Rct = 2*np.abs(-Zim[Zpeaks[1]])
        Rct = 2*np.abs(-Zim[Zpeaks[-1]]) #takes lowest frequency peak
        
        #Zpeak = Zpeaks[1] for the second peak
        print('Rct (R1/R1e): %.5f' % Rct)
        
        #wct is just the frequency of the first peak
        #wct = w[Zpeaks[1]]
        wct = w[Zpeaks[-1]] #takes lowest frequency peak
        
        #Calculate Cdl
        Cdl = 1/(Rct*wct)
        print('Cdl (C1/C1e): %.5f' % Cdl)
        
        R2e = Rsei
        C2e = Csei
        T2e = R2e*C2e
        
        R1e = Rct
        C1e = Cdl
        T1e = R1e*C1e
        
        #Add R1C1 and R2C2
        Zc1 = 1/(C1e*2*np.pi*f*1j)
        Zc2 = 1/(C2e*2*np.pi*f*1j)
        ZoutECM += 1/(1/R1e+1/Zc1) + 1/(1/(R2e+Zw)+1/Zc2)
    
    else: #No peaks assume R circuit
        Rct = 0
        Cdl = 0
    
    #output Rohm/R0, Rct/R1, Cdl/C1 and Z(f) of circuit model
    return Rohm, Rct, Cdl, ZoutECM

#Rinf is the min(Zre) where Zim < 0 input is (-Zim)
def Rinf(Zre, Zim):
    invalim = np.where(Zim > 0)[0] #get indices of invalid imaginary impedances
    return np.amin(np.delete(Zre, invalim)) #deletes invalid indices of Zre and finds minimum

#Add Warburg element? Returns sig = 0 
def Warburg(w, Zre, Zim):
    #Warburg
    #Determine R1
    #R1 is the maximum value of the real impedance
    R1s = np.max(Zre)
    #Find f1 (w1 = 2pif1)
    w1s = w[np.argmax(Zre)]
    
    #Take Arbitrary f2 (w2)
    #Take next lowest frequency
    #Aka take smallest freuqeuncy excluding w1
    w_arg = np.delete(w, np.argmax(Zre))
    w2s = np.min(w_arg)
    #Find R2
    R2s = float(Zre[np.argwhere(w == w2s)])
    
    y2 = -Zim[np.argmax(Zre)]
    y1 = -Zim[np.argwhere(w == w2s)]
    x2 = R1s #same as x value on Nyquist plot since its just the real part
    x1 = R2s
    al = math.atan2((y2 - y1),(x2 - x1))
    
    #sig = 1/Q =  Aw
    sig = (R1s - R2s)/(np.sqrt(w1s) - np.sqrt(w2s))
    #Q = 1/((R1s - R2s)/(np.power(w1s, al) - np.power(w2s, al))) #Warburg try testing later
    
    print('f1 = ' + str(w1s/(2*np.pi)) + ', f2 = ' + str(w2s/(2*np.pi)))
    print("Angle between low freq points: %.5f\u00B0" % math.degrees(al))
    print("Warburg Param Aw/\u03C3 = %.5f" % -sig)
    
    if 90 > math.degrees(al) and math.degrees(al) > 20:
        #Add Warburg (Semifinite)
        #return (1 - np.abs(math.degrees(al)-45)/math.degrees(al))*sig doesnt work well
        return sig
        #return 0
    else:
        return 0
