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
import math


def ECM_EXT_v3(f, Zre, Zim):
    #round to 15 decimals 32 bit -> 10 decimals (digits)
    Zim = np.round(Zim, 15)
    Zre = np.round(Zre, 15)
    
    print('-----------' + 'ECM_EXT_v3' + '-----------')
     
    #Output Impedance from ECM
    ZoutECM = np.zeros(len(f), dtype = "complex_")
    
    #Filter out Zre<
    invalre = np.where(Zre < 0)[0]
    Zre = np.delete(Zre, invalre)
    Zim = np.delete(Zim, invalre)
    ft = np.delete(f, invalre) #dont delete from f since need to generate ECM with same size and same spectrum
    #use ft to analyze filtered out data and f for ZoutECM
    
    #L inductive element
    L = 0
    XL = np.abs(np.min(-Zim)) 
    if XL >= 0 and XL == Zim[np.argmax(ft)]: #robust to frequencies that arent the maximum frequency
        #Find fL which is just the maximum frequency
        #(19) Calculate L incudance
        #L = XL/2pif = XL/w
        L = XL/(2*np.pi*np.max(ft))
    print('L: %.10f' % L)
      
    #Add Inductive Element
    ZoutECM += L*2*np.pi*f*1j
    
    #Removes 0<Zim based on last frequency
    if np.size(np.where(Zim > 0)[0]) != 0: #if no -Zim>0 measurements
        filim = np.max(np.where(Zim > 0)[0]) #take last element 
        Zre = np.delete(Zre, np.arange(0,filim+1))
        Zim = np.delete(Zim, np.arange(0,filim+1))
        ft = np.delete(ft, np.arange(0,filim+1)) 
    
    #RC Circuit Impedance
    ZrcECM = np.zeros(len(ft), dtype = "complex_")
    
    #Find Rohm = Rinf, take minimum real as Rohm
    Rohm = np.min(Zre)
    print('Rohm (Rinf, R0): %.5f' % Rohm)
    
    #Add Ohmic Resistance
    ZoutECM += Rohm
    ZrcECM += Rohm
    
    #Find peaks being R and C elements
    Zpeaks = (np.diff(np.sign(np.diff(-Zim)/np.diff(Zre))) < 0).nonzero()[0] + 1 # local max indices
    N = len(Zpeaks)
    print("Number of Z Peaks = %d" % N)
    Ztroughs = (np.diff(np.sign(np.diff(-Zim)/np.diff(Zre))) > 0).nonzero()[0] + 1 # local min indices
    Nt = len(Ztroughs)
    print("Number of Z Troughs = %d" % Nt)
    
    
    if N == 1:
        #(22, 23) Calculate Rct
        #0.5Rct = 0.5|Xct| = |Zimfctres| = |second[peak(-Zim)]|
        Rct = 2*np.abs(-Zim[Zpeaks[-1]])
        print('Rct (R1): %.5f' % Rct) #since theres only 1 peak
        
        #wct = 2piwct is just the frequency of this peak
        fct = ft[Zpeaks[-1]]
        
        #Calculate Cdl
        #C = 1/(Rw) (T = 1/w = 1/(2pif) = RC)
        Cdl = 1/(Rct*2*np.pi*fct)
        print('Cdl (C1): %.5f' % Cdl)
        #R3 = Rct and C3 = Cdl with R1C1=R2C2 being equivalent to R3C3 
        
        R1e = Rct
        C1e = Cdl
        T1e = R1e*C1e
        
        #Add R1C1 to RC Circuit ECM
        Zc1 = 1/(C1e*2*np.pi*ft*1j)
        
        ZrcECM += 1/(1/R1e+1/Zc1)
        
    elif N >= 2: #NEED MORE CONDITIONING FOR ROBUSTNESS, changed it to >=2
        #(20, 21) Calculate Rsei
        #0.5Rsei = 0.5|Xsei| = |Zimfseires| = |first[peak(-Zim)]|
        Rsei = 2*np.abs(-Zim[Zpeaks[-2]]) #takes second lowest frequency peak
        print('Rsei (R2): %.5f' % Rsei)
    
        #w = 2pif is just the frequency of at the first peak
        fsei = ft[Zpeaks[-2]] #takes second lowest frequency peak
        
        #Calculate Csei
        Csei = 1/(Rsei*2*np.pi*fsei)
        print('Csei (C2): %.5f' % Csei)
        
        #(22, 23) Calculate Rct
        #0.5Rct = 0.5|Xct| = |Zimfctres| = |second[peak(-Zim)]|
        Rct = 2*np.abs(-Zim[Zpeaks[-1]]) #takes lowest frequency peak
        print('Rct (R1/R1e): %.5f' % Rct)
        
        #wct is just the frequency of the first peak
        fct = ft[Zpeaks[-1]] #takes lowest frequency peak
        
        #Calculate Cdl
        Cdl = 1/(Rct*2*np.pi*fct)
        print('Cdl (C1/C1e): %.5f' % Cdl)
        
        R2e = Rsei
        C2e = Csei
        T2e = R2e*C2e
        
        R1e = Rct
        C1e = Cdl
        T1e = R1e*C1e
        
        #Add R1C1 and R2C2 to RC Circuit ECM
        Zc1 = 1/(C1e*2*np.pi*ft*1j)
        Zc2 = 1/(C2e*2*np.pi*ft*1j)
        ZrcECM += 1/(1/R1e+1/Zc1) + 1/(1/R2e+1/Zc2)
    
    else: #No peaks assume R circuit
        #CT and DL set to 0
        R1e = 0
        C1e = 0
        #SEI set to 0    
        R2e = 0
        C2e = 0
        
        #ZrcECM is set to 0s so warburg is robust
        #Return statement to handle no peaks and avoid warburg, need it with warburg? No?
        #return Rohm, Rct, Cdl, ZoutECM No b/c can be warburg
    
    #Warburg Parameter Estimation
    print('Warburg')
    
    #Impedance to generate Warburg Estimate 
    ZreW = Zre - ZrcECM.real
    ZimW = Zim - ZrcECM.imag
    
    #Two Lowest Frequencies Impedance Angle
    print('f1 = ' + str(ft[-1]) + ', f2 = ' + str(ft[-2]))
    ang = math.degrees(math.atan2((-ZimW[-1] - -ZimW[-2]),(ZreW[-1] - ZreW[-2])))
    print("Angle between f1 and f2 points: %.5f\u00B0" % ang)
    
    #Calculate Warburg paramter
    sig = (ZreW[-1] - ZreW[-2])/(1/np.sqrt(2*np.pi*ft[-1]) - 1/np.sqrt(2*np.pi*ft[-2]))
    
    if 90 > ang and ang > 0:
        #Warburg only applies a max 45 degree angle for low frequencies 
        if 45 > ang:
            sig = sig*ang/45 #scale warburg lower for less impact on low frequencies
    #Only apply Warburg if there is a clear component
    else:
        sig = 0
    print("Warburg Param Aw/\u03C3 = %.7f" % sig)
    Zw = sig*(np.sqrt(2/(1j*np.pi*f*2)))
    
    #Output Impedance from ECM with known Warburg  
    if N==1:
        #Add 1RC elements found with Warburg
        Zc1 = 1/(C1e*2*np.pi*f*1j)
        ZoutECM += 1/(1/(R1e+Zw)+1/Zc1)
    elif N >= 1:
        #Add 2RC elements found with Warburg to Rct/R1 in series
        Zc1 = 1/(C1e*2*np.pi*f*1j)
        Zc2 = 1/(C2e*2*np.pi*f*1j)
        ZoutECM += 1/(1/(R1e+Zw)+1/Zc1) + 1/(1/R2e+1/Zc2)
        
    #output Rohm/R0, Rct/R1, Cdl/C1, Rsei/R2, Csei/C2, sig/Aw and Z(f) of circuit model
    #Output Rohm/R0, Rct/R1, Cdl/C1 and Z(f) of circuit model (parameters that are most important)
    return Rohm, R1e, C1e, ZoutECM