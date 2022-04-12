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

#Circle fit Package: https://github.com/marian42/circle-fit

@author: avramkachura

"""
import matplotlib.pyplot as plt 
import numpy as np
import math
import circle_fit as cf #least_squares_circle or hyper_fit  algos
import scipy

def ECM_EXT_C_v2(f, Zre, Zim):
    print('-----------' + 'ECM_EXT_C_v2' + '-----------')
     
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
    
    #RC Circuit Impedance for Warburg
    ZrcECM = np.zeros(len(ft), dtype = "complex_")
    
    #Find Rohm = Rinf, take minimum real as Rohm, for robustness if no local mins, maxes and warburg
    '''Rohm = np.min(Zre)
    print('Rohm (Rinf, R0): %.5f' % Rohm)
    
    #Add Ohmic Resistance
    ZoutECM += Rohm
    ZrcECM += Rohm'''
    
    #Find peaks being R and C elements
    Zpeaks = (np.diff(np.sign(np.diff(-Zim)/np.diff(Zre))) < 0).nonzero()[0] + 1 # local max indices
    Np = len(Zpeaks)
    print("Number of Z Peaks = %d" % Np)
    Ztroughs = (np.diff(np.sign(np.diff(-Zim)/np.diff(Zre))) > 0).nonzero()[0] + 1 # local min indices
    Ztroughs1 = Ztroughs
    
    Nt = len(Ztroughs)
    print("Number of Z Troughs = %d" % Nt)
    
    #Remove Warburg impedance points in ZimC and ZreC for circle fit
    #Placehold to remove warburg to circle fit, No
    ZreC = Zre
    ZimC = Zim
    ftC = ft
    
    #Find initial angled points
    angCi = math.degrees(math.atan2((-ZimC[-1] - -ZimC[-2]),(ZreC[-1] - ZreC[-2])))
    #Loop from low frequency points
    '''for i in range(len(ft)-1, -1, -1):
        #Calculate Angle of 
        angC = math.degrees(math.atan2((-ZimC[-1] - -ZimC[-2]),(ZreC[-1] - ZreC[-2])))
        
        #print(angC)
        
        #Remove angled component or until reach first local min index
        #AND statement because if there are no local mins for when angle is gone
        if Nt == 0 and (80 > angC and angC > 10):
            #Remove points
            ZreC = np.delete(ZreC, i)
            ZimC = np.delete(ZimC, i)
            ftC = np.delete(ftC, i) 
        #elif np.max(Ztroughs, initial = -1) < i: dont need to add initial b/c first check
        elif Nt != 0 and Ztroughs[-1] < i and angCi > 0: #and statement for robustness against 1 circuit with min, use initial angle
            #Remove points
            ZreC = np.delete(ZreC, i)
            ZimC = np.delete(ZimC, i)
            ftC = np.delete(ftC, i) 
        elif Nt == 0 and angCi < 0 and Np == 0: #for plain R circuits
            Rohm = np.min(Zre)
            print('Rohm (Rinf, R0): %.5f' % Rohm)
    
            #Add Ohmic Resistance
            ZoutECM += Rohm*np.ones(len(f))
            
            R1e = 0
            print('Rct (R1): %.5f' % R1e)
            
            C1e = 0
            print('Cdl (C1): %.5f' % C1e)
            
            #Exit function no warburg and RC component
            return Rohm, R1e, C1e, ZoutECM
        else:
            #Remove lowest frequency local min from troughs if reached
            if Nt != 0 and Ztroughs[-1] == i: 
                Ztroughs = np.delete(Ztroughs, len(Ztroughs) - 1)
            break'''
    
    #Circle 1 from second last local min
    #2 cases: No local mins, >1 local min with updated troughs
    if len(Ztroughs) == 0: #or len(Zpeaks) == 1:
        print('Fitting 1 Circle')
        #Mirror points 
        ZremC = np.hstack((ZreC, ZreC))
        ZimmC = np.hstack((ZimC, -ZimC))
        ZtmC = np.transpose(np.array([ZremC, -ZimmC]))
        
        #fit with circle xc, yc, radius, var
        xc1, yc1, r1, s1 = cf.hyper_fit(ZtmC)
        
        #Find Rohm = Rinf, take fitted center Zre - radius
        Rohm = xc1-r1
        print('Rohm (Rinf, R0): %.5f' % Rohm)
        
        #Add Ohmic Resistance
        ZoutECM += Rohm 
        ZrcECM += Rohm
        
        #Take Rct as the 2x peak height of the circle fit
        R1e = r1*2
        print('Rct (R1): %.5f' % R1e)
    
        #Linear interpolate freq of R1e using circle fit pts, include all points
        #np interpolate for boundaries
        #Find closest real impedance freq correspondance to xcenter of circle
        #ft is y, Zre is x
        fR1e = np.interp(xc1, Zre, ft)
        
        #Calculate Cdl
        C1e = 1/(2*np.pi*R1e*fR1e)
        print('Cdl (C1): %.5f' % C1e)
        
        #Add R1C1 to RC Circuit ECM
        Zc1 = 1/(C1e*2*np.pi*ft*1j)
        
        ZrcECM += 1/(1/R1e+1/Zc1) 
        
    #Fit 2 circles never do this?
    else:
        print('Fitting 2 Circles')
        #lowest frequency troughs
        ind = Ztroughs[-1]
        
        #R1C1 fit with low frequencies
        #Circle 1 impedance, include local min no (+1)
        ZreC1 = ZreC[ind:]
        ZimC1 = ZimC[ind:]
        
        #Mirror points 
        ZremC1 = np.hstack((ZreC1, ZreC1))
        ZimmC1 = np.hstack((ZimC1, -ZimC1))
        ZtmC1 = np.transpose(np.array([ZremC1, -ZimmC1]))
        
        #fit with circle xc, yc, radius, var
        xc1, yc1, r1, s1 = cf.hyper_fit(ZtmC1)
        
        #Take Rct as the 2x peak height of the circle fit
        R1e = r1*2
        print('Rct (R1): %.5f' % R1e)
    
        #Linear interpolate freq of R1e using all pts
        #fR1e = np.interp(xc1, Zre, ft) #no extrapolation
        f_interp = scipy.interpolate.interp1d(Zre, np.log(ft), fill_value = "extrapolate") #extrapolation
        fR1e = np.exp(f_interp(xc1))
        print('fdl (fRC1): %.5f' % fR1e)
        
        #Calculate Cdl
        C1e = 1/(2*np.pi*R1e*fR1e)
        print('Cdl (C1): %.5f' % C1e)
        
        
        #R2C2 fit with high frequencies
        #Circle 2 impedance, include local min, (+1)
        ZreC2 = ZreC[:ind+1]
        ZimC2 = ZimC[:ind+1]
        
        #Mirror points 
        ZremC2 = np.hstack((ZreC2, ZreC2))
        ZimmC2 = np.hstack((ZimC2, -ZimC2))
        ZtmC2 = np.transpose(np.array([ZremC2, -ZimmC2]))
        
        #fit with circle xc, yc, radius, var
        xc2, yc2, r2, s2 = cf.hyper_fit(ZtmC2)
         
        #Find Rohm = Rinf, take fitted with circle 2 center Zre - radius
        Rohm = xc2-r2
        
        if Rohm < 0:
            Rohm = np.min(Zre) #make it the minimum
        
        print('Rohm (Rinf, R0): %.5f' % Rohm)
        
        #Add Ohmic Resistance
        ZoutECM += Rohm 
        ZrcECM += Rohm
        
        #Take Rsei as the 2x peak height of the circle fit
        R2e = r2*2
        print('Rsei (R2): %.5f' % R2e)
    
        #Linear interpolate freq of R2e using all pts
        #fR2e = np.interp(xc2, Zre, ft)
        fR2e = np.exp(f_interp(xc2))
        print('fsei (fRC2): %.5f' % fR2e)
        
        #Calculate Csei
        C2e = 1/(2*np.pi*R2e*fR2e)
        print('Csei (C2): %.5f' % C2e)
        
        #Add R1C1 and R2C2 to RC Circuit ECM
        Zc1 = 1/(C1e*2*np.pi*ft*1j)
        Zc2 = 1/(C2e*2*np.pi*ft*1j)
        ZrcECM += 1/(1/R1e+1/Zc1) + 1/(1/R2e+1/Zc2)
    
    
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
    '''if len(Ztroughs) == 0 or len(Zpeaks) == 1:
        #Add 1RC elements found with Warburg
        Zc1 = 1/(C1e*2*np.pi*f*1j)
        ZoutECM += 1/(1/(R1e+Zw)+1/Zc1)
    else:
        #Add 2RC elements found with Warburg to Rct/R1 in series
        Zc1 = 1/(C1e*2*np.pi*f*1j)
        Zc2 = 1/(C2e*2*np.pi*f*1j)
        ZoutECM += 1/(1/(R1e+Zw)+1/Zc1) + 1/(1/R2e+1/Zc2)'''
        
    #No warburg
    if len(Ztroughs) == 0: #or len(Zpeaks) == 1:
        #Add 1RC elements found with Warburg
        Zc1 = 1/(C1e*2*np.pi*f*1j)
        ZoutECM += 1/(1/R1e+1/Zc1)
    else:
        #Add 2RC elements found with Warburg to Rct/R1 in series
        Zc1 = 1/(C1e*2*np.pi*f*1j)
        Zc2 = 1/(C2e*2*np.pi*f*1j)
        ZoutECM += 1/(1/R1e+1/Zc1) + 1/(1/R2e+1/Zc2)
        
        
    #plotting for debugging
    '''plt.figure(3)
    plt.plot(Zre, -Zim, label = 'Actual', marker = 'o', color  = 'purple', linestyle = 'None') 
    plt.plot(ZreC, -ZimC, marker = 'x', linestyle = 'None', label = 'C points')
    
    plt.plot(ZrcECM.real, -ZrcECM.imag, marker = 'x', color = 'blue', label = 'RC Circuit')
    plt.plot(ZoutECM.real, -ZoutECM.imag, marker = 'x', color = 'green', label = 'ECM')
    
    plt.plot(Zre[Ztroughs1],-Zim[Ztroughs1], 'kv', label = 'local min')
    plt.plot(Zre[Zpeaks],-Zim[Zpeaks], 'r^', label = 'local max')
    
    plt.xlabel('Zreal (\u03A9)')
    plt.ylabel('-Zim (\u03A9)')
    plt.grid(b=True, which='both', linestyle='-')
    plt.legend()
    #plt.figure(2)
    cf.plot_data_circle(ZremC1,-ZimmC1,xc1,yc1,r1)
    #plt.figure(3)
    cf.plot_data_circle(ZremC2,-ZimmC2,xc2,yc2,r2)'''
    
    #cf.plot_data_circle(ZremC,-ZimmC,xc1,yc1,r1)
    
    #output Rohm/R0, Rct/R1, Cdl/C1, Rsei/R2, Csei/C2, sig/Aw and Z(f) of circuit model
    #Output Rohm/R0, Rct/R1, Cdl/C1 and Z(f) of circuit model (parameters that are most important)
    #return Rohm, R1e, C1e, R2e, C2e, ZoutECM
    return Rohm, R1e, C1e, ZoutECM
