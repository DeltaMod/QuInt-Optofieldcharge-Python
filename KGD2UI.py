# -*- coding: utf-8 -*-
"""
Created on Thu Nov  7 11:44:03 2019

@author: Vidar Flodgren
"""
import json


import OFunc
import numpy as np
from OFunc import FFTD as FFTD

from scipy.constants import c as c                     #m/s
from scipy.constants import epsilon_0 as eps_0         #m**-3kg**-1s**4A**2    - Permittivity of free space

def Calculate_Dispersion(D,mat):
    D['F_0x']      = np.linspace(0,2.5,1+int(( 2.5 / 0.05)))*10**10   #Vm**-1              - Optical Field
    ## -- Variable Parameters -- ##
    D['f_0x']     = 375*10**12              #Hz                - D['L']aser frequency (From our lab) 
    D['f_0y']     = 375*10**12*2            #Hz                - D['L']aser frequency (first Harmonic) 
    D['lambda_0'] = c/D['f_0x']                 # m
    D['w_0x']     = 2*np.pi*D['f_0x']            #Rad/s             - D['L']aser Frequency
    D['w_0y']     = 2*np.pi*D['f_0y']            #Rad/s             - D['L']aser Frequency
    D['Aeff']     = 2.3*10**-12             # m**2              - Effective area 
    D['Ncycx']    = 1.5                    # no unit          - Number of optical cycles in the pulse?s FWHM (Khurgin uses 1.7 in the paper)
    D['Ncycy']    = 1.5                    # no unit          - Number of optical cycles in the pulse?s FWHM (Khurgin uses 1.7 in the paper)
    D['F_a']      = 5.36*10**10             # Vm**-1            - Atomic Field
    D['L']        = 0.1                    # mm               - D['L']ength of material dispersion
    
    #Calculate the Material Properties of your Selected Material#          
    D[mat] = mat             
    D[mat] = OFunc.MaterialProperties(D['lambda_0'],mat)
    D[mat].DispCoeff(D['L'])
    
    
    D['ORD']     = 5 
    D['BPF']     = 0
    D['A_t']      = 1                      # no unit           - Amplitude
    D['A_w']      = 1                      # no unit           - Amplitude
    ## -- Defining the Electromagnetic Field Transient -- ##
    D['T_x'] = (D['Ncycx']*2*np.pi)/D['w_0x']             # seconds - Pulse Duration (downconverted)  
    D['T_y'] = (D['Ncycx']*2*np.pi)/D['w_0x']             # seconds - Pulse Duration (downconverted)  
    #Documentation claims that: ???t=4ln(2), that is W*T = 4*log(2) which means we get:
    D['W_x']    = 4*np.log(2)/D['T_x']
    D['W_y']    = 4*np.log(2)/D['T_y']
    D['t1']  = 0; D['t2'] = 15*D['T_x']; D['t_0'] = (D['t2']+D['t1'])/2
    D['N']    = 100
    D['t']    = np.linspace(D['t1'],D['t2'],D['N'])             # Seconds - Time
    Dt   = (D['t2']-D['t1'])                       # Sampling interval
    dt   = Dt/D['N']  
    D['w_t1'] = 2*np.pi/Dt; D['w_t2'] = 2*np.pi/dt   # Min-Max Frequency from time graph
    
    D['w1']   = D['w_0x'] - D['w_0x']*0.99; D['w2'] = D['w_0x'] + 20*D['w_0x'] 
    D['w']    = np.linspace(D['w1'],D['w2'],D['N'])
    D['Dw']   = D['w2']-D['w1'] 
    D['dw']   = D['Dw']/D['N']
    D['Dt_w'] = 2*np.pi/D['w1']; D['dt_w'] = 2*np.pi/D['w2']   
    D['t_w1'] = -D['Dt_w']/2; D['t_w2'] = D['Dt_w']/2
    
    D['Delt1']  = -50; D['Delt2'] = 50 #Note, this used to be given in fs already
    
    D['Delt'] = np.linspace(D['Delt1'],D['Delt2'],D['N'])
        
    
    #%%  Define separate Et and Ew series based on input variables #%% 
    D['theta'] = OFunc.Phase_Offset(0)
    D['phiw'] = OFunc.Dispersion_Factor(D[mat].phi,D['w'],D['w_0x'])
    
    f = open("SessionRestore.txt","w")
    f.write( str(D) )
    f.close()
    R
    Et = OFunc.TDGE(D['A_t'],D['t'],D['t_0'],D['T_x'],D['w_0x'],D['theta']) 
    Ew = OFunc.FDGE(D['A_w'],D['w'],D['W_x'],D['w_0x'],D['phiw'])
    
    #Note: These are related to each other by: \delta w*\delta t =4ln(2), that is W*T = 4*log(2) which means we get: W  = 4*log(2)/T
    #T = (D['Ncycx']*2*pi())/D['w_0x'] => W = 4*log(2)*D['w_0x']/(D['Ncycx']*2*pi()) 
    
    
    ## if PlotSelect == "Non Fourier Et(t)"
    # plot(t-D['t_0'],real(Et)); xlabel('time [s]'); ylabel('E_t(t)'); fprintf(['\n','Plotting [',PlotSelect{1},']'])
    # end
    
    ## -- E_t -> FFT = E_w -> Band Filter + IFFT = E_t -- ## 
    Estr = OFunc.FFTD(D['t'],Et,D['w_0x'],D[mat].phi) #This generates:  | Estr.ttt.Et & Estr.ttt.t | Estr.ttt.Ew & Estr.ttf.w | Estr.ttt.Etdisp & Estr.ttt.tdisp|
    Estr.ttt(D['theta'])
    
    #%%
    
    ##Bandpass Filter
    if D['BPF'] > 0:
    
        for n in range(len(Estr.Ew)):
            if abs(Estr.ttt.Ew(n))<max(abs(Estr.ttt.Ew))/D['BPF']:
                Estr.Ew[n] = 0
     
        Estr2 = FFTD(Estr.wf,Estr.Ew,D['w_0x'],D[mat].phi)
        Estr2.ftt(D['theta'])
        Estr = Estr2
        del(Estr2)
    
    ## if PlotSelect == "E_t FFT Plot"
    #plot(Estr.ttt.w,abs(Estr.ttt.Ew)); grid on; fprintf(['\n','Plotting [',PlotSelect{1},']']);
    #end
    
    ## if PlotSelect == "E_t Final Applied Dispersion"
    #plot(Estr.ttt.tdisp,real(Estr.ttt.Etdisp));grid on; fprintf(['\n','Plotting [',PlotSelect{1},']']);
    #end
    
    #%% -- Post Dispersion <a^2n+1> and Photoinduced Charge -- #%%
    #After we get back Etifft, we now have a new value for a(t), and thus need to calculate a new <a^2n+1>
    ## -- Calculating <a^2n+1> for a dispersed pulse -- ##
    a2disp = OFunc.Vec_Pot_Mom(D['w_0x'],Estr.tf,Estr.Etf,D['ORD']) 
    ## -- Equation 17 -- ##
    #fun_Q2 = @(D['F_0x'],D['F_a'],a2disp,D['Aeff']) eps_0*D['F_0x'].*(D['F_0x']/D['F_a'])**2.*(a2disp(1) +(D['F_0x']/D['F_a'])**2*a2disp(2)...
    #                                   +(D['F_0x']/D['F_a'])**4*a2disp(3)+(D['F_0x']/D['F_a'])**6*a2disp(4)+(D['F_0x']/D['F_a'])**8*a2disp(5))*D['Aeff']
    #Q = fun_Q2(F_0,D['F_a'],a2disp,D['Aeff'])
    Q = OFunc.Photoinduced_Charge(D['F_0x'],D['F_a'],a2disp,D['Aeff'],D['ORD'])
    Etf2n = [np.real(Estr.Etf)]+[np.real(Estr.Etf)**(2*n+1) for n in range(D['ORD'])]
    
    ## if PlotSelect == "E_t(t)^2n+1 After Dispersion"
    #plot(Estr.ttt.tdisp,Etf2n(:,:),'D['L']ineWidth',1.3)
    #end
    
    
    ## if PlotSelect == "Post Dispersion Photoinduced Charge"
    #plot(F_0,abs(Q))
    #figure(66)
    #hold on
    #plot(F_0,abs(Q))
    #grid on
    #xlabel('Optical Field [Vm^-^1]')
    #ylabel('Photinduced Charge [C]')
    # set(gca, 'YScale', 'log')
    # plot([0,2.5*10^10],    [1.6*10^-19 1.6*10^-19]'r--')
    # legend('D['Ncycx'] = 1.0', 'D['Ncycx'] = 1.1', 'D['Ncycx'] = 1.2', 'D['Ncycx'] = 1.3', 'D['Ncycx'] = 1.4', 'D['Ncycx'] = 1.5', 'D['Ncycx'] = 1.6', 'D['Ncycx'] = 1.7', 'D['Ncycx'] = 1.8', 'D['Ncycx'] = 1.9', 'D['Ncycx'] = 2.0', 'e')
    #end
    
    
    
    
    
    ## -- E_w -> IFFT = E_t, just to show that the dispersion applied is equivalent 
    Estw = OFunc.FFTD(D['w'],Ew,D['w_0x'],D[mat].phi)
    Estw.ftt(D['theta'])
    ## if PlotSelect ==  "E_w IFFT Dispersion"
    #plot(Estw.ftt.tdisp,real(Estw.ftt.Etdisp)); grid on; fprintf(['\n','Plotting [',PlotSelect{1},']']);
    #end
    
    ## -- Dual Pulses with Temporal Delay Delt -- ##
    #The driving pulse (field D['F_0x']) Ed(t)  Ed(t) (Oscillation direction is between plates)
    #And the Injection Pulse (field F_y)  Ei(t) (Oscillation direction is perpendicular to the plates)
    #Here: a(t) represents the Driving pulse (D['F_0x']) and a^2n represents the Injection pulse (F_y
    
    #w_d = D['w_0x']; w_i = D['w_0y']; # w_d, the weaker driving pulse (D['w_0x']), and w_i, the stronger injection pulse (D['w_0y'])
    #W_d = W  ; W_i = D['W_y']; #Respective FWHM for each
    
    ## Dual Polarisation Equation ##
    #fun_QD['Delt2'] = @(D['F_0x'],F_y,a2n,D['Aeff']) eps_0.*D['F_0x'].*(F_y/D['F_a'])**2*(1/3 * a2n(1)+...
    #                                                              1/5 * a2n(2).*(F_y./D['F_a'])**2+...
    #                                                              1/7 * a2n(3).*(F_y./D['F_a'])**4+...
    #                                                              1/9 QPol2(n) = fun_QDelt(F_0(end)/12.5,F_0(end)/1.25,D['F_a'],a2pol(n,:),D['Aeff'],D['ORD'])  * a2n(4).*(F_y./D['F_a'])**6+...
    #                                                              1/11* a2n(5).*(F_y./D['F_a'])**8)*D['Aeff']; 
                                                         
    D['N_Delt'] = int(D['Delt2'] - D['Delt1'] + 1)
    Deltn  = np.linspace(D['Delt1'],D['Delt2'],D['N_Delt']).round().astype(int)
    E_td = OFunc.FFTD(D['t'],OFunc.TDGE(D['A_t'],D['t'],D['t_0'],D['T_x'],D['w_0x'],0),D['w_0x'],D[mat].phi)
    E_td.ttt(0)
    E_ti = OFunc.FFTD(D['t'],OFunc.TDGE(D['A_t'],D['t'],D['t_0'],D['T_y'],D['w_0y'],0),D['w_0y'],D[mat].phi)
    E_ti.ttt(0)
    D['Delt'] = np.linspace(D['Delt1']*dt,D['Delt2']*dt,D['N_Delt'])
    QHarm = []
    QPol = []
    a2harm =  np.zeros([D['N_Delt'],D['ORD']])
    a2pol  =  np.zeros([D['N_Delt'],D['ORD']])
    for n in range(D['N_Delt']):
        Etdshift = np.roll(E_td.Etf,Deltn[n])
        tdshift  = np.roll(E_td.tf,Deltn[n]) 
        Edires = E_td.Etf+E_ti.Etf
    
        for m in range(D['ORD']):
            a2harm[n,m] = D['w_0x']*np.trapz(E_td.tf,np.real(Edires)**(2*m+1))
            a2pol[n,m]  = D['w_0x']*np.trapz(E_td.tf,np.real(Etdshift)*np.real(E_ti.Etf)**(2*m)) 
    
    #Consider restructuring this command to include all inside of Q.Sum[n] instead of Q[n].Sum, it will make your life easier.
        QHarm.append(OFunc.Photoinduced_Charge(D['F_0x'][-1],D['F_a'],a2harm[n,:],D['Aeff'],D['ORD']))
        QPol.append(OFunc.Delta_Photoinduced_Charge(D['F_0x'][-1],D['F_a'],D['F_a'],a2pol[n,:],D['Aeff'],D['ORD']))
    return() 
    ## if PlotSelect == "Temporal Overlap Induced Charge"
    #plot(Delt,real(QHarm(:,:))); fprintf(['\n','Plotting [',PlotSelect{1},']']);
    #legend(['T_E_i =' num2str((0.5*n-10)*T+T,'#.4g') ')'])
    #end
    
    ## if PlotSelect == "Dual Polarised Induced Charge"
    #plot(Delt,real(QPol(:,:)));fprintf(['\n','Plotting [',PlotSelect{1},']']);
    #legend(['T_E_i =' num2str((0.5*n-10)*T+T,'#.4g') ')'])
    #end
