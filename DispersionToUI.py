# -*- coding: utf-8 -*-
"""
Created on Thu Nov  7 11:44:03 2019

@author: Vidar Flodgren
"""
import OFunc
import numpy as np
from OFunc import FFTD as FFTD
import matplotlib.pyplot as plt

from scipy.constants import c as c                     #m/s
from scipy.constants import epsilon_0 as eps_0         #m**-3kg**-1s**4A**2    - Permittivity of free space
UIRUN = 0
if UIRUN == 0:
    F_0      = np.linspace(0,2.5,1+int(( 2.5 / 0.05)))*10**10   #Vm**-1              - Optical Field
    ## -- Variable Parameters -- ##
    f_0x     = 375*10**12              #Hz                - Laser frequency (From our lab) 
    f_0y     = 375*10**12*2            #Hz                - Laser frequency (first Harmonic) 
    lambda_0 = c/f_0x                 # m
    w_0x     = 2*np.pi*f_0x            #Rad/s             - Laser Frequency
    w_0y     = 2*np.pi*f_0y            #Rad/s             - Laser Frequency
    Aeff     = 2.3*10**-12             # m**2              - Effective area 
    Ncycx    = 1.5                    # no unit          - Number of optical cycles in the pulse?s FWHM (Khurgin uses 1.7 in the paper)
    Ncycy    = 1.5                    # no unit          - Number of optical cycles in the pulse?s FWHM (Khurgin uses 1.7 in the paper)
    F_a      = 5.36*10**10             # Vm**-1            - Atomic Field
    L        = 0.1                    # mm               - Length of material dispersion
    
    #Calculate the Material Properties of your Selected Material#                       
    SiO2 = OFunc.MaterialProperties(lambda_0,'SiO2')
    SiO2.DispCoeff(L)
    
    ORD     = 5 
    BPF     = 0
    A_t      = 1                      # no unit           - Amplitude
    A_w      = 1                      # no unit           - Amplitude
    ## -- Defining the Electromagnetic Field Transient -- ##
    T_x = (Ncycx*2*np.pi)/w_0x             # seconds - Pulse Duration (downconverted)  
    T_y = (Ncycx*2*np.pi)/w_0x             # seconds - Pulse Duration (downconverted)  
    #Documentation claims that: ???t=4ln(2), that is W*T = 4*log(2) which means we get:
    W_x    = 4*np.log(2)/T_x
    W_y    = 4*np.log(2)/T_y
    t1   = 0; t2 = 15*T_x; t_0 = (t2+t1)/2
    N    = 100
    t    = np.linspace(t1,t2,N)             # Seconds - Time
    Dt   = (t2-t1)                       # Sampling interval
    dt   = Dt/N  
    w_t1 = 2*np.pi/Dt; w_t2 = 2*np.pi/dt   # Min-Max Frequency from time graph
    
    w1   = w_0x - w_0x*0.99; w2 = w_0x + 20*w_0x 
    w    = np.linspace(w1,w2,N)
    Dw   = w2-w1
    dw   = Dw/N
    Dt_w = 2*np.pi/w1; dt_w = 2*np.pi/w2   
    t_w1 = -Dt_w/2; t_w2 = Dt_w/2
    
    Delt1  = -100*10**(-15); Delt2 = 100*10**(-15)
    
    Delt = np.linspace(Delt1,Delt2,N)
    
"""
Old code that needs to be moved into the "Photoinduced only" simulation
## New a_fun code @(t,T,w_0x) that is more universal and can be integrated
a =  (exp(-2*log(2).*t.^2/T^2).*cos(w_0x*t))
a = fun_a(t,T_x,w_0x)

#We define a new function handle a_fun2n1 to help us integrate it over time with the new ^2n+1
fun_a2n1 = @(t,T,n,w_0x) (exp(-2*log(2).*t.^2/T.^2).*cos(w_0x*t)).^(2*n+1)

for n = 1:ORD
a2(n) = w_0x*trapz(t,fun_a2n1(t,T_x,n,w_0x))
end
"""

## Define separate Et and Ew series based on input variables
theta = OFunc.Phase_Offset(0)
phiw = OFunc.Dispersion_Factor(SiO2.phi,w,w_0x)

Et = OFunc.TDGE(A_t,t,t_0,T_x,w_0x,theta) 
Ew = OFunc.FDGE(A_w,w,W_x,w_0x,phiw)

#Note: These are related to each other by: \delta w*\delta t =4ln(2), that is W*T = 4*log(2) which means we get: W  = 4*log(2)/T
#T = (Ncycx*2*pi())/w_0x => W = 4*log(2)*w_0x/(Ncycx*2*pi()) 


## if PlotSelect == "Non Fourier Et(t)"
# plot(t-t_0,real(Et)); xlabel('time [s]'); ylabel('E_t(t)'); fprintf(['\n','Plotting [',PlotSelect{1},']'])
# end

##-- Time to test conversions between the two --##
## -- E_t -> FFT = E_w -> Band Filter + IFFT = E_t -- ## 

Estr = OFunc.FFTD(t,Et,w_0x,SiO2.phi) #This generates:  | Estr.ttt.Et & Estr.ttt.t | Estr.ttt.Ew & Estr.ttf.w | Estr.ttt.Etdisp & Estr.ttt.tdisp|
Estr.ttt(theta)

#%%

##Bandpass Filter
if BPF > 0:

    for n = 1:length(Estr.ttt.Ew)
        if abs(Estr.ttt.Ew(n))<max(abs(Estr.ttt.Ew))/BPF
            Estr.ttt.Ew(n) = 0
 

    Estr2 = FFTD(Estr.wf,Estr.Ew,w_0,SiO2.phi,0)
    Estr2.ftt(theta)
    Estr = Estr2
    del(Estr2)



## if PlotSelect == "E_t FFT Plot"
#plot(Estr.ttt.w,abs(Estr.ttt.Ew)); grid on; fprintf(['\n','Plotting [',PlotSelect{1},']']);
#end
##


## if PlotSelect == "E_t Final Applied Dispersion"
#plot(Estr.ttt.tdisp,real(Estr.ttt.Etdisp));grid on; fprintf(['\n','Plotting [',PlotSelect{1},']']);
#end

## -- Post Dispersion <a^2n+1> and Photoinduced Charge -- ##
#After we get back Etifft, we now have a new value for a(t), and thus need to calculate a new <a^2n+1>

## -- Calculating <a^2n+1> for a dispersed pulse -- ##

a2disp = OFunc.Vec_Pot_Mom(w_0x,Estr.tf,Estr.Etf,ORD) 
## -- Equation 17 -- ##
#fun_Q2 = @(F_0x,F_a,a2disp,Aeff) eps_0*F_0x.*(F_0x/F_a).^2.*(a2disp(1) +(F_0x/F_a).^2*a2disp(2)...
#                                   +(F_0x/F_a).^4*a2disp(3)+(F_0x/F_a).^6*a2disp(4)+(F_0x/F_a).^8*a2disp(5))*Aeff
#Q = fun_Q2(F_0,F_a,a2disp,Aeff)
Q = fun_Q(F_0,F_a,a2disp,Aeff,ORD)  


Etf2n(1,:) = real(Estr.ttt.Etdisp)
for n = 1:4
Etf2n(n+1,:) = real(Estr.ttt.Etdisp).^(2*n+1)
end
## if PlotSelect == "E_t(t)^2n+1 After Dispersion"
#plot(Estr.ttt.tdisp,Etf2n(:,:),'LineWidth',1.3)
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
# legend('Ncycx = 1.0', 'Ncycx = 1.1', 'Ncycx = 1.2', 'Ncycx = 1.3', 'Ncycx = 1.4', 'Ncycx = 1.5', 'Ncycx = 1.6', 'Ncycx = 1.7', 'Ncycx = 1.8', 'Ncycx = 1.9', 'Ncycx = 2.0', 'e')
#end





## -- E_w -> IFFT = E_t, just to show that the dispersion applied is equivalent 

#Redundant Code
#phiw  = fun_phiw(SiO2.phi(1),SiO2.phi(2),SiO2.phi(3),SiO2.phi(4),w,w_0x); # Assigning the dispersion component constants, remember that SiO2.phi(1) represents \phi_0
#Ew    = fun_Ew(A_w,w,W,w_0x,phiw);   # Using the Ew function to define the frequency space dispersed wave
#nfft  = 2^nextpow2(length(Ew));     # Determining required length (nfft) of frequency axis: Padding by 2^P, where 2^(P-1)<length(Ew) and 2^P>lenght(Ew) (if length(Ew=1000), then P = 10 for instance since 2^10 = 1024>1000
#Ewfft = fftshift(ifft(Ew,nfft));    # Performing the fourier transform using nfft as the N=nfft number of points 
#PHITEXT = sprintf(['\\phi_0 = ', num2str(SiO2.phi(1),'#.4g'),'\n',...
#                   '\\phi_1 = ', num2str(SiO2.phi(2),'#.4g'),'\n',...
#                   '\\phi_2 = ', num2str(SiO2.phi(3),'#.4g'),'\n',...
#                   '\\phi_3 = ', num2str(SiO2.phi(4),'#.4g'),'\n']);
               
## Alternative Code ##
Eftt  = FFTD(w,Ew,w_0x,'ftt',SiO2.phi,0)

 

## Alternate Code
waitbar(0.45,PROGBAR,'Calculating $E_w(t)$')
Estw = FFTD(w,Ew,w_0x,'ftt',SiO2.phi,0)
## if PlotSelect ==  "E_w IFFT Dispersion"
#plot(Estw.ftt.tdisp,real(Estw.ftt.Etdisp)); grid on; fprintf(['\n','Plotting [',PlotSelect{1},']']);
#end

## -- Dual Pulses with Temporal Delay Delt -- ##
#To start off, the temporal delay \Delta t (We will call it Delt in the
#code) is the delay between two short pulses. Namely:
#The driving pulse (F_0 = Fx0) that we will call: Ed(t) (The oscillation direction is between plates)
#The Injection Pulse (F0y) that we will call:     Ei(t) (The oscillation direction is perpendicular to the plates)
#For all intents and purposes, these are called a(t), and were called this
#in reference in the photoinducedchargecode.
#The equation for the new vector potential <a^{2n+1}>(Delt) is now given by:

#<a^{2n+1}>(Delt) = w_0x integral(@(time) a(t)*(a(t-Delt))^2n,t(1),t(N))
#Note here, that a(t) represents the Driving pulse (F_0) and a^2n represents the Injection pulse F0y


#We start by defining the two pulses (To start off, let's make them identical)
##fun_Et(A_t,t,T,w_0x,theta)
##fun_Ew(A_w,w,W,w_0x,psi) - You need this->Et to do this I think...
#w_d = w_0x; w_i = w_0y; # w_d, the weaker driving pulse (w_0x or w_0xx), and w_i, the stronger injection pulse (w_0y, or w_i)
#W_d = W  ; W_i = W_y; #Respective FWHM for each




## Dual Polarisation Equation ##
#fun_QDelt2 = @(F_0x,F_0y,a2n,Aeff) eps_0.*F_0x.*(F_0y/F_a).^2*(1/3 * a2n(1)+...
#                                                              1/5 * a2n(2).*(F_0y./F_a).^2+...
#                                                              1/7 * a2n(3).*(F_0y./F_a).^4+...
#                                                              1/9 QPol2(n) = fun_QDelt(F_0(end)/12.5,F_0(end)/1.25,F_a,a2pol(n,:),Aeff,ORD)  * a2n(4).*(F_0y./F_a).^6+...
#                                                              1/11* a2n(5).*(F_0y./F_a).^8)*Aeff; 
                                                     


N_Delt = Delt2 - Delt1 + 1
Deltn  = linspace(Delt1,Delt2,N_Delt)
E_td = FFTD(t,fun_Et(A_t,t,t_0,T_x,w_0x,0),w_0x,'ttt',SiO2.phi,0)
E_ti = FFTD(t,fun_Et(A_t,t,t_0,T_y,w_0y,0),w_0y,'ttt',SiO2.phi,0)
Delt = linspace(Delt1*dt,Delt2*dt,N_Delt)
waitbar(0.5+0.5*0/N,PROGBAR,['Integrating to determine $\left<a^2n+1\right>$','   n/N$\_$Delt = ',num2str(0),'/',num2str(N_Delt)])

for n = 1:N_Delt
if rem(n,round(N_Delt/100,-1)) == 0
waitbar(0.5+0.5*n/N_Delt,PROGBAR,['Integrating to determine $\left<a^2n+1\right>$','   n/N$\_$Delt = ',num2str(n),'/',num2str(N_Delt)])
#if N/n = A*N/100 1000/200
end
E_td.ttt.Etdshift = circshift(E_td.ttt.Etdisp,Deltn(n)); E_td.ttt.tdshift = circshift(E_td.ttt.tdisp,Deltn(n)) 
Edires = E_td.ttt.Etdshift+E_ti.ttt.Etdisp

for m = 1:ORD
    a2harm(n,m) = w_0x*trapz(E_td.ttt.tdisp,real(Edires).^(2*m+1))
    a2pol(n,m)  = w_0x*trapz(E_td.ttt.tdisp,real(E_td.ttt.Etdshift).*real(E_ti.ttt.Etdisp).^(2*m)) 
end


QHarm(n) = fun_Q(F_0(end),F_a,a2harm(n,:),Aeff,ORD)
QPol(n)  = fun_QDelt(F_0x,F_0y,F_a,a2pol(n,:),Aeff,ORD)
 
## In case you want animation of the transformation here :)

anim = 1
if anim == 0
    figure(67); clf;
subplot(3,1,1)
plot(E_ti.ttt.tdisp,real(E_ti.ttt.Etdisp))
ylim([-1.3 1.3])
hold on
plot(E_td.ttt.tdisp,real(E_td.ttt.Etdshift))
ylim([-1.3 1.3])
subplot(3,1,2)
plot(E_td.ttt.tdisp,real(Edires(n,:)))
ylim([-1.3 1.3])
subplot(3,1,3)
plot(Delt(1:n),QHarm(1:n))
drawnow

figure(68)
plot(real(E_td.ttt.Etdshift))
drawnow
end
end


## if PlotSelect == "Temporal Overlap Induced Charge"
#plot(Delt,real(QHarm(:,:))); fprintf(['\n','Plotting [',PlotSelect{1},']']);
#legend(['T_E_i =' num2str((0.5*n-10)*T+T,'#.4g') ')'])
#end

## if PlotSelect == "Dual Polarised Induced Charge"
#plot(Delt,real(QPol(:,:)));fprintf(['\n','Plotting [',PlotSelect{1},']']);
#legend(['T_E_i =' num2str((0.5*n-10)*T+T,'#.4g') ')'])
#end

#Old code used to plot all of this IN-SCRIPT. Now we plot it in a separate
#script to make it easier to change/add things (without having to do
#everything twice

waitbar(1,PROGBAR,'DONE!'); #pause(0.5);
delete(PROGBAR)
evalin('base','UIGraphPlotter')



THESISGRAPHS = 0
if THESISGRAPHS == 1 #This is for when you wish to reproduce graphs for your paper, keep these constant so you don't need to type them in again

## Plotting different phi_n dispersion graphs, we do this to show exactly how the different components interact
fg = figure(6); clf; FIG.Name = 'E_w->E_t For Different \phi_n'; fg.Position = [550 270 700 520]
     #Note, the last number helps 
for n = 5
    DISP{1} = [0 0 0 0 0] 
    DISP{2} = [pi 0 0 0 1]
    DISP{3} = [0 20*(10^-15) 0 0 2]
    DISP{4} = [0 0 25*(10^-15)^2 0 3]
    DISP{5} = [0 0 0 60*(10^-15)^3 4]
DISPFAC = DISP{n}
DISPDISP = [pi 20 25 60]
Estr2 = FFTD(t,Et,w_0x,'ttt',DISPFAC,0)
if n == 1
PHITEXT{1} =      ['$\Phi = 0$']
elseif n == 2
    PHITEXT{n} =      ['$\phi_',num2str(n-2),' = ', num2str(DISPDISP(n-1),'%.4g'), '$']
else
    PHITEXT{n} =      ['$\phi_',num2str(n-2),' = ', num2str(DISPDISP(n-1),'%.4g'), '$ fs$^',num2str(n-1),' $']
end

cla reset
plot(Estr2.ttt.tdisp/10^-15,real(Estr2.ttt.Etdisp)); fprintf(['\n','Plotting [',PlotSelect{1},']'])
legend(PHITEXT{n},'Location','northwest')
xlabel('Time [fs]')
ylabel('Amplitude')
grid on
set(gca,'fontsize', 20)
set(gca,'Box','on')
LEFT = 0.13; BOTTOM = 0.15; RIGHT = 0.02; TOP = 0.05
InSet = [LEFT BOTTOM RIGHT TOP] #Fontsize 12
set(gca, 'Position', [InSet(1:2), 1-InSet(1)-InSet(3), 1-InSet(2)-InSet(4)])
xticks('auto')
Legg=legend(PHITEXT{n},'Location','northwest')
set(Legg,'FontSize',18)
axis tight
ylim([-1 1])
end    
end
