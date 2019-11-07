# -*- coding: utf-8 -*-
"""
FFTD - A function transfer from Matlab
"""

import numpy as np

def ErrCode(String):
    SERR = '\33[31m'  
    EERR  = '\33[0m'
    print(SERR+String+EERR)

def NoteCode(String):
    SNOTE = '\33[34m'
    ENOTE = '\33[0m'
    print(SNOTE+String+ENOTE)

def Phase_Offset(theta):
    return(theta*np.pi)

### Dispersion Factor Function   
def Dispersion_Factor(PHI,w,w_0): 
    """
    This function calculates the dispersion coefficient Phi used for applying dispersion in frequency space.
    
    Ew*exp(-1i.*Phi), in which Phi = phi_0+phi_1+phi_2+phi_3+...+phi_n
    
    phi_0 is Carrier Envelope Phase - Values of pi() 
    phi_1 is the Group Delay 
    phi_2 is the Group Velocity Dispersion GVD
    phi_3 is the Third Order Dispersion 
    
    Usage:
        Give N number of PHI = [phi_0,phi_1,phi_2,phi_3,...,phi_N]
        Give an Mx1/1xM numpy array of w = np.array([data])
        Give the centre frequency w_0
    """
    if type(w) == list:
        w = np.array(w)
        NoteCode('Converted w from list to array - try to use np.array in the future to speed up your computation')
    try:
        if len(PHI) <= 3:
            NoteCode('Running with suboptimal parameters - Your PHI only contains coefficients up to the GVD in phi_2 - consider adding a TOD term phi_3')
        
        #Actually the part that does it - it's pretty neat if you ask me (it used to be 5 lines of code)
        Phiw = []
        for n in range(len(PHI)):
            Phiw.append(PHI[n]*np.power((w-w_0),n)/np.math.factorial(n))

        return(sum(Phiw))
    except TypeError:
        if type(PHI) == int:   
            ErrCode('TypeError: PHI is type int, make sure that PHI = np.array(phi0,phi1,phi2,phi3)')
        else:
           raise
        
def TDGE(A_t,t,t_0,T,w_0x,theta):
    """
    Time Domain Gaussian Envelope
    """
    return(np.sqrt( A_t*np.exp(-np.log(2)*((2*(t-t_0))/(T))**2))*np.exp(-1.j*(w_0x*(t-t_0)+theta)))

def FDGE(A_w,w,W,w_0x,psi):
    """
    Frequency Domain Gaussian Envelope
    """
    return(np.sqrt(A_w*np.exp(-np.log(2)*((2*(w-w_0x))/(W))**2))*np.exp(-1.j*psi))


        
### FFTD Functions
"""
This function is written to handle transforms to, and application of 
# dispersion between time and frequency domain.
#  
#
#----------------------------------USAGE----------------------------------#
#                     self = FFTD.xxx(x,y,w_0,'type',Phi,Theta)
# x is your x axis (time or frequency axis, select a proper type to match)
#
# y is your intensity/amplitude 
#
# 'type' is a string that indicates which transform you want to do:
#                        |  xxx = ttf = time to frequency   |
#                        |  xxx = ttt = time to time        |
#                        |  xxx = ftt = frequency to time   |
#
#---------------------------External Variables----------------------------#
# Theta, is a single value phase shift in the time domain. Set to 0 usually
#
# Phi is a 4 length vector containing phi_0, phi_1, phi_2, phi_3. You can
# define them as (x,y[phi_0,phi_1,phi_2,phi_3])
#
# phi_0 changes pulse phase, phi_1 changes pulse time "location", 
# phi_2 changes pulse width, phi_3 alters pulse shape in time 
#
#
#
#------------------------------Output Format------------------------------#
#
# All variables will be in A, as per: [A] = FFTD(x,y,w_0,'type',Phi,Theta)
# the actual output will be given as:
# self.variable - e.g. self.Et, gives the non dispersed time axis in the
# type = ttt, time to time dispersion. For the dispersion, you must write
# self.ttt.Etdisp. The same logic goes for all variables:
# -----+-------------------+--------------------+---------------------#
#  ttf | self.Eti & self.t | self.Ewf & self.wf |                     | 
# -----+-------------------+--------------------+---------------------+
#  ttt | self.Eti & self.t | self.Ewf & self.wf | self.Etf & self.tf  |
# -----+-------------------+--------------------+---------------------+
#  ftt |                   | self.Ewi & self.Aw | self.Etf & self.tf  |
# -----+-------------------+--------------------+---------------------+
# Note: using ftt assumes you've already applied dispersion. Inputting Phi
# is just to print PHITEXT, and keep track of which dispersion you used!
"""
def ttf(self,t,I,w_0,Phi,Tht):
     # Setting up variables for storage and calculation 
     N = len(t)                      # Length of series
     self.t    = t                   # Time Axis
     self.t_0  = np.mean(self.t)     # Centre Time Axis (mean of axis)
     self.Eti = I/max(I)             # Normalise E(t) to a(t)
     t1        = self.t[0]           # Start Time
     t2        = self.t[-1]          # End Time 
     Dt        = (t2-t1)             # Sampling interval
     dt        = Dt/N                # Smallest Time interval
     w_t1      = 2*np.pi/Dt           # Min Frequency from time  
     w_t2      = 2*np.pi/dt           # Max Frequency from time
     self.phi  = Phi                                     # Store phi_n terms
     self.phiw = Dispersion_Factor(Phi,self.w,self.w_0)  # Get the dispersion factor 
     self.theta = Phase_Offset(Tht)                     # Get the phase offset
    
     ## Calculate post dispersion delay, and how many units to shift the pulse by to center it 
     #tD    = self.phi[1]                                             # The added time delay by phi_2
     #tcirc = round(N*(round(self.phi[1]/self.Dt)-self.phi[1]/Dt)) 
    
     ## ----- E_t -> E_w ----- ##
     # Perform fft and add dispersion
     self.Ewf      = np.fft.fft(self.Et)                   # fft with no padding
     self.wf       = np.linspace(w_t2,w_t1,N) # Defining the frequency axis - it's "backwards" because the fft goes high to low freq
     self.Ewf      = self.Ewf*np.exp(-1.j*self.phiw) #Apply dispersion
     # Save text of which dispersion factor was used, useful for legend entries
     self.PHITEXT     =         ('\\phi_0 = '+str(self.phi[0])+'\n'+
                                 '\\phi_1 = '+str(self.phi[1])+'\n'+
                                 '\\phi_2 = '+str(self.phi[2])+'\n'+
                                 '\\phi_3 = '+str(self.phi[3])+'\n')
    
def ftt(self,w,A,w_0,Phi,Tht):
     # Setting up variables for storage and calculation 
     N = len(w)                      # Length of series
     self.w    = w                   # Frequency Axis (iff type = ftt)
     self.w_0  = w_0                 # Centre Frequency
     self.Ew     = A/max(A) #Note: The current Ew equation includes the phi component, so you must change that to apply a different dispersion.
     w1 = self.w[0]                  # Highest Freq
     #w2 = self.w[-1]                 # Lowest Freq
     #Dw   = w2 - w1                  # Frequency Difference
     #dw   = Dw/N                     # Smallest Frequency inverval
     Dt_w = 2*np.pi/w1               # Determining Max time interval from frequency
     #dt_w = 2*np.pi/w2               # Determining Min time interval from frequency
     t_w1 = -Dt_w/2                  # Determining time axis start    
     t_w2 =  Dt_w/2                  # Determining time axis end variables
     self.phi  = Phi                                     # Store phi_n terms
     self.phiw = Dispersion_Factor(Phi,self.w,self.w_0)  # Get the dispersion factor 
     self.theta = Phase_Offset(Tht)                     # Get the phase offset
    
#     ## Calculate post dispersion delay, and how many units to shift the pulse by to center it 
#     tD    = self.phi[1]                                             # The added time delay by phi_2
#     tcirc = round(N*(round(self.phi[1]/Dt)-self.phi[1]/Dt)) 
    
     ## ---- E_w -> E_t ---- ##
     # Perform ifft 
     self.Etf  = np.fft.fftshift(np.fft.ifft(self.Ew))           # since we start with frequency, we need fftshift to get dispersed Et
     self.tf   = np.linspace(t_w1,t_w2,N) # Get time axis
   
     #In this function, inputting phi is only required to keep track of
     #which was used. The code treats the ftt example has ALREADY HAVING
     #DISPERSION - You can modify this code easily for yourself, but this
     #was done to keep it simple!
     self.PHITEXT     =          ('\\phi_0 = '+str(self.phi[0])+'\n'+
                                  '\\phi_1 = '+str(self.phi[1])+'\n'+
                                  '\\phi_2 = '+str(self.phi[2])+'\n'+
                                  '\\phi_3 = '+str(self.phi[3])+'\n')

def ttt(self,t,I,w_0,Phi,Tht):
     # Setting up variables for storage and calculation 
     N = len(t)                      # Length of series
     self.t     = t                   # Time Axis
     self.t_0   = np.mean(self.t)     # Centre Time Axis (mean of axis)
     self.Et = I/max(I)              # Normalise E(t) to a(t)
     t1         = self.t[0]           # Start Time
     t2         = self.t[-1]          # End Time 
     Dt         = (t2-t1)             # Sampling interval
     dt         = Dt/N                # Smallest Time interval
     w_t1       = 2*np.pi/Dt           # Min Frequency from time  
     w_t2       = 2*np.pi/dt           # Max Frequency from time
     self.phi   = Phi                                     # Store phi_n terms
     self.phiw  = Dispersion_Factor(Phi,self.w,self.w_0)  # Get    the dispersion factor 
     self.theta = Phase_Offset(Tht)                      # Get the phase offset
    
     ## Calculate post dispersion delay, and how many units to shift the pulse by to center it 
     tD         = self.phi[1]                                # The added time delay by phi_2
     tcirc      = round(N*(round(self.phi[1]/self.Dt)-self.phi[1]/Dt)) 
    
     ## ----- E_t -> E_w ----- ##
     # Perform fft and add dispersion
     self.Etf   = np.fft.fft(self.Et)          # fft with no padding
     self.wf    = np.linspace(w_t2,w_t1,N)     # Defining the frequency axis - it's "backwards" because the fft goes high to low freq
     #self.Etf   = self.Etf*exp(-1.j*self.phiw) #Apply dispersion
    
     ## ---- E_w -> E_t ---- ##
     self.Etf   = np.flip((np.fft.ifft(self.Etf)))     #We flip Ew before ifft, because otherwise it is backwards when transformed back
     Dt = 2*np.pi/self.wf[-1]                  #Extract max time
     #dt_w = 2*np.pi/self.wf[1];                #Extract min time    
     t1 =   -self.t_0                          # set t1
     t2 = Dt-self.t_0;                         # set t2
     self.tf = np.linspace(t1,t2,self.N);        # Get time axis (idea case => same as input)
     self.Etfc = np.roll(self.Etf,tcirc);       # Central pulse
     self.tfc  = self.tf + tD;          # Time shifted axis - for centralised pulse
     # Save text of which dispersion factor was used, useful for legend entries
     self.PHITEXT     =         ('\\phi_0 = '+str(self.phi[0])+'\n'+
                                 '\\phi_1 = '+str(self.phi[1])+'\n'+
                                 '\\phi_2 = '+str(self.phi[2])+'\n'+
                                 '\\phi_3 = '+str(self.phi[3])+'\n')

        
def RefrInd(LAMBD,MATERIAL):
    if MATERIAL == 'SiO2':
        RefInd = np.sqrt( ((0.6961663*(LAMBD*10**6)**2)/((LAMBD*10**6)**2 - 0.0684043**2)) +
                           ((0.4079426*(LAMBD*10**6)**2)/((LAMBD*10**6)**2 - 0.1162414**2)) +
                           ((0.8974794*(LAMBD*10**6)**2)/((LAMBD*10**6)**2 - 9.896161**2))  +  1)
    return(RefInd);
                       

%Take Derivatives of the function defining n, we need them for \Phi
syms WavLen
fun_nd1 = eval(['@(WavLen)' char(diff(fun_n(WavLen)))]);
nd1     = fun_nd1(lambda_0)*(10^-6);
fun_nd2 = eval(['@(WavLen)' char(diff(fun_nd1(WavLen)))]);
nd2     = fun_nd2(lambda_0)*(10^-6);
fun_nd3 = eval(['@(WavLen)' char(diff(fun_nd2(WavLen)))]);
nd3     = fun_nd3(lambda_0)*(10^-6);
clear WavLen
                       
%Graph of equation for the relevant range above - we don't really need this though
nind = n_ref; lambrange = linspace(0.21,6.7,1000); % \mu m

%Calculate Group dispersion, Group Velocity Dispersion and Third order Dispersion;
GD    = 1/((c/n_ref)/(1 - (lambda_0/n_ref) * nd1)) * 1000;                                 %units s/mm
GVD   = (lambda_0^3)/(2 * pi * c^2) * nd2 * 1000;                                          %units s^2/mm
TOD   = -((lambda_0)/(2 * pi() * c))^2 * 1/c * (3*lambda_0^2*nd2+lambda_0^3*nd3) * 1000;   %units s^3/mm

%Calculate terms dispersion factor
SiO2.n   = n_ref;        % Refractive Index  - CEP replaces this later. 
SiO2.CD  = L*(0.025+nd1)*10^6;  %                   - Cromatic dispersion
SiO2.GD  = L*GD;         % s^2               - Group Delay
SiO2.GVD = L*GVD;        % s^2               - Group Velocity dispersion
SiO2.TOD = L*TOD;        % s^3               - Third Order Dispersion

ORD     = 5; 
BPF     = 0;
        