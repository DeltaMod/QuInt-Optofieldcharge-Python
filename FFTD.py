# -*- coding: utf-8 -*-
"""
FFTD - A function transfer from Matlab
"""

from scipy import signal
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

#%% Dispersion Factor Function   
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
            
        Phiw = []
        for n in range(len(PHI)):
            Phiw.append(PHI[n]*np.power((w-w_0),n)/np.math.factorial(n))

        return(sum(Phiw))
    except TypeError:
        if type(PHI) == int:   
            ErrCode('TypeError: PHI is type int, make sure that PHI = np.array(phi0,phi1,phi2,phi3)')
        else:
           raise
        
        
        
#%% FFTD Functions
"""
This function is written to handle transforms to, and application of 
% dispersion between time and frequency domain.
%  
%
%----------------------------------USAGE----------------------------------%
%                     A = FFTD(x,y,w_0,'type',Phi,Theta)
% x is your x axis (time or frequency axis, select a proper type to match)
%
% y is your intensity/amplitude 
%
% 'type' is a string that indicates which transform you want to do:
%                        |  ttf = time to frequency   |
%                        |  ttt = time to time        |
%                        |  ftt = frequency to time   |
%
%---------------------------External Variables----------------------------%
% Theta, is a single value phase shift in the time domain. Set to 0 usually
%
% Phi is a 4 length vector containing phi_0, phi_1, phi_2, phi_3. You can
% define them as (x,y[phi_0,phi_1,phi_2,phi_3])
%
% phi_0 changes pulse phase, phi_1 changes pulse time "location", 
% phi_2 changes pulse width, phi_3 alters pulse shape in time 
%
%
%
%------------------------------Output Format------------------------------%
%
% All variables will be in A, as per: [A] = FFTD(x,y,w_0,'type',Phi,Theta)
% the actual output will be given as:
% A.type.variable - e.g. A.ttt.Et, gives the non dispersed time axis in the
% type = ttt, time to time dispersion. For the dispersion, you must write
% A.ttt.Etdisp. The same logic goes for all variables:
%  -----------------------------------------------------------------------%
%  | A.ttf.Et & A.ttf.t | A.ttf.Ew & A.ttf.w |                           | 
%  |--------------------|--------------------|---------------------------|
%  | A.ttt.Et & A.ttt.t | A.ttt.Ew & A.ttt.w | A.ttt.Etdisp & A.ttt.tdisp|
%  |--------------------|--------------------|---------------------------|
%  |                    | A.ftt.Ew & A.ftt.w | A.ftt.Etdisp & A.ftt.tdisp|
%  -----------------------------------------------------------------------
% Note: using ftt assumes you've already applied dispersion. Inputting Phi
% is just to print PHITEXT, and keep track of which dispersion you used!
"""
def FFTDINIT(self,x,y,Phi,Tht,selector):
    #Temporal Variables
    self.N    = length(x)                                       # N length vector  
    if selector == 'Temporal':                                                           
        self.t    = x                                           # Time Axis
        self.t1   = self.t[0]; self.t2 = self.t[-1]             # Start and end time
        self.t_0  = np.mean(self.t)                             # Centre Time Axis (mean of axis)
        self.Dt   = (t2-t1)                                     # Sampling interval
        self.dt   = self.Dt/self.N                              # Smallest Time interval
        self.w_t1 = 2*pi()/self.Dt; self.w_t2 = 2*pi()/self.dt  # Min-Max Frequency from time graph
    
    #Frequency Variables
    if selector == 'Frequency'
        self.w    = x;                                          # Frequency Axis (iff type = ftt)
        self.w_0  = w_0                                         # Centre Frequency
        self.w1 = self.w[0]; self.w2 = self.w[-1];              # Highest/Lowest Freq
        self.Dw   = self.w2 - self.w1                           # Frequency Difference
        self.dw   = self.Dw/self.N                              # Smallest Frequency inverval
        self.Dt_w = 2*pi()/self.w1; self.dt_w = 2*pi()/self.w2  # Determining Max/min time interval from frequency
        self.t_w1 = -self.Dt_w/2; self.t_w2 = self.Dt_w/2       # Determining time axis start/end variables
    
    
    fun_theta = @(theta) theta*pi();                                        # This sets phase offset
    fun_phiw = @(phi0, phi1, phi2, phi3,w,w_0)     phi0                   +...  #phi_0 is Carrier Envelope Phase - Values of pi() 
                                                   phi1.*((w-w_0))        +...  #phi_1 is the Group Delay 
                                                   phi2.*((w-w_0).^2)/2   +...  #phi_2 is the Group Velocity Dispersion GVD
                                                   phi3.*((w-w_0).^3)/6 ;       #phi_3 is the Third Order Dispersion 
    A.phi = Phi;                                                            # Get \phi values in a struct 
    A.phiw = fun_phiw(A.phi(1), A.phi(2), A.phi(3), A.phi(4),A.w,A.w_0);    # Get dispersion factor 
    A.theta = fun_theta(Theta);                                             # Get phase offset
    
    ## Calculate post dispersion delay, and how many units to shift the pulse by to center it 
    A.tD    = A.phi(2); % The added time delay by phi_2
    A.tcirc = round(A.N*(round(A.phi(2)/A.Dt)-A.phi(2)/A.Dt)); 

def ttf(t,A,w_0,Phi,Tht):
    
def ftt(w,A,w_0,Phi,Tht):
    
def ttt(t,A,w_0,Phi,Tht):
    
        

        