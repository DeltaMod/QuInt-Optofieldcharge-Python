# -*- coding: utf-8 -*-
'''
Created on Fri Nov  8 13:31:20 2019

@author: Vidar Flodgren
'''

# -*- coding: utf-8 -*-
'''
FFAST-MPEG - A quicker than usual way to make gifs, trim videos, split videos and so on!
'''
import json
import OFunc
from OFunc import FFTD as FFTD
import KGD2UI
import os, sys
import tkinter as tk
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import (
    FigureCanvasTkAgg, NavigationToolbar2Tk)
from matplotlib.figure import Figure
import numpy as np
import sympy as sp
from scipy.constants import c
from scipy.constants import epsilon_0 as eps_0
sp.init_printing(use_latex=False) 
#import matplotlib
#matplotlib.use('TKagg') 
#Define Colours to use for certain components
plt.rcParams['figure.dpi']         = 150
plt.rcParams['axes.grid']          = True
plt.rcParams['axes.axisbelow']     = True
plt.rcParams['axes.spines.left']   = True
plt.rcParams['axes.spines.right']  = True
plt.rcParams['axes.spines.bottom'] = True
plt.rcParams['axes.spines.top']    = True
plt.rcParams['axes.spines.top']    = True
plt.rcParams['ytick.left']         = True
plt.rcParams['ytick.labelleft']    = True
plt.rcParams['xtick.bottom']       = True
plt.rcParams['xtick.labelbottom']  = True

rootBG   = 'aliceblue'
PlotBG   = 'ghostwhite'
ButtonBG = 'WhiteSmoke'
ButtonABG = 'azure'
LabelBG  = 'aliceblue'
VarBG    = LabelBG 

SIMSEL = ['Simple Photoinduced Charge',
                    'Dispersion and Photoinduced Charge',
                    'Dispersion Grapher'];

PLTSEL =      [['Gaussian Laser Pulse',        
               '<a^2n+1>',                     
               'Photoinduced Charge' ,         
               '<a^2n+1> Term Contributions to Charge']];
PLTSEL.append(['Non Fourier Et(t)'           ,
               'E_t FFT Plot'                 ,
               'E_t Final Applied Dispersion' ,
               'E_t(t)^2n+1 After Dispersion' ,
               'Post Dispersion Photoinduced Charge',
               'E_w IFFT Dispersion'          ,
               'Current Pulse Overlap at given Delta t',
               'Temporal Overlap Induced Charge',
               'Dual Polarised Induced Charge'])
PLTSEL.append(['No Phi'                       ,
               'Phi0'                       ,  
               'Phi1'                        , 
               'Phi2'                        , 
               'Phi3'                        ,
               'All_Phi'])
# Make the lists readable through strings - This is done to reduce the number 
#of compute commands to "find" the correct index using an OptionsBox
SIMSEL = [str(n)+' - '+SIMSEL[n] for n in range(len(SIMSEL))]
             
for i in range(len(PLTSEL)):
    for j in range(len(PLTSEL[i])):
        PLTSEL[i][j] = str(j)+' - '+PLTSEL[i][j]

PLTSELLen = len(max(PLTSEL, key=len))
SIMSELLen = len(max(SIMSEL, key=len))
SelBoxLen = max([PLTSELLen,SIMSELLen])

VAR = {'Sim Select':0,'Plot Select':0,
            'f_0x':0,'f_0y':0, 'Ncycx':0, 'Ncycy':0, 'F_x':0, 'F_y':0, 't1':0, 't2':0, 'L':0, 'Aeff':0,
            'BPF':0,  'TDn':0, 'Delt1':0, 'Delt2':0,  'ORD':0,  'N':100,'Mat':'SiO2'}



        
class KGUI:
    LABEL_TEXT = ['FFAST-MPEG']
    def __init__(self, master):
        self.master = master 
        master.title('Khurgin Py-UI')
        XDIM = 1280
        YDIM = 720
        root.geometry('{}x{}'.format(XDIM,YDIM))
        root.config(bg='aliceblue')
        tk.Grid.rowconfigure(root, 0, weight=1)
        tk.Grid.columnconfigure(root, 0, weight=1)
        tk.Grid.columnconfigure(root, 2, weight=1)
        tk.Grid.rowconfigure(root, 1, weight=1)
        tk.Grid.rowconfigure(root, 2, weight=1)
        self.SIMSEL = SIMSEL
        self.PLTSEL = PLTSEL
        
        self.SavedFields = ['Sim Select','Plot Select','f_0x','f_0y', 'Ncycx', 'Ncycy', 'F_x', 'F_y',
            't1', 't2', 'L', 'Aeff', 'BPF',  'TDn', 'Delt1', 'Delt2',  'ORD',  'N','Mat']
        
        self.EntryList = ['f_0x','f_0y', 'Ncycx', 'Ncycy', 'F_x', 'F_y',
            't1', 't2', 'L', 'Aeff', 'BPF',  'TDn', 'Delt1', 'Delt2',  'ORD',  'N']
        if 'SessionRestore.dat' in os.listdir():
            with open("SessionRestore.dat") as SessionRestore:
                pyresponse = json.loads(SessionRestore.read())
                print(pyresponse)
                self.VAR = pyresponse 
            for entry in self.SavedFields:
                if entry not in list(self.VAR.keys()):
                    if entry=='Mat':
                        self.VAR[entry] = 'SiO2'
                    else:
                        self.VAR[entry] = '0'
                    
        else:
            self.VAR = {}
            for entry in self.SavedFields:
                self.VAR[entry] = '0'  
            
          
        #Creating Graph Frame
        self.GraphPlot = plt.figure(1)
        self.GraphPlot.subplots_adjust(bottom=None, top=None, left=None, right=None)
        self.GCanv = FigureCanvasTkAgg(self.GraphPlot,master=root)  # A tk.DrawingArea. #,bg='white',width=YDIM, height=YDIM, relief = 'raised'
        self.GCanv.draw()
        self.GCanv.get_tk_widget().grid(column=0,row=0,rowspan=4,sticky='nwes',padx=10,pady=10)
        
        self.TBFrame = tk.Frame(root,bg='ghostwhite',)
        self.TBFrame.grid(row=0,column=0,sticky='nw',padx=10,pady=10)
        self.GToolbar = NavigationToolbar2Tk(self.GCanv, self.TBFrame)
        self.GToolbar.update()
       
        #Creating Preset Box Frame
        self.PSFrame = tk.Frame(root,bg='ghostwhite')
        self.PSFrame.grid(row = 0, column = 1,sticky='nswe')
        
        
        #Creating Editbox Frame
        self.EBFrame = tk.Frame(root,bg='ghostwhite')
        self.EBFrame.grid(row = 1, column = 1,sticky = 'nwe')
        self.EBFrame.bind('<Configure>',self.MaintainAspect)
        
        
        #Creating Simulation Box Frame
        self.RSFrame = tk.Frame(root,bg='ghostwhite')
        self.RSFrame.grid(row = 2, column = 1,sticky='nwe')
        #Setting Weights of all root relevant columms to be 1
       
                 

        
         #%% Listing all of the location and span variables used to cotrol the appearance of the UI
         #In the future, consider instead creating a matrix that uses exec to automatically generate these boxes based upon some matrix storage, like:
         #A = [Name,Type,Frame,text,relief,anchor,bg,row,column,columnspan,sticky]
        #Preset Frame
        #Variable Frame
        ROWf_0x   = 4;  COLf_0x   = 2; DEntSpan = 3
        ROWf_0y   = 4;  COLf_0y   = 6; SEntSpan = 5
        ROWNcycx  = 6;  COLNcycx  = 2; 
        ROWNcycy  = 6;  COLNcycy  = 7; 
        ROWF_x   = 8;  COLF_x   = 2; DDISpan  = 1
        ROWF_y   = 8;  COLF_y   = 7; LabSpan  = 1
        ROWt1     = 10; COLt1     = 2; EWDT = 8
        ROWt2     = 10; COLt2     = 7; 
        ROWL      = 12; COLL      = 3; 
        ROWAeff   = 14; COLAeff   = 3; 
        ROWBPF    = 16; COLBPF    = 3;
        ROWORD    = 18; COLORD    = 3;
        ROWN      = 20; COLN      = 3;
        
        #Slider Frame
        ROWTDn    = 22; COLTDn    = 2;
        ROWDelt1  = 24; COLDelt1  = 1;
        ROWDelt2  = 24; COLDelt2  = 8;
        ROWSLD    = 24; COLSLD    = 3;
        #Sim Select Frame
        ROWSIMSEL = 0; COLSIMSEL = 4; 
        ROWPLTSEL = 2; COLPLTSEL = 4;
        ROWCONVRT = 4; COLCONVRT = 0;
        ROWCLOSE  = 4; COLCLOSE  = 10;
        ROWRESET  = 4; COLRESET  = 5;
        
        #Column for unit labels
        COLUNITS = 11;
        
        
        #root.update()
        #print(self.EBFrame.winfo_width())
        for i in range(0,ROWSLD,2):
            self.EBFrame.grid_rowconfigure(i, weight=1,uniform='fred')
        for j in range(0,9):
            if j not in [0,1,4,6]:
                self.EBFrame.grid_columnconfigure(j, weight=1,uniform='fred')

        
        #%% Preset Frame
        ##Run Simulation Button
        self.Convert_Button = tk.Button(self.PSFrame, text='Popout Graph', command=self.RunSim)
        self.Convert_Button.grid(row = ROWCONVRT,column = COLCONVRT,rowspan=2,columnspan=2,sticky='nwse')
        
        #Close Button
        self.close_button = tk.Button(self.PSFrame, text='Load Preset', command=self.close)
        self.close_button.grid(row = ROWCLOSE,column = COLCLOSE,sticky='w')
        
        #Reset Button
        self.save_preset_button = tk.Button(self.PSFrame, text='Save Preset', command=self.Save_Preset)
        self.save_preset_button.grid(row = ROWRESET,column = COLRESET,sticky='w')
        
        #%% Double Entries
        
        #Driving and Injection Labels
        tk.Label(self.EBFrame, text="Driving",relief='flat',justify='center',anchor='n',bg=LabelBG).grid(row=ROWf_0x-1,column=COLf_0x,columnspan=DEntSpan,sticky='we')
        tk.Label(self.EBFrame, text="Injection",relief='flat',anchor='n',bg=LabelBG).grid(row=ROWf_0y-1,column=COLf_0y,columnspan=DEntSpan,sticky='we')
        #x and y labels, just to make it absolutely clear:
        dROW = abs(ROWf_0x - ROWNcycx)
        xlbl = ['x','x','x','t1']
        ylbl = ['y','y','y','t2']
        mlbl = ['  ','  ','  ','  ']
        for i in range(4):
            tk.Label(self.EBFrame, text=xlbl[i],relief='flat',anchor='e',bg=LabelBG).grid(row=dROW*i+ROWf_0x,column=COLf_0x-1,columnspan=DDISpan,sticky='e')
            tk.Label(self.EBFrame, text=ylbl[i],relief='flat',anchor='e',bg=LabelBG).grid(row=dROW*i+ROWf_0y,column=COLf_0y-1,columnspan=DDISpan,sticky='e') 
            tk.Label(self.EBFrame, text=mlbl[i],relief='flat',anchor='w',bg=LabelBG).grid(row=dROW*i+ROWf_0x,column=COLL+2,columnspan=1,sticky='w')
        Rowunits = ['THz','No Units','V^Åm^-1','fs','mm','10^-12m^2','max(A_w)/BFP','No Units','No Units','fs','No Units']
        for i in range(len(Rowunits)):
            tk.Label(tk.Label(self.EBFrame, text=Rowunits[i],relief='flat',anchor='e',bg=LabelBG).grid(row=dROW*i+ROWf_0x,column=COLUNITS,columnspan=DDISpan,sticky='w'))
        
        self.vcmd = (root.register(self.validate),'%P')
        """
        %d Type of action: 1 for insert, 0 for delete, or -1 for focus, forced or textvariable validation.
        
        %i Index of char string to be inserted/deleted, if any, otherwise -1.
        
        %P The value of the entry if the edit is allowed. If you are configuring the entry widget to have a new textvariable, this will be the value of that textvariable.
        
        %s The current value of entry prior to editing.
        
        %S The text string being inserted/deleted, if any, {} otherwise.
        
        %v The type of validation currently set.
        
        %V The type of validation that triggered the callback (key, focusin, focusout, forced).
        
        %W The name of the entry widget.
        """
        self.VAR_SimulationCompleted   = tk.BooleanVar(root); self.VAR_SimulationCompleted.set(False)
        tk.Label(self.EBFrame, text="Laser Frequencies (f_0)",relief='flat',anchor='w',bg=LabelBG).grid(row=ROWf_0x,column=0,columnspan=LabSpan,sticky='we')
        
        self.VAR_f_0x   = tk.StringVar(root); self.VAR_f_0x.set(str(self.VAR['f_0x']))
        self.EBf_0x = tk.Entry(self.EBFrame,validate='key',validatecommand=(self.vcmd,'%P'),justify='center',width=EWDT,textvariable=self.VAR_f_0x)
        self.EBf_0x.grid(row=ROWf_0x,column=COLf_0x,columnspan=DEntSpan,sticky='we')
        
        self.VAR_f_0y   = tk.StringVar(root); self.VAR_f_0y.set(str(self.VAR['f_0y']))
        self.EBf_0y = tk.Entry(self.EBFrame,justify='center',validate='key',validatecommand=self.vcmd,width=EWDT,textvariable=self.VAR_f_0y)
        self.EBf_0y.grid(row=ROWf_0y,column=COLf_0y,columnspan=DEntSpan,sticky='we')
        
        tk.Label(self.EBFrame, text="Optical Cycles (f_0)",relief='flat',anchor='w',bg=LabelBG).grid(row=ROWNcycx,column=0 ,columnspan=LabSpan,sticky='we') 
        self.VAR_Ncycx   = tk.StringVar(root); self.VAR_Ncycx.set(str(self.VAR['Ncycx']))
        self.EBNcycx = tk.Entry(self.EBFrame,justify='center',validate='key',validatecommand=self.vcmd,width=EWDT,textvariable=self.VAR_Ncycx)
        self.EBNcycx.grid(row=ROWNcycx,column=COLNcycx,columnspan=DEntSpan,sticky='we')
        
        self.VAR_Ncycy   = tk.StringVar(root); self.VAR_Ncycy.set(str(self.VAR['Ncycy']))
        self.EBNcycy = tk.Entry(self.EBFrame,justify='center',validate='key',validatecommand=self.vcmd,width=EWDT,textvariable=self.VAR_Ncycy)
        self.EBNcycy.grid(row=ROWNcycy,column=COLNcycy,columnspan=DEntSpan,sticky='we')
        
        tk.Label(self.EBFrame, text="Optical Field Strength (F_)",relief='flat',anchor='w',bg=LabelBG).grid(row=ROWF_x,column=0 ,columnspan=LabSpan,sticky='we') 
        self.VAR_F_x   = tk.StringVar(root); self.VAR_F_x.set(str(self.VAR['F_x']))
        self.EBF_x = tk.Entry(self.EBFrame,justify='center',validate='key',validatecommand=self.vcmd,width=EWDT,textvariable=self.VAR_F_x)
        self.EBF_x.grid(row=ROWF_x,column=COLF_x,columnspan=DEntSpan,sticky='we')
        
        self.VAR_F_y   = tk.StringVar(root); self.VAR_F_y.set(str(self.VAR['F_y']))
        self.EBF_y = tk.Entry(self.EBFrame,justify='center',validate='key',validatecommand=self.vcmd,width=EWDT,textvariable=self.VAR_F_y)
        self.EBF_y.grid(row=ROWF_y,column=COLF_y,columnspan=DEntSpan,sticky='we')
        
        tk.Label(self.EBFrame, text="Time Range (t)",relief='flat',anchor='w',bg=LabelBG).grid(row=ROWt1,column=0,columnspan=LabSpan,sticky='we')
        self.VAR_t1   = tk.StringVar(root); self.VAR_t1.set(str(self.VAR['t1']))
        self.EBt1 = tk.Entry(self.EBFrame,justify='center',validate='key',validatecommand=self.vcmd,width=EWDT,textvariable=self.VAR_t1)
        self.EBt1.grid(row=ROWt1,column=COLt1,columnspan=DEntSpan,sticky='we')
        
        self.VAR_t2   = tk.StringVar(root); self.VAR_t2.set(str(self.VAR['t2']))
        self.EBt2 = tk.Entry(self.EBFrame,justify='center',validate='key',validatecommand=self.vcmd,width=EWDT,textvariable=self.VAR_t2)
        self.EBt2.grid(row=ROWt2,column=COLt2,columnspan=DEntSpan,sticky='we')
        
        #%% Single Entries
        tk.Label(self.EBFrame, text="Material Thickness (L)",relief='flat',anchor='w',bg=LabelBG).grid(row=ROWL,column=0,columnspan=LabSpan,sticky='we')
        self.VAR_L   = tk.StringVar(root); self.VAR_L.set(str(self.VAR['L']))
        self.EBL = tk.Entry(self.EBFrame,justify='center',validate='key',validatecommand=self.vcmd,width=EWDT,textvariable=self.VAR_L)
        self.EBL.grid(row=ROWL,column=COLL,columnspan=SEntSpan,sticky='we')
        
        tk.Label(self.EBFrame, text="Effective Area (Aeff)",relief='flat',anchor='w',bg=LabelBG).grid(row=ROWAeff,column=0,columnspan=LabSpan,sticky='we')
        self.VAR_Aeff   = tk.StringVar(root); self.VAR_Aeff.set(str(self.VAR['Aeff']))
        self.EBAeff = tk.Entry(self.EBFrame,justify='center',validate='key',validatecommand=self.vcmd,width=EWDT,textvariable=self.VAR_Aeff)
        self.EBAeff.grid(row=ROWAeff,column=COLAeff,columnspan=SEntSpan,sticky='we')
        
        tk.Label(self.EBFrame, text="Band Pass Filter",relief='flat',anchor='w',bg=LabelBG).grid(row=ROWBPF,column=0,columnspan=LabSpan,sticky='we')
        self.VAR_BPF = tk.StringVar(root); self.VAR_BPF.set(str(self.VAR['BPF']))
        self.EBBPF = tk.Entry(self.EBFrame,justify='center',validate='key',validatecommand=self.vcmd,width=EWDT,textvariable=self.VAR_BPF)
        self.EBBPF.grid(row=ROWBPF,column=COLBPF,columnspan=SEntSpan,sticky='we')
        
        tk.Label(self.EBFrame, text="Orders of <a^2n+1> (ORD)",relief='flat',anchor='w',bg=LabelBG).grid(row=ROWORD,column=0,sticky='we')
        self.VAR_ORD = tk.StringVar(root); self.VAR_ORD.set(str(self.VAR['ORD']))
        self.EBORD = tk.Entry(self.EBFrame,justify='center',validate='key',validatecommand=self.vcmd,width=EWDT,textvariable=self.VAR_ORD)
        self.EBORD.grid(row=ROWORD,column=COLORD,columnspan=SEntSpan,sticky='we')
        
        tk.Label(self.EBFrame, text="Sampling Points (N)",relief='flat',anchor='w',bg=LabelBG).grid(row=ROWN,column=0,sticky='w')
        self.VAR_N = tk.StringVar(root); self.VAR_N.set(str(self.VAR['N']))
        self.EBN = tk.Entry(self.EBFrame,justify='center',validate='key',validatecommand=self.vcmd,width=EWDT,textvariable=self.VAR_N)
        self.EBN.grid(row=ROWN,column=COLN,columnspan=SEntSpan,sticky='we')
        
         #Temporal Delay
        self.VAR_TBTDn = tk.StringVar(root); self.VAR_TBTDn.set(str(self.VAR['TDn']))
        
        tk.Label(self.EBFrame, text="Teporal Delay Δt                          ",relief='flat',anchor='w',bg=LabelBG).grid(row=ROWTDn,column=0,columnspan=LabSpan,sticky='we')
        tk.Label(self.EBFrame, text="=",relief='flat',anchor='c',bg=LabelBG).grid(row=ROWTDn,column =COLL+1,columnspan=3,sticky='we')
        tk.Label(self.EBFrame,justify='center',width=EWDT,textvariable=self.VAR_TBTDn).grid(row=ROWTDn,column=COLDelt2-1,columnspan=2,sticky='we')
        
        self.VAR_TDn = tk.StringVar(root); self.VAR_TDn.set(str(self.VAR['TDn']))
        self.EBTDn = tk.Entry(self.EBFrame,justify='center',validate='key',validatecommand=self.vcmd,width=EWDT,textvariable=self.VAR_TDn)
        self.EBTDn.grid(row=ROWTDn,column=COLTDn,columnspan=2,sticky='we')
        #Temporal Delay
        
        tk.Label(self.EBFrame, text="n_min",relief='flat',anchor='e',bg=LabelBG).grid(row=ROWDelt1,column=COLDelt1-1,sticky='swe')
        self.VAR_Delt1 = tk.StringVar(root); self.VAR_Delt1.set(str(self.VAR['Delt1']))
        self.EBDelt1 = tk.Entry(self.EBFrame,justify='center',validate='key',validatecommand=self.vcmd,width=EWDT,textvariable=self.VAR_Delt1)
        self.EBDelt1.grid(row=ROWDelt1,column=COLDelt1,columnspan=2,sticky='wse')
        
        ##Frame Slider 
        self.DeltSlider = tk.Scale(self.EBFrame,from_=-1, to=1, orient='horizontal',width=10,length=100,sliderlength=10)
        self.DeltSlider.grid(row=ROWSLD,column=COLSLD, columnspan=5, sticky='swe')
        
        ##Delt2
        self.VAR_Delt2 = tk.StringVar(root); self.VAR_Delt2.set(str(self.VAR['Delt2']))
        self.EBDelt2 = tk.Entry(self.EBFrame,justify='center',validate='key',validatecommand=self.vcmd,width=EWDT,textvariable=self.VAR_Delt2)
        self.EBDelt2.grid(row=ROWDelt2,column=COLDelt2,columnspan=2,sticky='swe')
        
        
        
        for entry in self.EntryList:
            self.__getattribute__('EB'+entry).bind('<FocusOut>',self.Entry_Focus_Out)
            self.__getattribute__('EB'+entry).bind('<Return>',self.Entry_Focus_Out)          
        #%% Populating the Simulation Selection Frame
        #Defining StringVar for SIMSEL Options Menu
        self.SSV = tk.StringVar(root)
        self.SSV.set(SIMSEL[VAR['Sim Select']])
        
        #Generating Simulation Selector options menu
        self.DDSIMS= tk.OptionMenu(self.RSFrame, self.SSV, *SIMSEL)
        tk.Label(self.RSFrame, text="Select a Simulation",relief='flat',anchor='w',bg=LabelBG).grid(row = ROWSIMSEL, column = 0,sticky='we')
        self.DDSIMS.grid(row = ROWSIMSEL, column =COLSIMSEL,columnspan=10,sticky='w')
        self.DDSIMS.config(width = SelBoxLen,bg=ButtonBG,activebackground =ButtonABG)
        self.DDSIMS["menu"].config(bg=ButtonABG)
        def PLT_DDGen(self):
            #Defining StringVar for PLTSEL Options Menu
            self.PSV = tk.StringVar(root)
            self.PSV.set(PLTSEL[int(self.SSV.get()[0])][int(str(self.VAR['Plot Select'])[0])])
            #Generating Plot Selector options menu
            self.DDPLT= tk.OptionMenu(self.RSFrame, self.PSV, *PLTSEL[int(self.SSV.get()[0])])
            tk.Label(self.RSFrame, text="Select a Plot",relief ='flat',anchor='w',bg=LabelBG ).grid(row = ROWPLTSEL, column = 0,sticky='we')
            self.DDPLT.grid(row = ROWPLTSEL, column =COLPLTSEL,columnspan=5,sticky='w')
            self.DDPLT.config(width = SelBoxLen,bg=ButtonBG,activebackground =ButtonABG)
            self.DDPLT["menu"].config(bg=ButtonABG)
            
        PLT_DDGen(self)   
        
        def SIM_DDChange(*args):
            self.VAR['Sim Select'] = int(self.SSV.get()[0])
            if int(self.SSV.get()[0]) == 0:
                self.VAR_L.set(0)
                self.EBL.config(state='disabled')
            elif int(self.SSV.get()[0]) != 0:
                self.VAR_L.set(self.VAR['L'])
                self.EBL.config(state='normal')
            self.VAR['Plot Select'] = 0
            self.DDPLT.destroy() #Destroy old OptionsMenu so it can be re-created (there might be a better way to do this)
            PLT_DDGen(self)      #Re generate PLT OptionsMenu
        def PLT_DDChange(*args):
            VAR['Plot Select'] = int(self.PSV.get()[0])
            print(self.VAR_SimulationCompleted.get())
            if self.VAR_SimulationCompleted.get() == True:
                print('Is this working?')
                self.Result_Plotter()
            
        # link function to change dropdown
        self.SSV.trace('w',SIM_DDChange)
        self.PSV.trace('w',PLT_DDChange)
        
      
                
        #VAR = {'Sim Select':0,'Plot Select':0,
        #   'f_0x':0,'f_0y':0, 'Ncycx':0, 'Ncycy':0, 'F_x':0, 'F_y':0,
        #   't1':0, 't2':0, 'L':0, 'Aeff':0, 'BPF':0,  'TDn':0, 'Delt1':0, 'Delt2':0,  'ORD':0,  'N':0}
         
        ##Run Simulation Button
        self.Convert_Button = tk.Button(self.RSFrame, text='Run Simulation', command=self.RunSim)
        self.Convert_Button.grid(row = ROWCONVRT,column = COLCONVRT,rowspan=2,columnspan=2,sticky='nwse')
         
        #Close Button
        self.close_button = tk.Button(self.RSFrame, text='Close', command=self.close)
        self.close_button.grid(row = ROWCLOSE,column = COLCLOSE,sticky='w')
        
        #Reset Button
        self.reset_button = tk.Button(self.RSFrame, text='Reset', command=self.reset)
        self.reset_button.grid(row = ROWRESET,column = COLRESET,sticky='w')
    def Entry_Focus_Out(self,event):
        if event.keysym == 'Return':
            root.focus_set()
        for entry in self.EntryList:
            if self.__getattribute__('VAR_'+entry).get()=='':
                self.__getattribute__('VAR_'+entry).set('0')
            try: 
                self.VAR[entry] = eval(self.__getattribute__('VAR_'+entry).get())
                self.__getattribute__('VAR_'+entry).set(str(self.VAR[entry]))
            except:
                self.__getattribute__('VAR_'+entry).set(str(self.VAR[entry]))
            
            self.VAR['Sim Select']  = self.SSV.get() 
            self.VAR['Plot Select'] = self.PSV.get()
        print('self.VAR = '+str(self.VAR))
        print('Sim Select = '+str(self.VAR['Sim Select']) + ' and Plot Select = '+str(self.VAR['Plot Select']))
           
    def validate(self, value_if_allowed):
        if value_if_allowed:
            if any(x not in ['0','1','2','3','4','5','6','7','8','9',' ','','+','-','*','**','^','/','.'] for x in value_if_allowed):
                print('That ain\'t it chief')
                return False

            else:
                return True
        else:
            return True



    def MaintainAspect(self,event):
        XDIM = int(root.winfo_width() - root.winfo_width() % 16)
        YDIM = int(XDIM*9/16)
        root.geometry('{}x{}'.format(XDIM,YDIM))
     

    def RunSim(self):
        class FAKEVENT(object):    
            def __init__(self,event):
                self.keysym = event
        TMPEVENT = FAKEVENT('Return')
        self.Entry_Focus_Out(TMPEVENT)
        KGD2UI.Calculate_Dispersion(self)
        self.Result_Plotter()
        
    def Result_Plotter(self):
        if self.SSV.get() == self.SIMSEL[0]: #'Simple Photoinduced Charge'
            if self.PSV.get()   == self.PLTSEL[0][0]: #'Gaussian Laser Pulse'
                self.Graph_Plotter(self.Res_t,self.Res_Et[-1],'Time','Normalised Field Amplitude',None,None)
            elif self.PSV.get() == self.PLTSEL[0][1]: #'<a^2n+1>'
                #If SIMSEL[0], should ensure that a range of Et's are sampled for this result-currently unimplemented
                self.Graph_Plotter(float(self.VAR_Ncycx.get()),self.Res_a2disp,'Number of Optical Cycles','Vector Potential Momenta',None,None) 
            elif self.PSV.get() == self.PLTSEL[0][2]: #'Photoinduced Charge'
                self.Graph_Plotter(self.Res_Et,self.Res_Q,'Field Strength (F_0x) [Vm^-^1]','Photoinduced Charge',None,None)
            elif self.PSV.get() == self.PLTSEL[0][3]: #'<a^2n+1> Term Contributions to Charge'
                self.Graph_Plotter(float(self.VAR_Ncycx.get()),self.Res_Q,'Number of Optical Cycles','Photoinduced Charge per <a^2n+1> Term',None,None)
        
        elif self.SSV.get() == self.SIMSEL[1]:# 'Dispersion and Photoinduced Charge'
            None
        elif self.SSV.get() == self.SIMSEL[2]:# 'Dispersion Grapher' 
            None
    def Graph_Plotter(self,x,y,xAxName,yAxName,Legend,AXLIM):
        DIMx = np.array(x.shape)
        DIMy = np.array(y.shape)
        self.GraphPlot = plt.clf()
        #We find the index of the non matching dimension, allowing us to get the number of graphs to plot
        #Since if y.shape = (5,100) and x.shape = (100), we want y.shape[0] to be the range of plots
        if len(DIMx) == 1 and len(DIMy) != 1:
            Yn = DIMy[abs(np.where(DIMy==DIMx[0])[0][0]-1)]
                
            for n in range(Yn):
                self.GraphPlot = plt.plot(x,np.real(y[n]))
            
        elif len(DIMy) == 1 and len(DIMx) != 1:
            Xn = DIMx[abs(np.where(DIMx==DIMy[0])[0][0]-1)]
            for n in range(Xn):
                self.GraphPlot = plt.plot(x[n],np.real(y))
        else: 
            self.GraphPlot = plt.plot(x,np.real(y))
        
        if Legend != None:
            plt.Legend(Legend,location='best')
        plt.xlabel(xAxName)
        plt.ylabel(yAxName)
        AXLIM = 'tight'
        if AXLIM == 'tight':
            None
            plt.tight_layout()
            #plt.axvline(0)
            #plt.axhline(0)
        elif AXLIM == None:
            None
        self.GCanv.draw()
        plt.tight_layout()
                
        """
        SIMSEL = ['Simple Photoinduced Charge',
                    'Dispersion and Photoinduced Charge',
                    'Dispersion Grapher'];

        PLTSEL =      [['Gaussian Laser Pulse',        
               '<a^2n+1>',                     
               'Photoinduced Charge' ,         
               '<a^2n+1> Term Contributions to Charge']];
        PLTSEL.append(['Non Fourier Et(t)'           ,
               'E_t FFT Plot'                 ,
               'E_t Final Applied Dispersion' ,
               'E_t(t)^2n+1 After Dispersion' ,
               'Post Dispersion Photoinduced Charge', 
               'E_w IFFT Dispersion'          ,
               'Current Pulse Overlap at given Delta t',
               'Temporal Overlap Induced Charge',
               'Dual Polarised Induced Charge'])
        PLTSEL.append(['No Phi'                       ,
               'Phi0'                       ,  
               'Phi1'                        , 
               'Phi2'                        , 
               'Phi3'                        ,
               'All_Phi'])Wat
        """
    def Save_Preset(self):
        Preset = tk.filedialog.asksaveasfilename(filetypes=[('Jsson','*.json'),('Text File','*.txt'),('Just Show me Whatever','*.*')], defaultextension="*.json") 
        if Preset is None: # asksaveasfile return `None` if dialog closed with "cancel".
            return
        json.dump(self.VAR, open(Preset,'w'))
        
    def reset(self,event):
        print('Make this reset in the future')
        self.VAR_SimulationCompleted.set(False)
    
    def GetEntries(self):
        return(self.EntryList)
    
    def close(self):
        reply = tk.messagebox.askyesnocancel("Save Session Before Quitting?", "Save Session Before Quitting?")
        if reply == True:
            json.dump(self.VAR, open("SessionRestore.dat",'w'))
            print('Your session has been spared this time!')
            root.destroy()
            plt.close('all')
            
            sys.exit()
            

            
        elif reply == False:
            print('Your session is now lost to time!')
            root.destroy()
            plt.close('all')
            sys.exit()

        else:
            print('Why did you press this button in the first place?')
            None
        
        
        
    def Window_Exit_Event(self):
        #if tk.messagebox.askokcancel("Quit", "Do you want to quit?"):
            json.dump(self.VAR, open("SessionRestore.dat",'w'))
            root.destroy()
            sys.exit()
        


root = tk.Tk()

KurgGUI= KGUI(root)
root.protocol("WM_DELETE_WINDOW", KurgGUI.Window_Exit_Event)
root.mainloop()
