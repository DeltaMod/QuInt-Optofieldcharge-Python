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
               'Photoinduced Charge at F_x' ,         
               '<a^2n+1> Term Contributions to Charge at Ncycx',
               '<a^2n+1> Term Contributions to Charge at F_0x']];
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

MATS = ['SiO2','GaN']
# Make the lists readable through strings - This is done to reduce the number 
#of compute commands to "find" the correct index using an OptionsBox
SIMSEL = [str(n)+' - '+SIMSEL[n] for n in range(len(SIMSEL))]
             
for i in range(len(PLTSEL)):
    for j in range(len(PLTSEL[i])):
        PLTSEL[i][j] = str(j)+' - '+PLTSEL[i][j]

PLTSELLen = len(max(PLTSEL, key=len))
SIMSELLen = len(max(SIMSEL, key=len))
SelBoxLen = max([PLTSELLen,SIMSELLen])+10

VAR = {'Sim_Select':'0 - Simple Photoinduced Charge','Plot_Select': '0 - Gaussian Laser Pulse',
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
        
        
        self.DDS = ['Sim_Select','Plot_Select','Mat']
        self.EntryList = ['Sim_Select','Plot_Select','f_0x','f_0y', 'Ncycx', 'Ncycy', 'F_x', 'F_y',
            't1', 't2', 'L', 'Aeff', 'BPF',  'TDn', 'Delt1', 'Delt2',  'ORD',  'N', 'Mat']
        
        if 'SessionRestore.dat' in os.listdir():
            with open("SessionRestore.dat") as SessionRestore:
                pyresponse = json.loads(SessionRestore.read())
                print(pyresponse)
                self.VAR = pyresponse 
            for entry in self.EntryList:
                if entry not in list(self.VAR.keys()):
                    if entry=='Mat':
                        self.VAR[entry] = 'SiO2'
                    else:
                        self.VAR[entry] = '0'
                    
        else:
            self.VAR = {}
            for entry in self.EntryList:
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
        self.EB_Frame = tk.Frame(root,bg='ghostwhite')
        self.EB_Frame.grid(row = 1, column = 1,sticky = 'nwe')
        self.EB_Frame.bind('<Configure>',self.MaintainAspect)
        
        
        #Creating Simulation Box Frame
        self.RSFrame = tk.Frame(root,bg='ghostwhite')
        self.RSFrame.grid(row = 2, column = 1,sticky='nwe')
        #Setting Weights of all root relevant columms to be 1
       
                 

        
         #%% Listing all of the location and span variables used to cotrol the appearance of the UI
         #In the future, consider instead creating a matrix that uses exec to automatically generate these boxes based upon some matrix storage, like:
         #A = [Name,Type,Frame,text,relief,anchor,bg,row,column,columnspan,sticky]
        #Preset Frame
        ROWPOPOUT = 0 ; COLPOPOUT = 2; PSETSPAN = 2;
        ROWSAVE   = 0 ; COLSAVE   = 4;
        ROWLOAD   = 0 ; COLLOAD   = 6;
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
        
        #Sim_Select Frame
        ROWSIMSEL = 0; COLSIMSEL = 4; 
        ROWPLTSEL = 2; COLPLTSEL = 4;
        ROWRESET  = 4; COLRESET  = 5;
        ROWCLOSE  = 4; COLCLOSE  = 8;
        ROWRUN    = 4; COLRUN    = 0;
        #Column for unit labels
        COLUNITS = 11;
        
        
        #root.update()
        #print(self.EB_Frame.winfo_width())
        for i in range(0,COLCLOSE):
            self.RSFrame.grid_columnconfigure(i, weight=1)
        for j in range(0,ROWCLOSE):
            self.RSFrame.grid_rowconfigure(j, weight=1)
        for i in range(0,ROWSLD,2):
            self.EB_Frame.grid_rowconfigure(i, weight=1,uniform='fred')
        for j in range(0,9):
            if j not in [0,1,4,6]:
                self.EB_Frame.grid_columnconfigure(j, weight=1,uniform='fred')

        
        #%% Preset Frame
        ##Popout Graph Button
        self.Popout_Button = tk.Button(self.PSFrame, text='Popout Graph', command=self.RunSim)
        self.Popout_Button.grid(row = ROWPOPOUT,column = COLPOPOUT,columnspan=PSETSPAN,rowspan=PSETSPAN,sticky='nwse')
        
        #Save Preset Button
        self.save_preset_button = tk.Button(self.PSFrame, text='Save Preset', command=self.Save_Preset)
        self.save_preset_button.grid(row = ROWSAVE,column = COLSAVE,columnspan=PSETSPAN,sticky='w')
        
        #Load Preset Button
        self.load_preset_button = tk.Button(self.PSFrame, text='Load Preset', command=self.Load_Preset)
        self.load_preset_button.grid(row = ROWLOAD,column = COLLOAD,columnspan=PSETSPAN, sticky='w')
        
        #Select Material Dropdown
        self.VAR_Mat = tk.StringVar(root); self.VAR_Mat.set(str(self.VAR['Mat']))
        self.EB_Mat  = tk.OptionMenu(self.PSFrame, self.VAR_Mat, *MATS)
        self.EB_Mat.grid(row = ROWSAVE+2, column =COLLOAD,columnspan=PSETSPAN,sticky='we')
        self.EB_Mat.config(bg=ButtonBG,activebackground =ButtonABG)
        self.EB_Mat["menu"].config(bg=ButtonABG)
        tk.Label(self.PSFrame, text="Material",relief='flat',anchor='w',bg=LabelBG).grid(row = ROWSAVE+2, column =COLSAVE,columnspan=PSETSPAN,sticky='we')
        
        #%% Double Entries
        
        #Driving and Injection Labels
        tk.Label(self.EB_Frame, text="Driving",relief='flat',justify='center',anchor='n',bg=LabelBG).grid(row=ROWf_0x-1,column=COLf_0x,columnspan=DEntSpan,sticky='we')
        tk.Label(self.EB_Frame, text="Injection",relief='flat',anchor='n',bg=LabelBG).grid(row=ROWf_0y-1,column=COLf_0y,columnspan=DEntSpan,sticky='we')
        #x and y labels, just to make it absolutely clear:
        dROW = abs(ROWf_0x - ROWNcycx)
        xlbl = ['x','x','x','t1']
        ylbl = ['y','y','y','t2']
        mlbl = ['  ','  ','  ','  ']
        for i in range(4):
            tk.Label(self.EB_Frame, text=xlbl[i],relief='flat',anchor='e',bg=LabelBG).grid(row=dROW*i+ROWf_0x,column=COLf_0x-1,columnspan=DDISpan,sticky='e')
            tk.Label(self.EB_Frame, text=ylbl[i],relief='flat',anchor='e',bg=LabelBG).grid(row=dROW*i+ROWf_0y,column=COLf_0y-1,columnspan=DDISpan,sticky='e') 
            tk.Label(self.EB_Frame, text=mlbl[i],relief='flat',anchor='w',bg=LabelBG).grid(row=dROW*i+ROWf_0x,column=COLL+2,columnspan=1,sticky='w')
        Rowunits = ['THz','No Units','V^Åm^-1','fs','mm','10^-12m^2','max(A_w)/BFP','No Units','No Units','fs','No Units']
        for i in range(len(Rowunits)):
            tk.Label(tk.Label(self.EB_Frame, text=Rowunits[i],relief='flat',anchor='e',bg=LabelBG).grid(row=dROW*i+ROWf_0x,column=COLUNITS,columnspan=DDISpan,sticky='w'))
        
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
        tk.Label(self.EB_Frame, text="Laser Frequencies (f_0)",relief='flat',anchor='w',bg=LabelBG).grid(row=ROWf_0x,column=0,columnspan=LabSpan,sticky='we')
        
        self.VAR_f_0x   = tk.StringVar(root); self.VAR_f_0x.set(str(self.VAR['f_0x']))
        self.EB_f_0x = tk.Entry(self.EB_Frame,validate='key',validatecommand=(self.vcmd,'%P'),justify='center',width=EWDT,textvariable=self.VAR_f_0x)
        self.EB_f_0x.grid(row=ROWf_0x,column=COLf_0x,columnspan=DEntSpan,sticky='we')
        
        self.VAR_f_0y   = tk.StringVar(root); self.VAR_f_0y.set(str(self.VAR['f_0y']))
        self.EB_f_0y = tk.Entry(self.EB_Frame,justify='center',validate='key',validatecommand=self.vcmd,width=EWDT,textvariable=self.VAR_f_0y)
        self.EB_f_0y.grid(row=ROWf_0y,column=COLf_0y,columnspan=DEntSpan,sticky='we')
        
        tk.Label(self.EB_Frame, text="Optical Cycles (f_0)",relief='flat',anchor='w',bg=LabelBG).grid(row=ROWNcycx,column=0 ,columnspan=LabSpan,sticky='we') 
        self.VAR_Ncycx   = tk.StringVar(root); self.VAR_Ncycx.set(str(self.VAR['Ncycx']))
        self.EB_Ncycx = tk.Entry(self.EB_Frame,justify='center',validate='key',validatecommand=self.vcmd,width=EWDT,textvariable=self.VAR_Ncycx)
        self.EB_Ncycx.grid(row=ROWNcycx,column=COLNcycx,columnspan=DEntSpan,sticky='we')
        
        self.VAR_Ncycy   = tk.StringVar(root); self.VAR_Ncycy.set(str(self.VAR['Ncycy']))
        self.EB_Ncycy = tk.Entry(self.EB_Frame,justify='center',validate='key',validatecommand=self.vcmd,width=EWDT,textvariable=self.VAR_Ncycy)
        self.EB_Ncycy.grid(row=ROWNcycy,column=COLNcycy,columnspan=DEntSpan,sticky='we')
        
        tk.Label(self.EB_Frame, text="Optical Field Strength (F_)",relief='flat',anchor='w',bg=LabelBG).grid(row=ROWF_x,column=0 ,columnspan=LabSpan,sticky='we') 
        self.VAR_F_x   = tk.StringVar(root); self.VAR_F_x.set(str(self.VAR['F_x']))
        self.EB_F_x = tk.Entry(self.EB_Frame,justify='center',validate='key',validatecommand=self.vcmd,width=EWDT,textvariable=self.VAR_F_x)
        self.EB_F_x.grid(row=ROWF_x,column=COLF_x,columnspan=DEntSpan,sticky='we')
        
        self.VAR_F_y   = tk.StringVar(root); self.VAR_F_y.set(str(self.VAR['F_y']))
        self.EB_F_y = tk.Entry(self.EB_Frame,justify='center',validate='key',validatecommand=self.vcmd,width=EWDT,textvariable=self.VAR_F_y)
        self.EB_F_y.grid(row=ROWF_y,column=COLF_y,columnspan=DEntSpan,sticky='we')
        
        tk.Label(self.EB_Frame, text="Time Range (t)",relief='flat',anchor='w',bg=LabelBG).grid(row=ROWt1,column=0,columnspan=LabSpan,sticky='we')
        self.VAR_t1   = tk.StringVar(root); self.VAR_t1.set(str(self.VAR['t1']))
        self.EB_t1 = tk.Entry(self.EB_Frame,justify='center',validate='key',validatecommand=self.vcmd,width=EWDT,textvariable=self.VAR_t1)
        self.EB_t1.grid(row=ROWt1,column=COLt1,columnspan=DEntSpan,sticky='we')
        
        self.VAR_t2   = tk.StringVar(root); self.VAR_t2.set(str(self.VAR['t2']))
        self.EB_t2 = tk.Entry(self.EB_Frame,justify='center',validate='key',validatecommand=self.vcmd,width=EWDT,textvariable=self.VAR_t2)
        self.EB_t2.grid(row=ROWt2,column=COLt2,columnspan=DEntSpan,sticky='we')
        
        #%% Single Entries
        tk.Label(self.EB_Frame, text="Material Thickness (L)",relief='flat',anchor='w',bg=LabelBG).grid(row=ROWL,column=0,columnspan=LabSpan,sticky='we')
        self.VAR_L   = tk.StringVar(root); self.VAR_L.set(str(self.VAR['L']))
        self.EB_L = tk.Entry(self.EB_Frame,justify='center',validate='key',validatecommand=self.vcmd,width=EWDT,textvariable=self.VAR_L)
        self.EB_L.grid(row=ROWL,column=COLL,columnspan=SEntSpan,sticky='we')
        
        tk.Label(self.EB_Frame, text="Effective Area (Aeff)",relief='flat',anchor='w',bg=LabelBG).grid(row=ROWAeff,column=0,columnspan=LabSpan,sticky='we')
        self.VAR_Aeff   = tk.StringVar(root); self.VAR_Aeff.set(str(self.VAR['Aeff']))
        self.EB_Aeff = tk.Entry(self.EB_Frame,justify='center',validate='key',validatecommand=self.vcmd,width=EWDT,textvariable=self.VAR_Aeff)
        self.EB_Aeff.grid(row=ROWAeff,column=COLAeff,columnspan=SEntSpan,sticky='we')
        
        tk.Label(self.EB_Frame, text="Band Pass Filter",relief='flat',anchor='w',bg=LabelBG).grid(row=ROWBPF,column=0,columnspan=LabSpan,sticky='we')
        self.VAR_BPF = tk.StringVar(root); self.VAR_BPF.set(str(self.VAR['BPF']))
        self.EB_BPF = tk.Entry(self.EB_Frame,justify='center',validate='key',validatecommand=self.vcmd,width=EWDT,textvariable=self.VAR_BPF)
        self.EB_BPF.grid(row=ROWBPF,column=COLBPF,columnspan=SEntSpan,sticky='we')
        
        tk.Label(self.EB_Frame, text="Orders of <a^2n+1> (ORD)",relief='flat',anchor='w',bg=LabelBG).grid(row=ROWORD,column=0,sticky='we')
        self.VAR_ORD = tk.StringVar(root); self.VAR_ORD.set(str(self.VAR['ORD']))
        self.EB_ORD = tk.Entry(self.EB_Frame,justify='center',validate='key',validatecommand=self.vcmd,width=EWDT,textvariable=self.VAR_ORD)
        self.EB_ORD.grid(row=ROWORD,column=COLORD,columnspan=SEntSpan,sticky='we')
        
        tk.Label(self.EB_Frame, text="Sampling Points (N)",relief='flat',anchor='w',bg=LabelBG).grid(row=ROWN,column=0,sticky='w')
        self.VAR_N = tk.StringVar(root); self.VAR_N.set(str(self.VAR['N']))
        self.EB_N = tk.Entry(self.EB_Frame,justify='center',validate='key',validatecommand=self.vcmd,width=EWDT,textvariable=self.VAR_N)
        self.EB_N.grid(row=ROWN,column=COLN,columnspan=SEntSpan,sticky='we')
        
         #Temporal Delay
        self.VAR_TBTDn = tk.StringVar(root); self.VAR_TBTDn.set(str(self.VAR['TDn']))
        
        tk.Label(self.EB_Frame, text="Teporal Delay Δt                          ",relief='flat',anchor='w',bg=LabelBG).grid(row=ROWTDn,column=0,columnspan=LabSpan,sticky='we')
        tk.Label(self.EB_Frame, text="=",relief='flat',anchor='c',bg=LabelBG).grid(row=ROWTDn,column =COLL+1,columnspan=3,sticky='we')
        tk.Label(self.EB_Frame,justify='center',width=EWDT,textvariable=self.VAR_TBTDn).grid(row=ROWTDn,column=COLDelt2-1,columnspan=2,sticky='we')
        
        self.VAR_TDn = tk.StringVar(root); self.VAR_TDn.set(str(self.VAR['TDn']))
        self.EB_TDn = tk.Entry(self.EB_Frame,justify='center',validate='key',validatecommand=self.vcmd,width=EWDT,textvariable=self.VAR_TDn)
        self.EB_TDn.grid(row=ROWTDn,column=COLTDn,columnspan=2,sticky='we')
        #Temporal Delay
        
        tk.Label(self.EB_Frame, text="n_min",relief='flat',anchor='e',bg=LabelBG).grid(row=ROWDelt1,column=COLDelt1-1,sticky='swe')
        self.VAR_Delt1 = tk.StringVar(root); self.VAR_Delt1.set(str(self.VAR['Delt1']))
        self.EB_Delt1 = tk.Entry(self.EB_Frame,justify='center',validate='key',validatecommand=self.vcmd,width=EWDT,textvariable=self.VAR_Delt1)
        self.EB_Delt1.grid(row=ROWDelt1,column=COLDelt1,columnspan=2,sticky='wse')
        
        ##Frame Slider 
        self.DeltSlider = tk.Scale(self.EB_Frame,from_=-1, to=1, orient='horizontal',width=10,length=100,sliderlength=10)
        self.DeltSlider.grid(row=ROWSLD,column=COLSLD, columnspan=5, sticky='swe')
        
        ##Delt2
        self.VAR_Delt2 = tk.StringVar(root); self.VAR_Delt2.set(str(self.VAR['Delt2']))
        self.EB_Delt2 = tk.Entry(self.EB_Frame,justify='center',validate='key',validatecommand=self.vcmd,width=EWDT,textvariable=self.VAR_Delt2)
        self.EB_Delt2.grid(row=ROWDelt2,column=COLDelt2,columnspan=2,sticky='swe')
        
        
        #%% Populating the Simulation Selection Frame
        #Defining StringVar for SIMSEL Options Menu
        self.VAR_Sim_Select = tk.StringVar(root)
        self.VAR_Sim_Select.set(self.VAR['Sim_Select'])
        
        #Generating Simulation Selector options menu
        self.EB_Sim_Select= tk.OptionMenu(self.RSFrame, self.VAR_Sim_Select, *SIMSEL)
        tk.Label(self.RSFrame, text="Select a Simulation",relief='flat',anchor='w',bg=LabelBG).grid(row = ROWSIMSEL, column = 0,sticky='we')
        self.EB_Sim_Select.grid(row = ROWSIMSEL, column =COLSIMSEL,columnspan=10,sticky='ew')
        self.EB_Sim_Select.config(width = SelBoxLen,bg=ButtonBG,activebackground =ButtonABG)
        self.EB_Sim_Select["menu"].config(bg=ButtonABG)
        def PLT_DDGen(self):
            #Defining StringVar for PLTSEL Options Menu
            self.VAR_Plot_Select = tk.StringVar(root)
            self.VAR_Plot_Select.set(PLTSEL[int(self.VAR_Sim_Select.get()[0])][int(str(self.VAR['Plot_Select'])[0])])
            #Generating Plot_Selector options menu
            self.EB_Plot_Select= tk.OptionMenu(self.RSFrame, self.VAR_Plot_Select, *PLTSEL[int(self.VAR_Sim_Select.get()[0])])
            tk.Label(self.RSFrame, text="Select a Plot",relief ='flat',anchor='w',bg=LabelBG ).grid(row = ROWPLTSEL, column = 0,sticky='we')
            self.EB_Plot_Select.grid(row = ROWPLTSEL, column =COLPLTSEL,columnspan=5,sticky='ew')
            self.EB_Plot_Select.config(width = SelBoxLen,bg=ButtonBG,activebackground =ButtonABG)
            self.EB_Plot_Select["menu"].config(bg=ButtonABG)
            
        PLT_DDGen(self)
        
        def SIM_DDChange(*args):
            self.VAR['Sim_Select'] = self.VAR_Sim_Select.get()
            if int(self.VAR_Sim_Select.get()[0]) == 0:
                self.VAR_L.set(0)
                self.EB_L.config(state='disabled')
            elif int(self.VAR_Sim_Select.get()[0]) != 0:
                self.VAR_L.set(self.VAR['L'])
                self.EB_L.config(state='normal')
            self.VAR['Plot_Select'] = 0
            self.EB_Plot_Select.destroy() #Destroy old OptionsMenu so it can be re-created (there might be a better way to do this)
            PLT_DDGen(self)      #Re generate PLT OptionsMenu
            self.Entry_Focus_Out(None)
        def PLT_DDChange(*args):
            VAR['Plot_Select'] = self.VAR_Plot_Select.get()
            self.Entry_Focus_Out(None)
            print(self.VAR_SimulationCompleted.get())
            if self.VAR_SimulationCompleted.get() == True:
                self.Result_Plotter()

            
        # link function to change dropdown
        self.VAR_Sim_Select.trace('w',SIM_DDChange)
        self.VAR_Plot_Select.trace('w',PLT_DDChange)
        
        
        for entry in self.EntryList:
            self.__getattribute__('EB_'+entry).bind('<FocusOut>',self.Entry_Focus_Out)
            self.__getattribute__('EB_'+entry).bind('<Return>',self.Entry_Focus_Out)          
      
                
        #VAR = {'Sim_Select':0,'Plot_Select':0,
        #   'f_0x':0,'f_0y':0, 'Ncycx':0, 'Ncycy':0, 'F_x':0, 'F_y':0,
        #   't1':0, 't2':0, 'L':0, 'Aeff':0, 'BPF':0,  'TDn':0, 'Delt1':0, 'Delt2':0,  'ORD':0,  'N':0}
         
        ##Run Simulation Button
        self.RunSim_Button = tk.Button(self.RSFrame, text='Run Simulation', command=self.RunSim)
        self.RunSim_Button.grid(row = ROWRUN,column = COLRUN,rowspan=2,columnspan=2,sticky='nwse')
         
        #Close Button
        self.close_button = tk.Button(self.RSFrame, text='Close', command=self.close)
        self.close_button.grid(row = ROWCLOSE,column = COLCLOSE,sticky='e')
        
        #Reset Button
        self.reset_button = tk.Button(self.RSFrame, text='Reset', command=self.reset)
        self.reset_button.grid(row = ROWRESET,column = COLRESET,sticky='w')
    def Entry_Focus_Out(self,event):
        if event != None:    
            if event.keysym == 'Return':
                root.focus_set()
        for entry in self.EntryList:
            if entry not in self.DDS:
                if self.__getattribute__('VAR_'+entry).get()=='':
                    self.__getattribute__('VAR_'+entry).set('0')
                try: 
                    self.VAR[entry] = eval(self.__getattribute__('VAR_'+entry).get())
                    self.__getattribute__('VAR_'+entry).set(str(self.VAR[entry]))
                except:
                    self.__getattribute__('VAR_'+entry).set(str(self.VAR[entry]))
            else:   
                self.VAR['Sim_Select']  = self.VAR_Sim_Select.get() 
                self.VAR['Plot_Select'] = self.VAR_Plot_Select.get()
                self.VAR['Mat']         = self.VAR_Mat.get()
        print('self.VAR = '+str(self.VAR))
        print('Sim_Select = '+str(self.VAR['Sim_Select']) + ' and Plot_Select = '+str(self.VAR['Plot_Select']))
           
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
        if self.VAR_Sim_Select.get() == self.SIMSEL[0]: #'Simple Photoinduced Charge'
            if self.VAR_Plot_Select.get()   == self.PLTSEL[0][0]: #'Gaussian Laser Pulse'
                self.Graph_Plotter(self.Res_t,self.Res_Et[-1],'time [s]','Normalised Field Amplitude',None,None)
                
            elif self.VAR_Plot_Select.get() == self.PLTSEL[0][1]: #'<a^2n+1>'
                self.Graph_Plotter(self.Ncycxarray,self.Res_a2disp,'Number of Optical Cycles','Vector Potential Momenta',None,None) 
                
            elif self.VAR_Plot_Select.get() == self.PLTSEL[0][2]: #'Photoinduced Charge'
                self.Graph_Plotter(self.Ncycxarray,self.Res_Q.Sum,'Number of Optical Cycles','Photoinduced Charge at F_0x = '+str(self.VAR['F_x'])+' [V/m]',None,None)
            
            elif self.VAR_Plot_Select.get() == self.PLTSEL[0][3]: #'<a^2n+1> Term Contributions to Charge at Ncycx'
                self.Graph_Plotter(self.F_0xarray,self.Res_Q.Terms_ord,'Optical Field Strength','Photoinduced Charge per <a^2n+1> at Ncycx = '+str(self.VAR['Ncycx']),None,None)    
            
            elif self.VAR_Plot_Select.get() == self.PLTSEL[0][4]: #'<a^2n+1> Term Contributions to Charge at Ncycx'
                self.Graph_Plotter(self.Ncycxarray,self.Res_Q.Terms_ord,'Number of Optical Cycles','Photoinduced Charge per <a^2n+1> at F_0x = '+str(self.Res_VAR['F_x']),None,None)
        
        elif self.VAR_Sim_Select.get() == self.SIMSEL[1]:# 'Dispersion and Photoinduced Charge'
             if self.VAR_Plot_Select.get()   == self.PLTSEL[1][0]:  #0 - 'Non Fourier Et(t)',

                 self.Graph_Plotter(self.Res_t,self.Res_Et,'time [s]','Non Dispersed $E_t(t)$',None,None)

             elif self.VAR_Plot_Select.get()   == self.PLTSEL[1][1]:  #1 - 'E_t FFT Plot' , 1 - 'E_t FFT Plot' , 
                 self.Graph_Plotter(self.Res_w,self.Res_Ew,'\omega [rad/s]','Amplitude',None,None)

             elif  self.VAR_Plot_Select.get()   == self.PLTSEL[1][2]: #2 - 'E_t Final Applied Dispersion'  
                 #xlabel('Time [fs]');  ylabel('Amplitude')
                 self.Graph_Plotter(self.Res_t,self.Res_Et,'\omega [rad/s]','Amplitude',None,None)
                 pass
             elif  self.VAR_Plot_Select.get()   ==   self.PLTSEL[1][3]: # 3 - 'E_t(t)^2n+1 After Dispersion' 
                 #xlabel('Optical Cycle $\left(\frac{t-t_0}{T}\right)$'); ylabel('Vector Potential (and $a^{2n+1}$)')
                 pass
             
             elif  self.VAR_Plot_Select.get()   ==   self.PLTSEL[1][4]:  # 4 - 'Post Dispersion Photoinduced Charge'
                 #xlabel('Optical Field [Vm$^{-1}$]'); ylabel('Photinduced Charge [fC]')    
                 pass
             
             elif  self.VAR_Plot_Select.get()   ==   self.PLTSEL[1][5]:  # 5 - 'E_w IFFT Dispersion' 
                 #xlabel('Time [s]'); ylabel('Amplitude')
                 pass
             
             elif  self.VAR_Plot_Select.get()   ==   self.PLTSEL[1][6]:  # 6 - 'Current Pulse Overlap at given Delta t',
                 #xlabel('$\Delta t$ [s]'); ylabel('Sum of Vector Potential')
                 pass
             
             elif  self.VAR_Plot_Select.get()   ==   self.PLTSEL[1][7]:  # 7 - Temporal Overlap Induced Charge',
                  #xlabel('$\Delta t$ [fs]'); ylabel('$Q(\Delta t$ [C])')   
                 pass
             
             elif  self.VAR_Plot_Select.get()   ==   self.PLTSEL[1][8]:  # 8 - 'Dual Polarised Induced Charge'
                 #xlabel('$\Delta t$ [s]'); ylabel('Vector Potential Momenta')
                 pass

        elif self.VAR_Sim_Select.get() == self.SIMSEL[2]:# 'Dispersion Grapher' 
             if self.VAR_Plot_Select.get()   ==   self.PLTSEL[2][0]:  # 0 - 'No Phi'
                 #xlabel('Time [fs]'); ylabel('Amplitude')
                 pass
             
             if self.VAR_Plot_Select.get()   ==   self.PLTSEL[2][1]:  # 1 - 'Phi 0'
                 pass
             
             if self.VAR_Plot_Select.get()   ==   self.PLTSEL[2][2]:  # 2 - 'Phi 1'
                 pass
             
             if self.VAR_Plot_Select.get()   ==   self.PLTSEL[2][3]:  # 3 - 'Phi 2'
                 pass
             
             if self.VAR_Plot_Select.get()   ==   self.PLTSEL[2][4]:  # 4 - 'Phi 3'
                 pass
             
             if self.VAR_Plot_Select.get()   ==   self.PLTSEL[2][5]:  # 5 - 'All_Phi' 
                 pass
            
    def Graph_Plotter(self,x,y,xAxName,yAxName,Legend,AXLIM):
        x = np.asarray(x)
        y = np.asarray(y)
        
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
        Preset = tk.filedialog.asksaveasfilename(filetypes=[('Json','*.json'),('Text File','*.txt'),('Just Show me Whatever','*.*')], defaultextension="*.json") 
        if Preset is None: # asksaveasfile return `None` if dialog closed with "cancel".
            return
        json.dump(self.VAR, open(Preset,'w'))
    
    def Load_Preset(self):
        FILENAME = tk.filedialog.askopenfilename(initialdir = os.getcwd(),filetypes=[('Json','*.json'),('Text File','*.txt'),('Just Show me Whatever','*.*')], defaultextension="*.json")

        with open(FILENAME) as SessionLoad:
            pyresponse = json.loads(SessionLoad.read())
            print(pyresponse)
            self.VAR = pyresponse 
        for entry in self.EntryList:
            if entry not in list(self.VAR.keys()):
                if entry=='Mat':
                    self.VAR[entry] = 'SiO2'
                else:
                    self.VAR[entry] = '0'
            self.__getattribute__('VAR_'+entry).set(self.VAR[entry])
            
                
    def reset(self,event):
        print('Make this reset in the future')
        self.VAR_SimulationCompleted.set(False)
    
    def GetEntries(self):
        return(self.EntryList)
    
    def close(self):
        reply = tk.messagebox.askyesnocancel("Save Session Before Quitting?", "Save Session Before Quitting?")
        print(reply)
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
