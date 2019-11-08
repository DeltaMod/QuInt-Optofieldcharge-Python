# -*- coding: utf-8 -*-
'''
Created on Fri Nov  8 13:31:20 2019

@author: Vidar Flodgren
'''

# -*- coding: utf-8 -*-
'''
FFAST-MPEG - A quicker than usual way to make gifs, trim videos, split videos and so on!
'''
import os, sys
import tkinter as tk
import matplotlib
from matplotlib.backends.backend_qt4agg import NavigationToolbar2QT
from matplotlib.figure import Figure
import numpy as np

#Define Colours to use for certain components
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

VARSTORE = {'Sim Select':0,'Plot Select':0,'F_x':0,'F_y':0}



        
class FFAST_MPEGUI:
    LABEL_TEXT = ['FFAST-MPEG']
    def __init__(self, master):
        self.master = master 
        master.title('Khurgin Py-UI')
        XDIM = 1280
        YDIM = 720
        root.geometry('{}x{}'.format(XDIM,YDIM))
        root.config(bg='aliceblue')
        #Creating Graph Frame
        self.GCanvas = tk.Frame(root, bg='white', width=YDIM, height=YDIM, relief = 'raised') # , 
        self.GCanvas.grid(row = 0, column = 0,  sticky='nwse')
        #Creating Variables Frame
        self.GFrame = tk.Frame(root,bg='ghostwhite',width=XDIM-YDIM,height=YDIM)
        self.GFrame.grid(row = 0, column = 1)
        self.GFrame.bind('<Configure>',self.MaintainAspect)
        #Setting Weights of all root relevant columms to be 1
        tk.Grid.rowconfigure(root, 0, weight=1)
        tk.Grid.columnconfigure(root, 0, weight=1)
        tk.Grid.columnconfigure(root, 1, weight=1)
        
                
        ##Convert Button
        self.Convert_Button = tk.Button(self.GFrame, text='Convert', command=self.convert)
        self.Convert_Button.grid(row = 1,column = 1,sticky='nwse')
      
        #Defining StringVar for SIMSEL Options Menu
        self.SSV = tk.StringVar(root)
        self.SSV.set(SIMSEL[VARSTORE['Sim Select']])
        #Generating the actual OptionsMenu
        self.DDSIMS= tk.OptionMenu(self.GFrame, self.SSV, *SIMSEL)
        tk.Label(self.GFrame, text="Select a Simulation",relief='flat',bg=LabelBG).grid(row = 9, column = 0)
        self.DDSIMS.grid(row = 10, column =0)
        self.DDSIMS.config(width = SelBoxLen,bg=ButtonBG,activebackground =ButtonABG)
        self.DDSIMS["menu"].config(bg=ButtonABG)
        def PLT_DDGen(self):
            #Defining StringVar for PLTSEL Options Menu
            self.PSV = tk.StringVar(root)
            self.PSV.set(PLTSEL[int(self.SSV.get()[0])][VARSTORE['Plot Select']])
            #Generating the actual OptionsMenu
            self.DDPLT= tk.OptionMenu(self.GFrame, self.PSV, *PLTSEL[int(self.SSV.get()[0])])
            tk.Label(self.GFrame, text="Select a Plot",relief ='flat',bg=LabelBG ).grid(row = 9, column = 1)
            self.DDPLT.grid(row = 10, column =1)
            self.DDPLT.config(width = SelBoxLen,bg=ButtonBG,activebackground =ButtonABG)
            self.DDPLT["menu"].config(bg=ButtonABG)
            
        PLT_DDGen(self)   
        
        def SIM_DDChange(*args):
            VARSTORE['Sim Select'] = int(self.SSV.get()[0])
            VARSTORE['Plot Select'] = 0
            self.DDPLT.destroy() #Destroy old OptionsMenu so it can be re-created (there might be a better way to do this)
            PLT_DDGen(self)      #Re generate PLT OptionsMenu
        def PLT_DDChange(*args):
            VARSTORE['Plot Select'] = int(self.PSV.get()[0])
            
        # link function to change dropdown
        self.SSV.trace('w',SIM_DDChange)
        self.PSV.trace('w',PLT_DDChange)
        
        
#        
        #Close Button
        self.close_button = tk.Button(self.GFrame, text='Close', command=self.close)
        self.close_button.grid(row = 11,column = 5,sticky='w')
        
    
    def MaintainAspect(self,event):
        XDIM = int(root.winfo_width() - root.winfo_width() % 16)
        YDIM = int(XDIM*9/16)
        root.geometry('{}x{}'.format(XDIM,YDIM))
     
     
    def convert(self,event):
        print(root.winfo_height())


    def close(self):
        print('Bye!')
        root.destroy()

root = tk.Tk()

FFASTGUI = FFAST_MPEGUI(root)
root.mainloop()