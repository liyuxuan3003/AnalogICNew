# !%config InlineBackend.figure_format = 'svg'

from PyLTSpice import *
from cycler import cycler
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import os
import copy

import numpy as np

def LDiff(x,y):
    d=np.diff(y)/np.diff(x)
    d=np.append(d,d[-1])
    return d

def LDCross(y1,y2):
    return np.argmin(np.abs(y1-y2))

def SpiceRun(nl):
    raw,log=runner.run_now(nl,run_filename="Temp.net")
    rdat=RawRead(raw)
    return rdat

def SpiceWave(rdat,wave):
    return rdat.get_trace(wave).get_wave(0)

folder="./build"
filename="Chapter03D_04"

fileGhost=os.path.join(folder,filename+".fig.pdf")
fileASC=filename+".fig.asc"
fileNET=filename+".fig.net"
def fileExport(id): return os.path.join(folder,filename+"_"+str(id)+".fig.pdf")

xNMOS=r".model xNMOS NMOS(LEVEL=1,VTO=0.7,KP=110e-6,GAMMA=0.4,LAMBDA=0.04,PHIF=0.7,WD=0.0u,LD=0.016u,CGSO=220e-12,CGDO=220e-12,CGBO=700e-12,CJ=770e-6,CJSW=380e-12,MJ=0.5,MJSW=0.38)"
xPMOS=r".model xPMOS PMOS(LEVEL=1,VTO=-0.7,KP=50e-6,GAMMA=0.57,LAMBDA=0.05,PHIF=0.8,WD=0.0u,LD=0.015u,NFS=6e11,CGSO=220e-12,CGDO=220e-12,CGBO=700e-12,CJ=560e-6,CJSW=350e-12,MJ=0.5,MJSW=0.38)"

xColorList=['lightcoral','sandybrown','gold','yellowgreen','mediumaquamarine','skyblue','mediumpurple']
xColorCycle=cycler(color=xColorList)

xArrowProp=dict(facecolor='black', shrink=0.05,width=0.02,headlength=6,headwidth=3)

os.makedirs(folder,exist_ok=True)
open(fileGhost,"w")

plt.rc('text',usetex=True)
plt.rc('text.latex',preamble=" ".join([
    r'\usepackage{ctex}',
    r'\usepackage{amsmath}',
    r'\usepackage{physics}',
    r'\usepackage{siunitx}',
    r'\usepackage{xcolor}']))

plt.rc('xtick', labelsize=6)
plt.rc('ytick', labelsize=6)

plt.rcParams["axes.prop_cycle"]=xColorCycle

diffConst=0.001

runner=SimRunner(output_folder=folder,simulator=LTspice)
netlist=SpiceEditor(fileASC)
netlist.add_instruction(xNMOS)
netlist.add_instruction(xPMOS)
netlist.set_component_value("Vout","0")
netlist.set_component_value("IREF1","0.06m")
netlist.set_component_value("IREF2","0.06m")
netlist.set_element_model("M1","xNMOS W=0.8u L=0.8u")
netlist.set_element_model("M2","xNMOS W=0.8u L=0.8u")
netlist.set_element_model("M3","xNMOS W=0.8u L=3.2u")
netlist.set_element_model("M4","xNMOS W=0.8u L=0.8u")
netlist.set_element_model("M5","xNMOS W=0.8u L=0.8u")
netlist.add_instruction(r".dc Vout 0 5 0.01")

graphNum=4
idIV=0
idVD=1
idVSAT2=2
idVSAT4=3

for i in range(graphNum):
    plt.figure(i,figsize=(5,5*0.618))

xVT0=0.7
xgamma=0.4
xphi=0.7

def xVT(xVS):
    return xVT0+xgamma*(np.sqrt(xphi+xVS)-np.sqrt(xphi))

d={}

nl=copy.deepcopy(netlist)
rdat=SpiceRun(nl)
d["Vout"]=SpiceWave(rdat,"Vout")
d["Iout"]=SpiceWave(rdat,"Id(M2)")
d["Iref"]=SpiceWave(rdat,"Id(M1)")
d["VD1"]=SpiceWave(rdat,"V(VD1)")
d["VD2"]=SpiceWave(rdat,"V(VD2)")
d["VD3"]=SpiceWave(rdat,"V(VD3)")
d["VD4"]=SpiceWave(rdat,"V(VD4)")
d["VD5"]=SpiceWave(rdat,"V(VD5)")
d["VGT2"]=d["VD5"]-xVT0
d["VDS2"]=d["VD2"]
d["VGT4"]=d["VD3"]-d["VD2"]-xVT(d["VD2"])
d["VDS4"]=d["VD4"]-d["VD2"]

plt.figure(idIV)
plt.plot(d["Vout"],d["Iref"],c='r',ls='dashed',lw='0.75',label="$I_{REF}$")
plt.plot(d["Vout"],d["Iout"],c='k',label="$i_{OUT}$")

plt.figure(idVD)

plt.plot(d["Vout"],d["VD3"],c='r',label="$v_{D3}$",ls="dashed")
plt.plot(d["Vout"],d["VD1"],c='b',label="$v_{D1}$",ls="dashed")
plt.plot(d["Vout"],d["VD5"],c='r',label="$v_{D5}$",ls="dashdot")
plt.plot([],[],c='w',label="$~$")
plt.plot(d["Vout"],d["VD4"],c='r',label="$v_{D4}$")
plt.plot(d["Vout"],d["VD2"],c='b',label="$v_{D2}$")

plt.figure(idVSAT2)
plt.plot(d["Vout"],d["VDS2"],c='green',label="$v_{DS2}$")
plt.plot(d["Vout"],d["VGT2"],c='purple',label="$v_{GS2}-V_{T}$")

plt.figure(idVSAT4)
plt.plot(d["Vout"],d["VDS4"],c='green',label="$v_{DS4}$")
plt.plot(d["Vout"],d["VGT4"],c='purple',label="$v_{GS4}-V_{T}$")

for id in range(graphNum):
    plt.figure(id)
    axes=plt.gca()
    axes.grid(linewidth=0.25)
    axes.tick_params(labeltop=True,labelright=True,top=True,right=True,direction="in",width=0.5)
    if id in[idIV,idVD,idVSAT2,idVSAT4]:
        axes.set_xlim(-0.1,5.1)
        axes.xaxis.set_major_locator(ticker.MultipleLocator(0.5))
        axes.set_xlabel(r"$v_{OUT}~(\si{V})$")
        axes.vlines(2.0,-10,+10,colors='gray',lw=0.8,ls='dashed')
        if id in[idIV]:
            axes.set_ylim(-0.005e-3,0.085e-3)
            axes.yaxis.set_major_locator(ticker.MultipleLocator(0.01e-3))
            axes.yaxis.set_major_formatter(lambda x, pos:"$"+"{:.2f}".format(x/1e-3)+"$")
            axes.set_ylabel(r"$i_{OUT}~(\si{mA})$")
            axes.text(2.1,0.075e-3,r"$V_{\min}=2V_{ON}=\SI{2.0}{V}$",ha="left",va="center")
            axes.legend(loc="lower right")
        if id in[idVD,idVSAT2,idVSAT4]:
            axes.set_ylim(-0.1,5.1)
            axes.yaxis.set_major_locator(ticker.MultipleLocator(0.5))
            axes.set_ylabel(r"$v~(\si{V})$")
            if id in[idVD]:
                axes.legend(loc="upper left",ncol=3)
            if id in[idVSAT2,idVSAT4]:
                axes.legend(loc="upper left")
    plt.savefig(fileExport(id),bbox_inches ='tight')

os.remove(fileNET)