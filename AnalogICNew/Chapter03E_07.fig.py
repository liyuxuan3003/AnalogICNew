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
filename="Chapter03E_07"

fileGhost=os.path.join(folder,filename+".fig.pdf")
fileASC=filename+".fig.asc"
fileNET=filename+".fig.net"
def fileExport(id): return os.path.join(folder,filename+"_"+str(id)+".fig.pdf")

xNMOS=r".model xNMOS NMOS(LEVEL=1,VTO=0.7,KP=110e-6,GAMMA=0.4,LAMBDA=0.04,PHI=0.7,WD=0.0u,LD=0.016u,CGSO=220e-12,CGDO=220e-12,CGBO=700e-12,CJ=770e-6,CJSW=380e-12,MJ=0.5,MJSW=0.38)"
xPMOS=r".model xPMOS PMOS(LEVEL=1,VTO=-0.7,KP=50e-6,GAMMA=0.57,LAMBDA=0.05,PHI=0.8,WD=0.0u,LD=0.015u,CGSO=220e-12,CGDO=220e-12,CGBO=700e-12,CJ=560e-6,CJSW=350e-12,MJ=0.5,MJSW=0.38)"

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
LTspice.create_netlist(fileASC)
netlist=SpiceEditor(fileNET)
netlist.add_instruction(xNMOS)
netlist.add_instruction(xPMOS)
netlist.set_component_value("VDD","0")
netlist.set_element_model("M1","xNMOS W=0.8u L=0.8u")
netlist.set_element_model("M2","xNMOS W=0.8u L=0.8u")
netlist.set_element_model("M3","xPMOS W=0.8u L=0.8u")
netlist.set_element_model("M4","xPMOS W=0.8u L=0.8u")
netlist.set_component_value("R1","200k")
netlist.add_instruction(r".dc VDD 0 8 0.01")

graphNum=3
idVx=0
idIx=1
idSVV=2

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
d["VDD"]=SpiceWave(rdat,"V(VDD)")
d["VREF"]=SpiceWave(rdat,"V(VREF)")
d["VD1"]=SpiceWave(rdat,"V(VD1)")
d["VD2"]=SpiceWave(rdat,"V(VD2)")
d["ID1"]=SpiceWave(rdat,"Id(M1)")
d["ID2"]=SpiceWave(rdat,"Id(M2)")

d["dVREF"]=LDiff(d["VDD"],d["VREF"])
d["SVV"]=(d["VDD"]/d["VREF"])*d["dVREF"]

plt.figure(idVx)
plt.plot(d["VDD"],d["VDD"])
plt.plot(d["VDD"],d["VREF"])
plt.plot(d["VDD"],d["VD1"])
plt.plot(d["VDD"],d["VD2"])

plt.figure(idIx)
plt.plot(d["VDD"],d["ID1"])
plt.plot(d["VDD"],d["ID2"])

plt.figure(idSVV)
plt.plot(d["VDD"],d["SVV"])

for id in range(graphNum):
    plt.figure(id)
    axes=plt.gca()
    axes.grid(linewidth=0.25)
    axes.tick_params(labeltop=True,labelright=True,top=True,right=True,direction="in",width=0.5)
    axes.set_xlim(4.5,5.5)
    if id in [idSVV]:
        axes.set_ylim(0,0.1)
    if id in [idVx]:
        axes.set_ylim(0,2)
    plt.savefig(fileExport(id),bbox_inches ='tight')

os.remove(fileNET)