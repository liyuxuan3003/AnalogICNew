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
filename="Chapter04A_01"

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
netlist.set_component_value("VDD","5")
netlist.set_component_value("Vin","0")
netlist.set_element_model("M1","xNMOS W=1u L=1u")
netlist.set_element_model("M2","xPMOS W=1u L=1u")
netlist.add_instruction(r".dc Vin 0 5 0.001")

graphNum=1
idVV=0

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
d["Vin"]=SpiceWave(rdat,"V(Vin)")
d["Vout"]=SpiceWave(rdat,"V(Vout)")

plt.figure(idVV)
plt.plot(d["Vin"],d["Vout"])

for id in range(graphNum):
    plt.figure(id)
    axes=plt.gca()
    axes.grid(linewidth=0.25)
    axes.tick_params(labeltop=True,labelright=True,top=True,right=True,direction="in",width=0.3)
    if id in[idVV]:
        axes.set_xlim(-0.2,5.2)
        axes.set_ylim(-0.2,5.2)
    plt.savefig(fileExport(id),bbox_inches ='tight')

os.remove(fileNET)