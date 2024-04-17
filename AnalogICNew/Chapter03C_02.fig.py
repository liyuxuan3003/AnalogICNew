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
filename="Chapter03C_02"

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
netlist.set_component_value("VG1","1.7")
netlist.set_component_value("VG2","2.7")
netlist.set_element_model("M1","xNMOS W=0.8u L=0.8u")
netlist.set_element_model("M2","xNMOS W=0.8u L=0.8u")
netlist.add_instruction(r".dc Vout 0 5 0.01")

graphNum=2
idIV=0
idVmid=1

for i in range(graphNum):
    plt.figure(i,figsize=(5,5*0.618))

d={}

nl=copy.deepcopy(netlist)
rdat=SpiceRun(nl)
d["Vout"]=SpiceWave(rdat,"V(Vout)")
d["Vmid"]=SpiceWave(rdat,"V(Vmid)")
d["Iout"]=SpiceWave(rdat,"Id(M1)")

plt.figure(idIV)
plt.plot(d["Vout"],d["Iout"],c='k')

plt.figure(idVmid)
plt.plot(d["Vout"],d["Vout"],c="k",label="$v_{D2}$")
plt.plot(d["Vout"],d["Vmid"],c='gray',label="$v_{D1}$")

for id in range(graphNum):
    plt.figure(id)
    axes=plt.gca()
    axes.grid(linewidth=0.25)
    axes.tick_params(labeltop=True,labelright=True,top=True,right=True,direction="in",width=0.3)
    axes.vlines(2.0,-10,+10,colors='gray',lw=0.8,ls='dashed')
    if id in [idIV,idVmid]:
        axes.set_xlim(-0.1,5.1)
        axes.xaxis.set_major_locator(ticker.MultipleLocator(0.5))
        axes.set_xlabel(r"$v_{OUT}~(\si{V})$")
        if id in [idIV]:
            axes.set_ylim(-0.005e-3,0.085e-3)
            axes.yaxis.set_major_locator(ticker.MultipleLocator(0.01e-3))
            axes.yaxis.set_major_formatter(lambda x, pos:"$"+"{:.2f}".format(x/1e-3)+"$")
            axes.text(2.1,0.075e-3,r"$v_{OUT}\geq 2V_{ON}$",ha="left",va="center")
            axes.set_ylabel(r"$i_{OUT}~(\si{mA})$")
        if id in [idVmid]:
            axes.set_ylim(-0.1,5.1)
            axes.yaxis.set_major_locator(ticker.MultipleLocator(0.5))
            axes.set_ylabel(r"$v_{D}~(\si{V})$")
            axes.legend(loc="upper left")
    plt.savefig(fileExport(id),bbox_inches ='tight')

os.remove(fileNET)