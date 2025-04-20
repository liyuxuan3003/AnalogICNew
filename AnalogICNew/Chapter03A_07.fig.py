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
filename="Chapter03A_07"

fileGhost=os.path.join(folder,filename+".fig.pdf")
fileASC=filename+".fig.asc"
fileNET=filename+".fig.net"
def fileExport(id): return os.path.join(folder,filename+"_"+str(id)+".fig.pdf")

xNMOS=r".model xNMOS NMOS(LEVEL=3,VTO=0.7,UO=660,DELTA=2.4,ETA=0.1,KAPPA=0.15,THETA=0.1,NSUB=3e16,TOX=140e-10,XJ=0.2u,WD=0u,LD=0.016u,NFS=7e11,CGSO=220e-12,CGDO=220e-12,CGBO=700e-12,CJ=770e-6,CJSW=380e-12,MJ=0.5,MJSW=0.38)"
xPMOS=r".model xPMOS PMOS(LEVEL=3,VTO=-0.7,UO=210,DELTA=1.25,ETA=0.1,KAPPA=2.5,THETA=0.1,NSUB=6e16,TOX=140e-10,XJ=0.2u,WD=0u,LD=0.015u,NFS=6e11,CGSO=220e-12,CGDO=220e-12,CGBO=700e-12,CJ=560e-6,CJSW=350e-12,MJ=0.5,MJSW=0.38)"

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
netlist.set_component_value("Vin",str(diffConst))
netlist.set_component_value("Vphi","0")
netlist.set_component_value("BVout","V=0")
netlist.set_element_model("M1","xNMOS W=0.8u L=0.8u")
netlist.add_instruction(r".dc Vphi 0 5 0.05")

graphNum=1
idRVphi=0

for i in range(graphNum):
    plt.figure(i,figsize=(5,5*0.618))

d={}

nl=copy.deepcopy(netlist)
rdat=SpiceRun(nl)
d["Vphi"]=SpiceWave(rdat,"V(Vphi)")
d["Id"]=SpiceWave(rdat,"Id(M1)")
d["RON"]=1/(d["Id"]/diffConst)

plt.figure(idRVphi)
plt.plot(d["Vphi"],d["RON"],c='k',label='$v_{out}=v_{in}=0$')

for id in range(graphNum):
    plt.figure(id)
    axes=plt.gca()
    axes.grid(linewidth=0.25)
    axes.tick_params(labeltop=True,labelright=True,top=True,right=True,direction="in",width=0.5)
    if id in[idRVphi]:
        axes.set_xlim(0,5)
        axes.xaxis.set_major_locator(ticker.MultipleLocator(0.5))
        axes.set_xlabel(r"$v_{\phi}~(\si{V})$")
        axes.set_ylim(-1e3,26e3)
        axes.yaxis.set_major_locator(ticker.MultipleLocator(5e3))
        axes.yaxis.set_major_formatter(lambda x, pos:"$"+"{:.0f}".format(x/1e3)+"$")
        axes.set_ylabel(r"$r_{ON}~(\si{k\ohm})$")
        axes.legend(loc="upper right")
    plt.savefig(fileExport(id),bbox_inches ='tight')

os.remove(fileNET)