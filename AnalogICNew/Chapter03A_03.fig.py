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
filename="Chapter03A_03"

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

runner=SimRunner(output_folder=folder,simulator=LTspice)
netlist=SpiceEditor(fileASC)
netlist.add_instruction(xNMOS)
netlist.add_instruction(xPMOS)
netlist.set_component_value("Vphi","0")
netlist.set_component_value("Vin","1")
netlist.set_component_value("C1","200f ic=1")
netlist.set_element_model("M1","xNMOS W=0.8u L=0.8u")

graphNum=2
idVphi=0
idVc=1

for i in range(graphNum):
    plt.figure(i,figsize=(5,5*0.618))

simconfig=\
{
    "FAST":
    {
        "simcmd":r".tran 0.001n 16n",
        "pulse":r"PULSE(5 0 2n 0.2n 0.2n 20n 50n 1)"
    },
    "SLOW":
    {
        "simcmd":r".tran 0.001n 16n",
        "pulse":r"PULSE(5 0 2n 10n 10n 20n 50n 1)"
    }
}

dNameList=["I","Vphi","Vin","Vc","t"]
d={simname:{} for simname,sim in simconfig.items()}

for simname,sim in simconfig.items():
    nl=copy.deepcopy(netlist)
    nl.add_instruction(sim["simcmd"])
    nl.set_component_value("Vphi",sim["pulse"])
    rdat=SpiceRun(nl)
    dt=d[simname]
    dt["I"]=SpiceWave(rdat,"Id(M1)")
    dt["Vphi"]=SpiceWave(rdat,"V(Vphi)")
    dt["Vin"]=SpiceWave(rdat,"V(Vin)")
    dt["Vc"]=SpiceWave(rdat,"V(Vc)")
    dt["t"]=SpiceWave(rdat,"time")
    d[simname]=dt
    
plt.figure(idVphi)
plt.plot(d["FAST"]["t"],d["FAST"]["Vin"],c='gray',label='$v_{in}$')
plt.plot(d["FAST"]["t"],d["FAST"]["Vphi"],c='r',label='$v_{\phi,fast}$')
plt.plot(d["SLOW"]["t"],d["SLOW"]["Vphi"],c='b',ls='dashed',label='$v_{\phi,slow}$')

plt.figure(idVc)
plt.plot(d["FAST"]["t"],d["FAST"]["Vin"],c='gray',label='$v_{in}$')
plt.plot(d["FAST"]["t"],d["FAST"]["Vc"],c='r',label='$v_{out,fast}$')
plt.plot(d["SLOW"]["t"],d["SLOW"]["Vc"],c='b',ls='dashed',label='$v_{out,slow}$')

for id in range(graphNum):
    plt.figure(id)
    axes=plt.gca()
    axes.grid(linewidth=0.25)
    axes.tick_params(labeltop=True,labelright=True,top=True,right=True,direction="in",width=0.5)
    if id in[idVphi,idVc]:
        axes.set_xlim(-0.5e-9,16.5e-9)
        axes.xaxis.set_major_locator(ticker.MultipleLocator(2e-9))
        axes.xaxis.set_major_formatter(lambda x, pos:"$"+"{:.0f}".format(x/1e-9)+"$")
        axes.set_xlabel(r"$t~(\si{ns})$")
        axes.set_ylabel(r"$v~(\si{V})$")
    if id in [idVphi]:
        axes.set_ylim(-0.2,5.2)
        axes.yaxis.set_major_locator(ticker.MultipleLocator(1))
        axes.legend(loc="upper right")
    if id in [idVc]:
        axes.set_ylim(0.978,1.012)
        axes.yaxis.set_major_locator(ticker.MultipleLocator(0.005))
        axes.legend(loc="upper right")
    plt.savefig(fileExport(id),bbox_inches ='tight')

os.remove(fileNET)