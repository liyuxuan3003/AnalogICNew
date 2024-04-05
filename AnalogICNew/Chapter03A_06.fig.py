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
filename="Chapter03A_06"

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

simconfig=\
{
    "FA":"V=V(Vin)",
    "SA":"V=0"
}

runner=SimRunner(output_folder=folder,simulator=LTspice)
netlist=SpiceEditor(fileASC)
netlist.add_instruction(xNMOS)
netlist.add_instruction(xPMOS)
netlist.set_component_value("Vin","0")
netlist.set_component_value("Vindiff","0")
netlist.set_component_value("Vphi","5")
netlist.set_component_value("BVout","V=0")
netlist.set_element_model("M1","xNMOS W=0.8u L=0.8u")
netlist.set_element_model("M2","xPMOS W=2.4u L=0.8u")
netlist.add_instruction(r".dc Vin 0 5 0.001")

graphNum=1
idRVin=0

for i in range(graphNum):
    plt.figure(i,figsize=(5,5*0.618))

d={simname:{} for simname,sim in simconfig.items()}

for simname,sim in simconfig.items():
    nl=copy.deepcopy(netlist)
    nl.set_component_value("BVout",sim)
    rdat=SpiceRun(nl)

    nldiff=copy.deepcopy(nl)
    nldiff.set_component_value("Vindiff",str(diffConst))
    rdatdiff=SpiceRun(nldiff)

    dt=d[simname]
    dt["I1"]=SpiceWave(rdat,"Id(M1)")
    dt["I2"]=SpiceWave(rdat,"Is(M2)")
    dt["I1x"]=SpiceWave(rdatdiff,"Id(M1)")
    dt["I2x"]=SpiceWave(rdatdiff,"Is(M2)")
    dt["I"]=dt["I1"]+dt["I2"]
    dt["Ix"]=dt["I1x"]+dt["I2x"]
    dt["Vin"]=SpiceWave(rdat,"V(Vin)")

    dt['RON']=1/((dt["Ix"]-dt["I"])/diffConst)
    dt['RON1']=1/((dt["I1x"]-dt["I1"])/diffConst)
    dt['RON2']=1/((dt["I2x"]-dt["I2"])/diffConst)

    d[simname]=dt

plt.figure(idRVin)
plt.plot(d["SA"]["Vin"],d["SA"]["RON"],c='gray',label=r'$v_{out}=0$')
plt.plot(d["FA"]["Vin"],np.abs(d["FA"]["RON1"]),c='b',ls='dashed',label=r'$v_{out}=v_{in},~r_{on1}$')
plt.plot(d["FA"]["Vin"],np.abs(d["FA"]["RON2"]),c='r',ls='dashed',label=r'$v_{out}=v_{in},~r_{on2}$')
plt.plot(d["FA"]["Vin"],d["FA"]["RON"],c='k',label=r'$v_{out}=v_{in}$')


for id in range(graphNum):
    plt.figure(id)
    axes=plt.gca()
    axes.grid(linewidth=0.25)
    axes.tick_params(labeltop=True,labelright=True,top=True,right=True,direction="in",width=0.3)
    if id in[idRVin]:
        axes.set_xlim(0,5)
        axes.xaxis.set_major_locator(ticker.MultipleLocator(0.5))
        axes.set_xlabel(r"$v_{in}~(\si{V})$")
        axes.set_ylim(-0.2e3,6.2e3)
        axes.yaxis.set_major_locator(ticker.MultipleLocator(1e3))
        axes.yaxis.set_major_formatter(lambda x, pos:"$"+"{:.0f}".format(x/1e3)+"$")
        axes.set_ylabel(r"$r_{ON}~(\si{k\ohm})$")
        axes.legend(loc="upper right")
    plt.savefig(fileExport(id),bbox_inches ='tight')

os.remove(fileNET)