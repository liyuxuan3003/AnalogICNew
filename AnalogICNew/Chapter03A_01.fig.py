# !%config InlineBackend.figure_format = 'svg'

from PyLTSpice import *
from cycler import cycler
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import os

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
filename="Chapter03A_01"

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

xVTHN0=+0.7
xVTHP0=-0.7
xVDD=5

runner=SimRunner(output_folder=folder,simulator=LTspice)
netlist=SpiceEditor(fileASC)
netlist.add_instruction(xNMOS)
netlist.add_instruction(xPMOS)
netlist.set_component_value("VG",str(0))
netlist.set_component_value("VS",str(0))
netlist.set_component_value("VD",str(0))
netlist.set_element_model("M1","xNMOS W=0.8u L=0.8u")

graphNum=6
idIVD=0
idRVD=1
idIVG=2
idRVG=3
idIVDext=4
idRVDext=5

for i in range(graphNum):
    plt.figure(i,figsize=(5,5*0.618))

simconfig=\
{
    "VD":
    {
        "simcmd":r".dc VD 0 2 0.01",
        "para":"VG",
        "paraval":[5,2.5],
        "paralatex":r"V_G",
        "col":["#000000","#777777"]
    },
    "VG":
    {
        "simcmd":r".dc VG 0 5 0.01",
        "para":"VD",
        "paraval":[2,1,0.5],
        "paralatex":r"V_D",
        "col":["#000000","#777777","#BBBBBB"]
    },
    "VDext":
    {
        "simcmd":r".dc VD 0 5 0.01",
        "para":"VG",
        "paraval":[5],
        "paralatex":r"V_G",
        "col":["#000000"]
    }
}

dNameList=["ID","VG","VD","VS","RON"]
d=\
{
    simname:
    {
        paraval:
        {
            dName:None for dName in dNameList
        } for paraval in sim["paraval"]
    } for simname,sim in simconfig.items()
}

diffconst=0.01
    

for simname,sim in simconfig.items():
    nl=netlist
    nl.add_instruction(sim["simcmd"])
    for paraval in sim["paraval"]:
        nl.set_component_value(sim["para"],str(paraval))
        rdat=SpiceRun(nl)
        rdatDiff=None
        if simname=="VG":
            nlDiff=nl
            nlDiff.set_component_value(sim["para"],str(paraval+diffconst))
            rdatDiff=SpiceRun(nlDiff)
        # dt -> dtemp
        dt=d[simname][paraval]
        dt["ID"]=SpiceWave(rdat,"Id(M1)")
        dt["VG"]=SpiceWave(rdat,"V(VG)")
        dt["VD"]=SpiceWave(rdat,"V(VD)")
        dt["VS"]=SpiceWave(rdat,"V(VS)")
        if simname in ["VD","VDext"]:
            dt["RON"]=1/LDiff(dt["VD"],dt["ID"])
        if simname in ["VG"]:
            dt["RON"]=1/((SpiceWave(rdatDiff,"Id(M1)")-dt["ID"])/diffconst)
        d[simname][paraval]=dt
    
for paraval,col in zip(simconfig["VD"]["paraval"],simconfig["VD"]["col"]):
    dt=d["VD"][paraval]
    plt.figure(idIVD)
    plt.plot(dt["VD"],dt["ID"],c=col)
    plt.figure(idRVD)
    plt.plot(dt["VD"],dt["RON"],c=col)

for paraval,col in zip(simconfig["VG"]["paraval"],simconfig["VG"]["col"]):
    dt=d["VG"][paraval]
    plt.figure(idIVG)
    plt.plot(dt["VG"],dt["ID"],c=col)
    plt.figure(idRVG)
    plt.plot(dt["VG"],dt["RON"],c=col)

for paraval,col in zip(simconfig["VDext"]["paraval"],simconfig["VDext"]["col"]):
    dt=d["VDext"][paraval]
    plt.figure(idIVDext)
    plt.plot(dt["VD"],dt["ID"],c=col)
    plt.figure(idRVDext)
    plt.plot(dt["VD"],dt["RON"],c=col)

for id in range(graphNum):
    plt.figure(id)
    axes=plt.gca()
    axes.grid(linewidth=0.25)
    axes.tick_params(labeltop=True,labelright=True,top=True,right=True,direction="in",width=0.5)
    if id in[idIVD,idRVD]:
        axes.set_xlim(0,2)
        axes.xaxis.set_major_locator(ticker.MultipleLocator(0.25))
        axes.set_xlabel(r"$V_{D}~(\si{V})$")
    if id in[idIVDext,idRVDext]:
        axes.set_xlim(0,5)
        axes.xaxis.set_major_locator(ticker.MultipleLocator(0.50))
        axes.set_xlabel(r"$V_{D}~(\si{V})$")
    if id in[idIVG,idRVG]:
        axes.set_xlim(0,5)
        axes.xaxis.set_major_locator(ticker.MultipleLocator(0.50))
        axes.set_xlabel(r"$V_{G}~(\si{V})$")
    if id in[idRVD,idRVG,idRVDext]:
        axes.set_ylim(-0.5e3,20.5e3)
        axes.yaxis.set_major_locator(ticker.MultipleLocator(2e3))
        axes.yaxis.set_major_formatter(lambda x, pos:"$"+"{:.0f}".format(x/1e3)+"$")
        axes.set_ylabel(r"$r_{ON}~(\si{k\ohm})$")
    if id in[idIVD,idIVG,idIVDext]:
        axes.set_ylim(-0.05e-3,1.25e-3)
        axes.yaxis.set_major_locator(ticker.MultipleLocator(0.2e-3))
        axes.yaxis.set_major_formatter(lambda x, pos:"$"+"{:.1f}".format(x/1e-3)+"$")
        axes.set_ylabel(r"$I_{D}~(\si{mA})$")
    plt.savefig(fileExport(id),bbox_inches ='tight')

os.remove(fileNET)