# !%config InlineBackend.figure_format = 'svg'

from PyLTSpice import *
from cycler import cycler
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np
import os
import copy

def MPLInit():
    xColorList=['lightcoral','sandybrown','gold','yellowgreen','mediumaquamarine','skyblue','mediumpurple']
    xColorCycle=cycler(color=xColorList)
    plt.rc('text',usetex=True)
    plt.rc('text.latex',preamble=" ".join([
        r'\usepackage{ctex}',
        r'\usepackage{amsmath}',
        r'\usepackage{physics}',
        r'\usepackage{siunitx}',
        r'\usepackage{xcolor}']))
    plt.rc('xtick',labelsize=6)
    plt.rc('ytick',labelsize=6)
    plt.rc("axes",prop_cycle=xColorCycle)

def MPLNewFigure(graphNum):
    for i in range(graphNum):
        plt.figure(i,figsize=(5,5*0.618))

def MPLAxesInit():
    axes=plt.gca()
    axes.grid(linewidth=0.25,which='major')
    axes.grid(linewidth=0.10,which='minor')
    axes.tick_params(labeltop=True,labelright=True,top=True,right=True,width=0.5,direction="in",which='major')
    axes.tick_params(labeltop=True,labelright=True,top=True,right=True,width=0.2,direction="in",which='minor')
    return axes

def MPLDrawPoints(x,y,iloc,cout,cin,size):
    dic=\
    {
        "marker":'o',
        "markersize":size,
        "c":cin,"mec":cout,
        "ls":''
    }
    plt.plot(x[iloc],y[iloc],**dic)

def ExportInit(dirBuild,fileName):
    fileGhost=os.path.join(dirBuild,fileName+".fig.pdf")
    os.makedirs(dirBuild,exist_ok=True)
    f=open(fileGhost,"w")
    f.close()

def ExportNameGen(dirBuild,fileName,id):
    return os.path.join(dirBuild,fileName+"_"+str(id)+".fig.pdf")

def NPDiff(x,y):
    d=np.diff(y)/np.diff(x)
    d=np.append(d,d[-1])
    return d

def NPCross(y1,y2):
    return np.argmin(np.abs(y1-y2))

def SpiceRunner(dirBuild):
    runner=SimRunner(output_folder=dirBuild,simulator=LTspice)
    return runner

def SpiceNetlist(fileName):
    LTspice.create_netlist(fileName+".fig.asc")
    netlist=SpiceEditor(fileName+".fig.net")
    return netlist

def SpiceDelNet(fileName):
    os.remove(fileName+".fig.net")

def SpiceRun(runner,netlist):
    raw,log=runner.run_now(netlist,run_filename="Temp.net")
    rdat=RawRead(raw)
    return rdat

def SpiceWave(rdat,name):
    return rdat.get_trace(name).get_wave(0)

def SpiceWaveList(rdat,name):
    waveList=[]
    i=0
    while True:
        wave=rdat.get_trace(name).get_wave(i)
        if list(wave)==[]:
            return waveList
        waveList.append(wave)
        i=i+1

#--------------------------------#

modelNMOSL1=r".model xNMOS NMOS(LEVEL=1,VTO=0.7,KP=110e-6,GAMMA=0.4,LAMBDA=0.04,PHI=0.7,WD=0.0u,LD=0.016u,NFS=7e11,CGSO=220e-12,CGDO=220e-12,CGBO=700e-12,CJ=770e-6,CJSW=380e-12,MJ=0.5,MJSW=0.38)"
modelPMOSL1=r".model xPMOS PMOS(LEVEL=1,VTO=-0.7,KP=50e-6,GAMMA=0.57,LAMBDA=0.05,PHI=0.8,WD=0.0u,LD=0.015u,NFS=6e11,CGSO=220e-12,CGDO=220e-12,CGBO=700e-12,CJ=560e-6,CJSW=350e-12,MJ=0.5,MJSW=0.38)"

xVT0={}
xgamma={}
xphi={}

xVT0["N"]=0.7
xVT0["P"]=0.7
xgamma["N"]=0.40
xgamma["P"]=0.57
xphi["N"]=0.7
xphi["P"]=0.8

def xVT(t,xVBS):
    return xVT0[t]+xgamma[t]*(np.sqrt(xphi[t]-xVBS)-np.sqrt(xphi[t]))

def SatAnalyze(t,d,suffix,xVG,xVD,xVS,xVB):
    sym=+1 if t=="N" else -1
    d["VGS"+suffix]=sym*(xVG-xVS)
    d["VDS"+suffix]=sym*(xVD-xVS)
    d["VBS"+suffix]=sym*(xVB-xVS)
    d["VT"+suffix]=xVT(t,d["VBS"+suffix])
    d["VGT"+suffix]=d["VGS"+suffix]-xVT(t,d["VBS"+suffix])

#--------------------------------#

dirBuild="./build"
fileName="Chapter04E_01"

ExportInit(dirBuild,fileName)

graphNum=3
idVV=0
idVA=1
idVI=2
MPLInit()
MPLNewFigure(graphNum)

xVDD=+2.5
xVSS=-2.5
xVG=0.0
xRL=10e3
xCL=1000e-12

xVinMin=xVSS
xVinMax=xVDD
xVinStep=0.0001

xV1Min=xVSS
xV1Max=xVDD
xV1Step=0.01

xV2Min=xVSS
xV2Max=xVDD
xV2Step=0.5

lenV2=int(1+(xV2Max-xV2Min)/xV2Step)
listV2=[float(i)*xV2Step+xV2Min for i in range(lenV2)]
listV2str=["{:.1f}".format(x) for x in listV2]

runner=SpiceRunner(dirBuild)
nl=SpiceNetlist(fileName)
nl.add_instruction(modelNMOSL1)
nl.add_instruction(modelPMOSL1)
nl.set_component_value("VDD",str(xVDD))
nl.set_component_value("VSS",str(xVSS))
nl.set_component_value("VG",str(xVG))
nl.set_component_value("RL",str(xRL))
nl.set_component_value("CL",str(xCL))
nl.set_component_value("Vin","0")
nl.set_element_model("M1","xNMOS W=1u L=1u")
nl.set_element_model("M2","xPMOS W=2u L=1u")

nl.add_instruction(" ".join([".dc","Vin",str(xVinMin),str(xVinMax),str(xVinStep)]))

d={}
rdat=SpiceRun(runner,nl)

d["Vin"]=SpiceWave(rdat,"V(Vin)")
d["Vout"]=SpiceWave(rdat,"V(Vout)")
d["VG"]=SpiceWave(rdat,"V(VG)")
d["Av"]=NPDiff(d["Vin"],d["Vout"])

d["ID1"]=SpiceWave(rdat,"Id(M1)")
d["ID2"]=SpiceWave(rdat,"Is(M2)")
d["IR"]=SpiceWave(rdat,"I(RL)")

plt.figure(idVV)
plt.plot(d["Vin"],d["Vout"],c="k",label=r"$v_{OUT}$")

plt.figure(idVA)
plt.plot(d["Vin"],d["Av"],c="k",label=r"$A_v$")

plt.figure(idVI)
plt.plot(d["Vin"],d["ID1"],c="r")
plt.plot(d["Vin"],d["ID2"],c="b")
plt.plot(d["Vin"],d["IR"],c="k",ls="dashed")

for id in range(graphNum):
    plt.figure(id)
    axes=MPLAxesInit()
    if id in[idVV,idVA]:
        axes.set_xlim(-2.7,2.7)
        axes.xaxis.set_major_locator(ticker.MultipleLocator(0.5))
        axes.set_xlabel(r"$v_{IN}~(\si{V})$")
    if id in[idVV]:
        axes.set_ylim(-2.7,2.7)
        axes.yaxis.set_major_locator(ticker.MultipleLocator(0.5))
        axes.set_ylabel(r"$v_{OUT}~(\si{V})$")
        axes.legend(loc="upper right")
    if id in[idVA]:
        axes.set_ylim(-37,2)
        axes.yaxis.set_major_locator(ticker.MultipleLocator(5))
        axes.set_ylabel(r"$A_v$")
        axes.legend(loc="lower right")
    plt.savefig(ExportNameGen(dirBuild,fileName,id),bbox_inches ="tight")

SpiceDelNet(fileName)