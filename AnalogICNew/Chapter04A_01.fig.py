# !%config InlineBackend.figure_format = 'svg'

from PyLTSpice import *
from cycler import cycler
import matplotlib.pyplot as plt
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
    netlist=SpiceEditor(fileName+".fig.asc")
    return netlist

def SpiceDelNet(fileName):
    os.remove(fileName+".fig.net")

def SpiceRun(runner,netlist):
    raw,log=runner.run_now(netlist,run_filename="Temp.net")
    rdat=RawRead(raw)
    return rdat

def SpiceWave(rdat,wave):
    return rdat.get_trace(wave).get_wave(0)

#--------------------------------#

modelNMOSL1=r".model xNMOS NMOS(LEVEL=1,VTO=0.7,KP=110e-6,GAMMA=0.4,LAMBDA=0.04,PHIF=0.7,WD=0.0u,LD=0.016u,CGSO=220e-12,CGDO=220e-12,CGBO=700e-12,CJ=770e-6,CJSW=380e-12,MJ=0.5,MJSW=0.38)"
modelPMOSL1=r".model xPMOS PMOS(LEVEL=1,VTO=-0.7,KP=50e-6,GAMMA=0.57,LAMBDA=0.05,PHIF=0.8,WD=0.0u,LD=0.015u,NFS=6e11,CGSO=220e-12,CGDO=220e-12,CGBO=700e-12,CJ=560e-6,CJSW=350e-12,MJ=0.5,MJSW=0.38)"

xVT0N=0.7
xVT0P=0.7
xgammaN=0.40
xgammaP=0.57
xphiN=0.7
xphiP=0.8
def xVT(xVBS,xVT0,xgamma,xphi):
    return xVT0+xgamma*(np.sqrt(xphi-xVBS)-np.sqrt(xphi))

#--------------------------------#

dirBuild="./build"
fileName="Chapter04A_01"

ExportInit(dirBuild,fileName)

graphNum=4
idVV=0
idVA=1
idM1=2
idM2=3
MPLInit()
MPLNewFigure(graphNum)

xVDD=5
xV1Min=0
xV1Max=xVDD
xV1Step=0.0002

xV2Min=0
xV2Max=xVDD
xV2Step=0.5

runner=SpiceRunner(dirBuild)
nl0=SpiceNetlist(fileName)
nl1=SpiceNetlist(fileName+"a")
nl2=SpiceNetlist(fileName+"b")
for nl in [nl0,nl1,nl2]:
    nl.add_instruction(modelNMOSL1)
    nl.add_instruction(modelPMOSL1)
    nl.set_component_value("VDD",str(xVDD))
    nl.set_component_value("Vin","0")
    nl.set_element_model("M1","xNMOS W=1u L=1u")
    nl.set_element_model("M2","xPMOS W=1u L=1u")

for nl in [nl0]:
    nl.add_instruction(" ".join([".dc","Vin",str(xV1Min),str(xV1Max),str(xV1Step)]))

for nl in [nl1,nl2]:
    nl.set_component_value("Vout","0")
    nl.add_instruction(" ".join([".dc","Vout",str(xV1Min),str(xV1Max),str(xV1Step),".dc","Vin",str(xV2Min),str(xV2Max),str(xV2Step)]))


d={}
rdat=SpiceRun(runner,nl0)
rdat1=SpiceRun(runner,nl1)
rdat2=SpiceRun(runner,nl2)

d["Vin"]=SpiceWave(rdat,"V(Vin)")
d["Vout"]=SpiceWave(rdat,"V(Vout)")
d["Vout"][d["Vin"]<0.01]=d["Vout"][d["Vin"]==0.01]
# 滤掉Vin较小时Vout的异常波动

d["Av"]=NPDiff(d["Vin"],d["Vout"])

d["VBS1"]=0
d["VDS1"]=d["Vout"]
d["VGS1"]=d["Vin"]
d["VGT1"]=d["VGS1"]-xVT(d["VBS1"],xVT0N,xgammaN,xphiN)

d["VBS2"]=0
d["VDS2"]=xVDD-d["Vout"]
d["VGS2"]=xVDD-d["Vout"]
d["VGT2"]=d["VGS2"]-xVT(d["VBS2"],xVT0P,xgammaP,xphiP)

iM1o=NPCross(d["VGT1"],0)
iM1o=iM1o+50    # 微调一下位置，让点的位置更漂亮

iM1s=NPCross(d["VGT1"],d["VDS1"])

plt.figure(idVV)
plt.plot(d["Vin"],d["Vout"],c="k")
plt.plot(d["Vin"],d["VGT1"],c="b",ls='dashed',lw=0.5)
MPLDrawPoints(d["Vin"],d["Vout"],[iM1s,iM1o],"k","w",4)

plt.figure(idVA)
plt.plot(d["Vin"],d["Av"],c="k")
MPLDrawPoints(d["Vin"],d["Av"],[iM1s,iM1o],"k","w",4)

plt.figure(idM1)
plt.plot(d["Vin"],d["VDS1"],c="purple")
plt.plot(d["Vin"],d["VGT1"],c="green")
MPLDrawPoints(d["Vin"],d["VGT1"],[iM1s,iM1o],"k","w",4)

plt.figure(idM2)
plt.plot(d["Vin"],d["VDS2"],c="purple")
plt.plot(d["Vin"],d["VGT2"],c="green")

for id in range(graphNum):
    plt.figure(id)
    axes=MPLAxesInit()
    if id in[idVV,idM1,idM2,idVA]:
        axes.set_xlim(-0.2,5.2)
    if id in[idVV,idM1,idM2]:
        axes.set_ylim(-0.2,5.2)
    plt.savefig(ExportNameGen(dirBuild,fileName,id),bbox_inches ="tight")

SpiceDelNet(fileName)