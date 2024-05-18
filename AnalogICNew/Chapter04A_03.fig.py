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
    netlist=SpiceEditor(fileName+".fig.asc")
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
fileName="Chapter04A_03"

ExportInit(dirBuild,fileName)

graphNum=5
idVV=0
idVA=1
idM1=2
idM2=3
idLO=4
MPLInit()
MPLNewFigure(graphNum)

xVDD=5.0
xVG=2.5

xVinMin=0.0
xVinMax=xVDD
xVinStep=0.0001

xV1Min=0.0
xV1Max=xVDD
xV1Step=0.01

xV2Min=0.0
xV2Max=xVDD
xV2Step=0.5

lenV2=int(1+(xV2Max-xV2Min)/xV2Step)
listV2=[float(i)*xV2Step+xV2Min for i in range(lenV2)]
listV2str=["{:.1f}".format(x) for x in listV2]

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
    nl.set_element_model("M2","xPMOS W=2u L=1u")

for nl in [nl0]:
    nl.add_instruction(" ".join([".dc","Vin",str(xVinMin),str(xVinMax),str(xVinStep)]))

for nl in [nl1,nl2]:
    nl.set_component_value("Vout","0")
    nl.add_instruction(" ".join([".dc","Vout",str(xV1Min),str(xV1Max),str(xV1Step),"Vin",str(xV2Min),str(xV2Max),str(xV2Step)]))


d={}
rdat=SpiceRun(runner,nl0)
rdat1=SpiceRun(runner,nl1)
rdat2=SpiceRun(runner,nl2)

d["Vin"]=SpiceWave(rdat,"V(Vin)")
d["Vout"]=SpiceWave(rdat,"V(Vout)")

d["Av"]=NPDiff(d["Vin"],d["Vout"])

d["VBS1"]=0
d["VDS1"]=d["Vout"]
d["VGS1"]=d["Vin"]
d["VGT1"]=d["VGS1"]-xVT(d["VBS1"],xVT0N,xgammaN,xphiN)

d["VBS2"]=0
d["VDS2"]=xVDD-d["Vout"]
d["VGS2"]=xVDD-d["Vin"]
d["VGT2"]=d["VGS2"]-xVT(d["VBS2"],xVT0P,xgammaP,xphiP)

d["lVout"]=SpiceWave(rdat1,"V(Vout)")
d["lIdM1"]=SpiceWaveList(rdat1,"Id(M1)")
d["lIdM2"]=SpiceWaveList(rdat2,"Is(M2)")

iM1o=NPCross(d["VGT1"],0)
iM1s=NPCross(d["VGT1"],d["VDS1"])
iM2o=NPCross(d["VGT2"],0)
iM2s=NPCross(d["VGT2"],d["VDS2"])

iLoad=[]
for i in range(lenV2):
    iLoad.append(NPCross(d["lIdM1"][i],d["lIdM2"][i]))

plt.figure(idVV)
plt.plot(d["Vin"],d["Vout"],c="k",label=r"$v_{OUT}$")
plt.plot(d["Vin"],d["VGT1"],c="b",ls='dashed',lw=0.5,label=r"$M_1,sat$")
plt.plot(d["Vin"],xVDD-d["VGT2"],c="r",ls='dashed',lw=0.5,label=r"$M_2,sat$")
MPLDrawPoints(d["Vin"],d["Vout"],[iM1s,iM1o,iM2s,iM2o],"k","w",4)

plt.figure(idVA)
plt.plot(d["Vin"],d["Av"],c="k",label=r"$A_v$")
MPLDrawPoints(d["Vin"],d["Av"],[iM1o,iM2o],"k","w",4)

plt.figure(idM1)
plt.plot(d["Vin"],d["VDS1"],c="green",label=r"$v_{DS1}$")
plt.plot(d["Vin"],d["VGT1"],c="purple",label=r"$v_{GS1}-V_T$")
MPLDrawPoints(d["Vin"],d["VGT1"],[iM1s,iM1o],"k","w",4)

plt.figure(idM2)
plt.plot(d["Vin"],d["VDS2"],c="green",label=r"$v_{DS2}$")
plt.plot(d["Vin"],d["VGT2"],c="purple",label=r"$v_{GS2}-V_T$")
MPLDrawPoints(d["Vin"],d["VGT2"],[iM2s,iM2o],"k","w",4)

plt.figure(idLO)
cola1=(0.00,0.00,1.00)
colb1=(0.95,0.95,1.00)
cola2=(1.00,0.00,0.00)
colb2=(1.00,0.95,0.95)
for i in range(lenV2):
    r=float(i)/float(lenV2)
    ci1=tuple(x1*r+x2*(1-r) for x1,x2 in zip(cola1,colb1))
    ci2=tuple(x1*r+x2*(1-r) for x1,x2 in zip(cola2,colb2))
    label1={"label":r"$i_D(M_1)$"} if i==lenV2-1 else {}
    label2={"label":r"$i_D(M_2)$"} if i==lenV2-1 else {}
    plt.plot(d["lVout"],d["lIdM1"][i],c=ci1,lw=0.5,**label1)
    plt.plot(d["lVout"],d["lIdM2"][i],c=ci2,lw=0.5,**label2)
    
for i in range(lenV2):
    MPLDrawPoints(d["lVout"],d["lIdM2"][i],iLoad[i],"k","w",3) 
     
for id in range(graphNum):
    plt.figure(id)
    axes=MPLAxesInit()
    if id in[idVV,idM1,idM2,idVA,idLO]:
        axes.set_xlim(-0.2,5.2)
        axes.xaxis.set_major_locator(ticker.MultipleLocator(0.5))
        axes.set_xlabel(r"$v_{IN}~(\si{V})$")
    if id in[idVV,idM1,idM2]:
        axes.set_ylim(-0.2,5.2)
        axes.yaxis.set_major_locator(ticker.MultipleLocator(0.5))
        if id in[idVV]:
            axes.set_ylabel(r"$v_{OUT}~(\si{V})$")
        if id in[idM1,idM2]:
            axes.set_ylabel(r"$v~(\si{V})$")
        axes.legend(loc="upper right")
    if id in[idVA]:
        axes.set_ylim(-32,2)
        axes.yaxis.set_major_locator(ticker.MultipleLocator(5))
        axes.set_ylabel(r"$A_v$")
        axes.legend(loc="lower right")
    if id in[idLO]:
        axes.set_ylim(-0.02e-3,0.62e-3)
        axes.yaxis.set_major_locator(ticker.MultipleLocator(0.1e-3))
        axes.yaxis.set_major_formatter(lambda x, pos:"$"+"{:.1f}".format(x/1e-3)+"$")
        axes.set_ylabel(r"$i_D~(\si{mA})$")
        axes.legend(loc="upper left")
    plt.savefig(ExportNameGen(dirBuild,fileName,id),bbox_inches ="tight")

SpiceDelNet(fileName)
SpiceDelNet(fileName+"a")
SpiceDelNet(fileName+"b")