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

modelNMOSL1=r".model xNMOS NMOS(LEVEL=1,VTO=0.7,KP=110e-6,GAMMA=0.4,LAMBDA=0.04,PHI=0.7,WD=0.0u,LD=0.016u,CGSO=220e-12,CGDO=220e-12,CGBO=700e-12,CJ=770e-6,CJSW=380e-12,MJ=0.5,MJSW=0.38)"
modelPMOSL1=r".model xPMOS PMOS(LEVEL=1,VTO=-0.7,KP=50e-6,GAMMA=0.57,LAMBDA=0.05,PHI=0.8,WD=0.0u,LD=0.015u,CGSO=220e-12,CGDO=220e-12,CGBO=700e-12,CJ=560e-6,CJSW=350e-12,MJ=0.5,MJSW=0.38)"

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
fileName="Chapter04B_04"

ExportInit(dirBuild,fileName)

graphNum=8
idVV=0
idIV=1
idM={}
idM[1]=2
idM[2]=3
idM[3]=4
idM[4]=5
idM[5]=6
idAV=7
MPLInit()
MPLNewFigure(graphNum)

xVDD=5.0
xVICmin=0
xVICmax=xVDD
xVICstep=0.0005
xIREF=0.06e-3

runner=SpiceRunner(dirBuild)
nl=SpiceNetlist(fileName)
nl.add_instruction(modelNMOSL1)
nl.add_instruction(modelPMOSL1)
nl.set_component_value("VDD",str(xVDD))
nl.set_component_value("VID","0")
nl.set_component_value("VIC","0")
nl.set_component_value("IREF",str(xIREF))
nl.set_element_model("M1","xNMOS W=1u L=1u")
nl.set_element_model("M2","xNMOS W=1u L=1u")
nl.set_element_model("M3","xPMOS W=2u L=1u")
nl.set_element_model("M4","xPMOS W=2u L=1u")
nl.set_element_model("M5","xNMOS W=1u L=1u")
nl.set_element_model("M6","xNMOS W=1u L=1u")
nl.set_element_model("M7","xPMOS W=4u L=1u")
nl.add_instruction(" ".join([".dc","VIC",str(xVICmin),str(xVICmax),str(xVICstep)]))

d={}
rdat=SpiceRun(runner,nl)

d["Vin1"]=SpiceWave(rdat,"V(Vin1)")
d["Vin2"]=SpiceWave(rdat,"V(Vin2)")
d["Vin1a"]=SpiceWave(rdat,"V(Vin1a)")
d["Vin2a"]=SpiceWave(rdat,"V(Vin2a)")
d["VIC"]=d["Vin1"]
d["VID"]=d["Vin1"]-d["Vin1a"]
d["Vp"]=SpiceWave(rdat,"V(Vp)")
d["Vout1"]=SpiceWave(rdat,"V(Vout1)")
d["Vout2"]=SpiceWave(rdat,"V(Vout2)")
d["Vm3"]=SpiceWave(rdat,"V(Vm3)")
d["Vm5"]=SpiceWave(rdat,"V(Vm5)")
d["VDD"]=SpiceWave(rdat,"V(VDD)")
d["GND"]=np.zeros(len(d["VDD"]),dtype=np.float32)
d["ID1"]=SpiceWave(rdat,"Id(M1)")
d["ID2"]=SpiceWave(rdat,"Id(M2)")
d["ISS"]=SpiceWave(rdat,"Id(M5)")
d["IREF"]=SpiceWave(rdat,"Id(M6)")
d["AV"]=NPDiff(d["VIC"],d["Vout1"])

# G D S B
SatAnalyze("N",d,"1",d["Vin1"],d["Vout1"],d["Vp"],d["GND"])
SatAnalyze("N",d,"2",d["Vin2"],d["Vout2"],d["Vp"],d["GND"])
SatAnalyze("P",d,"3",d["Vm3"],d["Vout1"],d["VDD"],d["VDD"])
SatAnalyze("P",d,"4",d["Vm3"],d["Vout2"],d["VDD"],d["VDD"])
SatAnalyze("N",d,"5",d["Vm5"],d["Vp"],d["GND"],d["GND"])

iM1o=NPCross(d["VGT1"],0)
iM2o=NPCross(d["VGT2"],0)

iM1s=NPCross(d["VGT1"],d["VDS1"])
iM2s=NPCross(d["VGT2"],d["VDS2"])
iM3s=NPCross(d["VGT3"],d["VDS3"])
iM4s=NPCross(d["VGT3"],d["VDS4"])
iM5s=NPCross(d["VGT5"],d["VDS5"])

plt.figure(idVV)
plt.plot(d["VIC"],d["Vout1"],c='k',label="$v_{OUT1}=v_{OUT2}$")
plt.plot(d["VIC"],d["Vp"],c="gray",ls='dashed',label="$v_{P}$")
MPLDrawPoints(d["VIC"],d["Vout1"],[iM1o,iM1s],"k","w",4)
MPLDrawPoints(d["VIC"],d["Vout1"],[iM5s,iM3s],"k","w",2)

plt.figure(idIV)
plt.plot(d["VIC"],d["ID1"],c='k',label="$i_{D1}=i_{D2}$")
plt.plot(d["VIC"],d["ISS"],c='gray',label="$i_{SS}$")
MPLDrawPoints(d["VIC"],d["ISS"],[iM1o,iM1s],"k","w",4)
MPLDrawPoints(d["VIC"],d["ISS"],[iM5s,iM3s],"k","w",2)
MPLDrawPoints(d["VIC"],d["ID1"],[iM1o,iM1s],"k","w",4)
MPLDrawPoints(d["VIC"],d["ID1"],[iM5s,iM3s],"k","w",2)


for i in [1,2,3,4,5]:
    plt.figure(idM[i])
    plt.plot(d["VIC"],d["VDS"+str(i)],c="green",label="$v_{DS"+str(i)+"}$")
    plt.plot(d["VIC"],d["VGT"+str(i)],c="purple",label="$v_{GS"+str(i)+"}-V_T$")

plt.figure(idM[1])
MPLDrawPoints(d["VIC"],d["VGT1"],[iM1o,iM1s],"k","w",4)

plt.figure(idM[2])
MPLDrawPoints(d["VIC"],d["VGT2"],[iM2o,iM2s],"k","w",4)

plt.figure(idM[3])
MPLDrawPoints(d["VIC"],d["VGT3"],[iM3s],"k","w",4)

plt.figure(idM[4])
MPLDrawPoints(d["VIC"],d["VGT4"],[iM4s],"k","w",4)

plt.figure(idM[5])
MPLDrawPoints(d["VIC"],d["VGT5"],[iM5s],"k","w",4)

plt.figure(idAV)
plt.plot(d["VIC"],d["AV"],c="k",label="$A_{vc}$")


for id in range(graphNum):
    plt.figure(id)
    axes=MPLAxesInit()
    axes.set_xlim(-0.2,5.2)
    axes.xaxis.set_major_locator(ticker.MultipleLocator(0.5))
    axes.set_xlabel(r"$v_{IC}~(\si{V})$")
    if id in [idVV]+list(idM.values()):
        axes.set_ylim(-0.2,5.2)
        axes.yaxis.set_major_locator(ticker.MultipleLocator(0.5))
        axes.set_ylabel(r"$v~(\si{V})$")
        if id in [idVV]:
            axes.legend(loc="upper right")
        else:
            axes.legend(loc="upper left")
    if id in [idIV]:
        axes.set_ylim(-0.005e-3,0.085e-3)
        axes.yaxis.set_major_locator(ticker.MultipleLocator(0.01e-3))
        axes.yaxis.set_major_formatter(lambda x, pos:"$"+"{:.2f}".format(x/1e-3)+"$")
        axes.set_ylabel(r"$i~(\si{mA})$")
        axes.legend(loc="upper left")
    if id in [idAV]:
        axes.set_ylim(-0.85,0.05)
        axes.yaxis.set_major_locator(ticker.MultipleLocator(0.1))
        axes.set_ylabel(r"$A_{vc}$")
        axes.legend(loc="lower right")
    plt.savefig(ExportNameGen(dirBuild,fileName,id),bbox_inches ="tight")

SpiceDelNet(fileName)