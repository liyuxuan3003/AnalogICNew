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
fileName="Chapter04B_02"

ExportInit(dirBuild,fileName)

graphNum=9
idVV=0
idIV=1
idM={}
idM[1]=2
idM[2]=3
idM[3]=4
idM[4]=5
idM[5]=6
idVOD=7
idAV=8
MPLInit()
MPLNewFigure(graphNum)

xVDD=5.0
xVIC=3.0
xVIDmin=-1.0
xVIDmax=+1.0
xVIDstep=0.0002
xIREF=0.06e-3

runner=SpiceRunner(dirBuild)
nl=SpiceNetlist(fileName)
nl.add_instruction(modelNMOSL1)
nl.add_instruction(modelPMOSL1)
nl.set_component_value("VDD",str(xVDD))
nl.set_component_value("VID","0")
nl.set_component_value("VIC",str(xVIC))
nl.set_component_value("IREF",str(xIREF))
nl.set_element_model("M1","xNMOS W=1u L=1u")
nl.set_element_model("M2","xNMOS W=1u L=1u")
nl.set_element_model("M3","xPMOS W=2u L=1u")
nl.set_element_model("M4","xPMOS W=2u L=1u")
nl.set_element_model("M5","xNMOS W=1u L=1u")
nl.set_element_model("M6","xNMOS W=1u L=1u")
nl.set_element_model("M7","xPMOS W=4u L=1u")
nl.add_instruction(" ".join([".dc","VID",str(xVIDmin),str(xVIDmax),str(xVIDstep)]))

d={}
rdat=SpiceRun(runner,nl)

d["Vin1"]=SpiceWave(rdat,"V(Vin1)")
d["Vin2"]=SpiceWave(rdat,"V(Vin2)")
d["Vin1a"]=SpiceWave(rdat,"V(Vin1a)")
d["Vin2a"]=SpiceWave(rdat,"V(Vin2a)")
d["VID"]=d["Vin1"]-d["Vin1a"]
d["VIC"]=d["Vin1a"]
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
d["VOD"]=d["Vout1"]-d["Vout2"]
d["VOC"]=(d["Vout1"]+d["Vout2"])/2
d["AV"]=NPDiff(d["VID"],d["VOD"])
d["AV1"]=NPDiff(d["VID"],d["Vout1"])
d["AV2"]=NPDiff(d["VID"],d["Vout2"])

# G D S B
SatAnalyze("N",d,"1",d["Vin1"],d["Vout1"],d["Vp"],d["GND"])
SatAnalyze("N",d,"2",d["Vin2"],d["Vout2"],d["Vp"],d["GND"])
SatAnalyze("P",d,"3",d["Vm3"],d["Vout1"],d["VDD"],d["VDD"])
SatAnalyze("P",d,"4",d["Vm3"],d["Vout2"],d["VDD"],d["VDD"])
SatAnalyze("N",d,"5",d["Vm5"],d["Vp"],d["GND"],d["GND"])

iMs={}

for i in [1,2,3,4,5]:
    iMs[i]=NPCross(d["VGT"+str(i)],d["VDS"+str(i)])

iMs5mirror=len(d["VID"])-iMs[5]

plt.figure(idVV)
plt.plot(d["VID"],d["Vout1"],c='r',label="$v_{OUT1}$")
plt.plot(d["VID"],d["Vout2"],c='b',label="$v_{OUT2}$")
plt.plot(d["VID"],d["Vp"],c="gray",ls='dashed',label="$v_{P}$")
plt.plot(d["VID"],d["VOC"],c="black",ls='dotted',label="$v_{OC}$")
plt.plot(d["VID"],d["VIC"]-xVT0["N"],c="b",ls='dashed',lw=0.5,label=r"$M_{1,2},sat$")
plt.plot(d["VID"],d["VDD"]-0.7,c="r",ls='dashed',lw=0.5,label=r"$M_{3,4},sat$")
MPLDrawPoints(d["VID"],d["Vout1"],[iMs[i] for i in [1,2]],"k","w",4)
MPLDrawPoints(d["VID"],d["Vout2"],[iMs[i] for i in [1,2]],"k","w",4)
MPLDrawPoints(d["VID"],d["Vout1"],[iMs[i] for i in [3,4]],"k","w",2)
MPLDrawPoints(d["VID"],d["Vout2"],[iMs[i] for i in [3,4]],"k","w",2)

plt.figure(idIV)
plt.plot(d["VID"],d["ID1"],c='r',label="$i_{D1}$")
plt.plot(d["VID"],d["ID2"],c='b',label="$i_{D2}$")
plt.plot(d["VID"],d["ISS"],c='gray',label="$i_{SS}$")
MPLDrawPoints(d["VID"],d["ISS"],[iMs[5],iMs5mirror],"k","w",4)

for i in [1,2,3,4,5]:
    plt.figure(idM[i])
    plt.plot(d["VID"],d["VDS"+str(i)],c="green",label="$v_{DS"+str(i)+"}$")
    plt.plot(d["VID"],d["VGT"+str(i)],c="purple",label="$v_{GS"+str(i)+"}-V_T$")
    MPLDrawPoints(d["VID"],d["VGT"+str(i)],[iMs[i]],"k","w",4)
    if i==5:
        MPLDrawPoints(d["VID"],d["VGT"+str(i)],[iMs5mirror],"k","w",4)

plt.figure(idVOD)
plt.plot(d["VID"],d["VOD"],c="k",label="$v_{OD}$")
MPLDrawPoints(d["VID"],d["VOD"],[iMs[i] for i in [1,2]],"k","w",4)

plt.figure(idAV)
plt.plot(d["VID"],d["AV"],c="k",label="$A_{vd}$")
# plt.plot(d["VID"],d["AV2"],c="k",label="$A_{V2}$")

for id in range(graphNum):
    plt.figure(id)
    axes=MPLAxesInit()
    axes.set_xlim(-1.04,1.04)
    axes.xaxis.set_major_locator(ticker.MultipleLocator(0.2))
    axes.set_xlabel(r"$v_{ID}~(\si{V})$")
    if id in [idVV]+list(idM.values()):
        axes.set_ylim(-0.2,5.2)
        axes.yaxis.set_major_locator(ticker.MultipleLocator(0.5))
        axes.set_ylabel(r"$v~(\si{V})$")
        if id in [idVV]:
            axes.legend(loc="lower right",ncol=1)
        else:
            axes.legend(loc="upper left")
    if id in [idIV]:
        axes.set_ylim(-0.005e-3,0.085e-3)
        axes.yaxis.set_major_locator(ticker.MultipleLocator(0.01e-3))
        axes.yaxis.set_major_formatter(lambda x, pos:"$"+"{:.2f}".format(x/1e-3)+"$")
        axes.set_ylabel(r"$i~(\si{mA})$")
        axes.legend(loc="upper left",ncol=2)
    if id in [idVOD]:
        axes.set_ylim(-5.4,5.4)
        axes.yaxis.set_major_locator(ticker.MultipleLocator(1.0))
        axes.set_ylabel(r"$v_{OD}~(\si{V})$")
        axes.legend(loc="upper right")
    if id in [idAV]:
        axes.set_ylim(-83,3)
        axes.yaxis.set_major_locator(ticker.MultipleLocator(10))
        axes.set_ylabel(r"$A_{vd}$")
        axes.legend(loc="lower right")
    plt.savefig(ExportNameGen(dirBuild,fileName,id),bbox_inches ="tight")

SpiceDelNet(fileName)