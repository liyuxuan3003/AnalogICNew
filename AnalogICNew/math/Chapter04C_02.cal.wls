#!/usr/bin/env wolframscript
(* ::Package:: *)

ClearAll["Global`*"]


f=10^(-15);
p=10^(-12);
n=10^(-9);
u=10^(-6);
m=10^(-3);
k=10^(3);
M=10^(6);


reqVDD=5;
reqAv=-50;
reqPdiss=1m;
reqCL=5p;
reqSR=10u^(-1);
reqVOUTmin=1.5;
reqVOUTmax=4.0;


vKx["n"]=110u;
vKx["p"]=50u;
vlambda["n"]=0.04;
vlambda["p"]=0.05;
vVT["n"]=0.7;
vVT["p"]=0.7;


MT[1]="n";
MT[2]="n";
MT[3]="p";


xKx=(vKx[MT[#]])&
xlambda=(vlambda[MT[#]])&
xVT=(vVT[MT[#]])&


xAv=-gm[1] Rout;
xRout=gds[3]^(-1);
xSR=ISS/CL;
xPdiss=ISS VDD;
xVOUTmin=VON[1]+VON[2];
xVOUTmax=VDD-VON[3];
xVON=(Sqrt[2 ID[#]/(Kx[#]WLr[#])])&;
xID=(ISS)&;
xgm=(Sqrt[2Kx[#]WLr[#] ID[#]])&;
xgds=(lambda[#]ID[#])&;


vISSbySR=ISS//.Solve[(xSR==reqSR)//.{CL->reqCL},ISS][[1]]


vISSbyPdiss=ISS//.Solve[(xPdiss==reqPdiss)//.{VDD->reqVDD},ISS][[1]]


N/@{vISSbySR,vISSbyPdiss}/u


desISS=100u;


desWLr[3]=WLr[3]//.Solve[
	(xVOUTmax==reqVOUTmax)//.
	{VDD->reqVDD,VON->xVON,ID->xID,Kx->xKx,ISS->desISS},WLr[3]][[1]]


desVON[3]=xVON[3]//.{ID->xID,Kx->xKx,WLr->desWLr,ISS->desISS}


desWLr[1]=WLr[1]//.Solve[
	(xAv==reqAv)//.
	{gm->xgm,gds->xgds,ID->xID,Rout->xRout,Kx->xKx,lambda->xlambda,ISS->desISS},WLr[1]][[1]]


desVON[1]=xVON[1]//.{ID->xID,Kx->xKx,WLr->desWLr,ISS->desISS}


desWLr[2]=WLr[2]//.Solve[
	(xVOUTmin==reqVOUTmin)//.
	{VDD->reqVDD,VON->xVON,ID->xID,Kx->xKx,ISS->desISS,WLr[1]->desWLr[1]},WLr[2]][[1]]


desVON[2]=xVON[2]//.{ID->xID,Kx->xKx,WLr->desWLr,ISS->desISS}


desVON[4]=desVON[1]+desVON[2]


desWLr[4]=WLr[4]//.Solve[((WLr[4]/WLr[2])==(VON[4]/VON[2])^(-2))//.{WLr[2]->desWLr[2],VON->desVON}][[1]]
