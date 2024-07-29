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
reqAv=100;
reqPdiss=1m;
reqFm3dB=100k
reqCL=5p;
reqSR=10u^(-1);
reqVINmin=1.0;
reqVINmax=4.5;


vKx["n"]=110u;
vKx["p"]=50u;
vlambda["n"]=0.04;
vlambda["p"]=0.05;
vVT["n"]=0.7;
vVT["p"]=0.7;


MT[1]="n";
MT[2]="n";
MT[3]="p";
MT[4]="p";
MT[5]="n";


xKx=(vKx[MT[#]])&
xlambda=(vlambda[MT[#]])&
xVT=(vVT[MT[#]])&


xAv=gm[2] Rout;
xRout=(gds[2]+gds[4])^(-1);
xSR=ISS/CL;
xPdiss=ISS VDD;
xFm3dB=1/(2 Pi Rout CL);
xVINmin=VON[2]+VON[5]+VT[2];
xVINmax=VDD-VON[4]-VT[4]+VT[2];
xVON=(Sqrt[2 ID[#]/(Kx[#]WLr[#])])&;
xID=(If[#==5,ISS,ISS/2])&;
xgm=(Sqrt[2Kx[#]WLr[#] ID[#]])&;
xgds=(lambda[#]ID[#])&;


vISSbySR=ISS//.Solve[(xSR==reqSR)//.{CL->reqCL},ISS][[1]]


vISSbyPdiss=ISS//.Solve[(xPdiss==reqPdiss)//.{VDD->reqVDD},ISS][[1]]


vISSbyFm3dB=ISS//.Solve[(xFm3dB==reqFm3dB)//.{VDD->reqVDD,Rout->xRout,CL->reqCL,gds->xgds,lambda->xlambda,ID->xID},ISS][[1]]


N/@{vISSbySR,vISSbyPdiss,vISSbyFm3dB}/u


desISS=100u;


desWLr[4]=WLr[4]//.Solve[
	(xVINmax==reqVINmax)//.
	{VDD->reqVDD,VON->xVON,ID->xID,Kx->xKx,ISS->desISS,VT->xVT},WLr[4]][[1]]


desVON[4]=xVON[4]//.{ID->xID,Kx->xKx,WLr->desWLr,ISS->desISS}


desWLr[2]=WLr[2]//.Solve[
	(xAv==reqAv)//.
	{gm->xgm,gds->xgds,ID->xID,Rout->xRout,Kx->xKx,lambda->xlambda,ISS->desISS},WLr[2]][[1]]


desVON[2]=xVON[2]//.{ID->xID,Kx->xKx,WLr->desWLr,ISS->desISS}


desWLr[5]=WLr[5]//.Solve[
	(xVINmin==reqVINmin)//.
	{VDD->reqVDD,VON->xVON,ID->xID,Kx->xKx,ISS->desISS,VT->xVT,WLr[2]->desWLr[2]},WLr[5]][[1]]


desVON[5]=xVON[5]//.{ID->xID,Kx->xKx,WLr->desWLr,ISS->desISS}
