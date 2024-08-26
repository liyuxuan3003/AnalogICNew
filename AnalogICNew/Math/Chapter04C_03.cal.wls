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
MT[4]="n";
MT[5]="p";


xKx=(vKx[MT[#]])&;
xlambda=(vlambda[MT[#]])&;
xVT=(vVT[MT[#]])&;


autoval={VDD->reqVDD,Pdiss->reqPdiss,SR->reqSR,CL->reqCL,Av->reqAv,VOUTmin->reqVOUTmin,VOUTmax->reqVOUTmax,
	Kx->xKx,lambda->xlambda,VT->xVT,
	ID->vID,VON->vVON,S->vS,gm->vgm,gds->vgds
};


xID1byVDD=Pdiss VDD^(-1)//.autoval//N


xID1byVDD=SR CL//.autoval//N


vID[1]=100u//N


vID[2]=ID[1]//.autoval//N
vID[3]=ID[1]//.autoval//N
vID[4]=ID[1]//.autoval//N
vID[5]=ID[1]//.autoval//N


vgds[1]=lambda[1] ID[1]//.autoval//N
vgds[2]=lambda[2] ID[2]//.autoval//N
vgds[3]=lambda[3] ID[3]//.autoval//N
vgds[4]=lambda[4] ID[4]//.autoval//N
vgds[5]=lambda[5] ID[5]//.autoval//N


vVON[3]=VDD-VOUTmax//.autoval//N
vS[3]=2ID[3]/(Kx[3]VON[3]^2)//.autoval//N
vgm[3]=2ID[3]/VON[3]//.autoval//N


vgm[1]=-Av gds[3]//.autoval//N
vS[1]=gm[1]^2/(2ID[1]Kx[1])//.autoval//N
vVON[1]=2ID[1]/gm[1]//.autoval//N


vVON[2]=VOUTmin-VON[1]//.autoval//N
vS[2]=2ID[2]/(Kx[2]VON[2]^2)//.autoval//N
vgm[2]=2ID[2]/VON[2]//.autoval//N


vVON[5]=VON[3]//.autoval//N
vS[5]=S[3]//.autoval//N
vgm[5]=gm[3]//.autoval//N


vVON[4]=VON[1]+VON[2]//.autoval//N
vS[4]=2ID[4]/(Kx[4]VON[4]^2)//.autoval//N
vgm[4]=2ID[4]/VON[4]//.autoval//N



