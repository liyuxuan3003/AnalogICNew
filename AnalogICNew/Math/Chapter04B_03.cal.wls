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
reqfc=100k;
reqCL=5p;
reqSR=10u^(-1);
reqVICmin=1.0;
reqVICmax=4.5;


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
MT[6]="n";


xKx=(vKx[MT[#]])&;
xlambda=(vlambda[MT[#]])&;
xVT=(vVT[MT[#]])&;


autoval={VDD->reqVDD,Pdiss->reqPdiss,SR->reqSR,CL->reqCL,Av->reqAv,VICmin->reqVICmin,VICmax->reqVICmax,fc->reqfc,
	Kx->xKx,lambda->xlambda,VT->xVT,
	ID->vID,VON->vVON,S->vS,gm->vgm,gds->vgds
};


xID5byVDD=Pdiss VDD^(-1)//.autoval//N


xID5byVDD=SR CL//.autoval//N


xID5byVDD=4Pi CL fc/(lambda[2]+lambda[4])//.autoval//N


vID[5]=100u//N


vID[2]=ID[5]/2//.autoval//N
vID[4]=ID[5]/2//.autoval//N
vID[6]=ID[5]//.autoval//N


vgds[2]=lambda[2] ID[2]//.autoval//N
vgds[4]=lambda[4] ID[4]//.autoval//N
vgds[5]=lambda[5] ID[5]//.autoval//N
vgds[6]=lambda[6] ID[6]//.autoval//N


vVON[4]=VDD-VICmax//.autoval//N
vS[4]=2ID[4]/(Kx[4]VON[4]^2)//.autoval//N
vgm[4]=2ID[4]/VON[4]//.autoval//N


vgm[2]=Av(gds[2]+gds[4])//.autoval//N
vS[2]=gm[2]^2/(2ID[2]Kx[2])//.autoval//N
vVON[2]=2ID[2]/gm[2]//.autoval//N


vVON[5]=VICmin-VON[2]-VT[2]//.autoval//N
vS[5]=2ID[5]/(Kx[5]VON[5]^2)//.autoval//N
vgm[5]=2ID[5]/VON[5]//.autoval//N


vVON[6]=VON[5]//.autoval//N
vS[6]=S[5]//.autoval//N
vgm[6]=gm[5]//.autoval//N



