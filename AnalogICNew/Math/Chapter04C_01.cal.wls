#!/usr/bin/env wolframscript
(* ::Package:: *)

ClearAll["Global`*"]


fmlvd1=gm1 vgs1+gds1 vds1-gm2 vgs2-gbs2 vbs2-gds2 vds2;
fmlvout=gm2 vgs2+gbs2 vbs2+gds2 vds2+gds3 vds3;


rv={
	vgs1->vin,vgs2->-vd1,vgs3->0,
	vbs1->0,vbs2->-vd1,vbs3->0,
	vds1->vd1,vds2->vout-vd1,vds3->vout
};


eqAv={fmlvd1==0,fmlvout==0}//.rv;
eqRout={fmlvd1==0,fmlvout==iout}//.rv//.{vin->0};
eqRd1={fmlvd1==id1,fmlvout==0}//.rv//.{vin->0};
eqAvs={fmlvd1+s Cd1 vd1+s Cind1(vd1-vin)==0,fmlvout+s Cout vout==0}//.rv;


sAv=Solve[eqAv,{vd1,vout}];
sRout=Solve[eqRout,{vd1,vout}];
sRd1=Solve[eqRd1,{vd1,vout}];
sAvs=Solve[eqAvs,{vd1,vout}];


xAv=FullSimplify[(vout//.sAv[[1]])/vin]
xAd1=FullSimplify[(vd1//.sAv[[1]])/vin]
xRout=FullSimplify[(vout//.sRout[[1]])/iout]
xRd1=FullSimplify[(vd1//.sRd1[[1]])/id1]


xAvs=FullSimplify[(vout//.sAvs[[1]])/vin]


FullSimplify[xAvs/xAv]
