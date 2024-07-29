#!/usr/bin/env wolframscript
(* ::Package:: *)

eq1=gm1 vin1+(gm3+gds1+gds3)vout1==0;
eq2=gm2 vin2+(gds2+gds4)vout2+gm4 vout1==0;


s=Solve[({eq1,eq2}//.{vin1->+vid/2,vin2->-vid/2}),{vout2,vout1}]


FullSimplify[(vout2//.s[[1]])/vid]
