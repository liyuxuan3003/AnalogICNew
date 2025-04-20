#!/usr/bin/env wolframscript
(* ::Package:: *)

ClearAll["Global`*"]


eq1=s Cm(vout-vinx)+s Cout vout+Rout^-1(vout-Av vinx)==0;
eq2=s Cm(vinx-vout)+s Cin vinx+Rs^-1(vinx-vin)==0;


sol=FullSimplify@Solve[{eq1,eq2},{vout,vinx}]


xAvs=FullSimplify[(vout//.sol[[1]])/vin]
