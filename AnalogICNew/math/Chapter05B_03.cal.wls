#!/usr/bin/env wolframscript
(* ::Package:: *)

ClearAll["Global`*"]


FullSimplify[(s Cc Rz^(-1))/(s Cc+Rz^(-1))]


eq1=Gm1 vin+s C1 vc+R1^-1 vc+(s Cc)(1+s Cc Rz)^(-1)(vc-vout)==0;
eq2=Gm2 vc+s C2 vout+R2^-1 vout+(s Cc)(1+s Cc Rz)^(-1)(vout-vc)==0;


sol=FullSimplify@Solve[{eq1,eq2},{vc,vout}]
