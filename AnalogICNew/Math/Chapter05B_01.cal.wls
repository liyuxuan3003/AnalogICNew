#!/usr/bin/env wolframscript
(* ::Package:: *)

ClearAll["Global`*"]


eq1=Gm1 vin+s C1 vc+R1^-1 vc+s Cc(vc-vout)==0;
eq2=Gm2 vc+s C2 vout+R2^-1 vout+s Cc(vout-vc)==0;


sol=Solve[{eq1,eq2},{vc,vout}]
