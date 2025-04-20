#!/usr/bin/env wolframscript
(* ::Package:: *)

ClearAll["Global`*"]


Solve[s^2/wn^2-2zeta s/wn+1==0,s]


FullSimplify@Integrate[Exp[alpha tau]Sin[k tau],{tau,0,t}]


Plot[{ArcTan[(2zeta)/Sqrt[Sqrt[4zeta^4+1]-2zeta^2]],ArcCos[Sqrt[4zeta^4+1]-2zeta^2]},{zeta,-2,2},Frame->True]


Plot[{2 Sin[t]-3Cos[t],Sqrt[2^2+3^2]Sin[t-ArcTan[3/2]]},{t,-5,5}]


N[ArcTan[2]]


N[ArcTan[-2]]


Arg[1/(-3 I-2)]
