#!/usr/bin/env wolframscript
(* ::Package:: *)

ClearAll["Global`*"]


Solve[s^2/wn^2-2zeta s/wn+1==0,s]


FullSimplify@Integrate[Exp[alpha tau]Sin[k tau],{tau,0,t}]


Plot[{ArcTan[(2zeta)/Sqrt[Sqrt[4zeta^4+1]-2zeta^2]],ArcCos[Sqrt[4zeta^4+1]-2zeta^2]},{zeta,-2,2},Frame->True]
