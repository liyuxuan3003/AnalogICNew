#!/usr/bin/env wolframscript
(* ::Package:: *)

ClearAll["Global`*"]


SetDirectory[NotebookDirectory[]];
Needs["MaTeX`"];


nbFileName="Chapter03E_08.fig.wls";
nbPDFName="build/"<>StringReplace[nbFileName,"wls"->"pdf"];
nbPDFNameList=(StringReplace[nbPDFName,".fig"->#<>".fig"]&)/@Alphabet[];
Export[nbPDFName,""]


xLaTeX=(MaTeX[
	#1,Magnification->#2,
	"Preamble"->{"\\usepackage{xcolor}","\\usepackage{siunitx}","\\usepackage{physics}","\\usepackage{upgreek}"}])&;
xLabelFont={FontColor->Black,FontFamily->"Latin Modern",FontSize->5};


xFigConfig={
	Frame->True,GridLines->Automatic,ImageSize->Automatic->300,
	FrameStyle->BlackFrame,PlotStyle->ColorData[10],BaseStyle->xLabelFont};
xFigLabel=(FrameLabel->{xLaTeX[#1,0.8],xLaTeX[#2,0.8]})&;


xFig3DConfig={
	ViewPoint->{-100,-80,65},PlotPoints -> 50,MaxRecursion -> 5,
	ExclusionsStyle->{None,{Black,Thin}},ColorFunction->ColorData["DarkRainbow"],
	BaseStyle->xLabelFont};
xFig3DLabel=(AxesLabel->{xLaTeX[#1,0.8],xLaTeX[#2,0.8],xLaTeX[#3,0.8]})&;


xTicks=(Table[{x,xLaTeX[ToString[PaddedForm[x,{#2[[1]],#2[[2]]}]],0.5]},{x,#1[[1]],#1[[2]],#1[[3]]}])&;
xTicksProp=(Table[{x #3,xLaTeX[ToString[PaddedForm[x,{#2[[1]],#2[[2]]}]],0.5]},{x,#1[[1]],#1[[2]],#1[[3]]}])&;
xTicksLog=(Table[{10^x,xLaTeX["10^{"<>ToString[x]<>"}",0.5]},{x,#1[[1]],#1[[2]],#1[[3]]}])&;
xTicksLog1=(Table[{10^x,xLaTeX["1\\times 10^{"<>ToString[x]<>"}",0.5]},{x,#1[[1]],#1[[2]],#1[[3]]}])&;
xTicksLog2=(Table[{2 10^x,xLaTeX["2\\times 10^{"<>ToString[x]<>"}",0.5]},{x,#1[[1]],#1[[2]],#1[[3]]}])&;
xTicksLog5=(Table[{5 10^x,xLaTeX["5\\times 10^{"<>ToString[x]<>"}",0.5]},{x,#1[[1]],#1[[2]],#1[[3]]}])&;


xGrids=(Table[x,{x,#1[[1]],#1[[2]],#1[[3]]}])&;
xGridsLog=(Flatten[Table[Table[{If[#2[[1]]<=(n 10^x)<=#2[[2]],(n 10^x),{}]},{n,1,9,1}],{x,#1[[1]],#1[[2]],#1[[3]]}]])&;


xText=Style[#1,FontSize->#2,FontFamily->"Noto Serif CJK SC"]&;
xTextConsolas=Style[#1,FontSize->#2,FontFamily->"Consolas"]&;


xPoint=({Black,PointSize[0.015],Point[#],White,PointSize[0.010],Point[#],Black})&;


f=10^(-15);
p=10^(-12);
n=10^(-9);
u=10^(-6);
m=10^(-3);
k=10^(3);
M=10^(6);


xkB=1.38064852 10^(-23);
xq=1.602176634 10^-19;
xm0=9.109 383 56*10^-31;
xT=300;


xVt=kB T/q;
xKx=110 10^(-6);
xVT=0.7;


rval={kB->xkB,q->xq,m0->xm0,T->xT,Vt->xVt,Kx->xKx,VT->xVT};


xImos=Piecewise[{{Kx WLr (VREF-VT)^2/2,VREF>VT},{0,VREF<=VT}}];
xImossat=Kx WLr (VREF-VT)^2/2;
xIr0=(VREF)/R0;


xVON=0.8;
xWLr=2;


xVREF=VT+1/(Kx WLr R0) (Sqrt[2 Kx WLr R0 VT+1]+1);


xVREF=FullSimplify[VREF//.Solve[xImossat==xIr0,{VREF}][[2]]]


FindRoot[(xVREF-VT==VON)//.Join[rval,{WLr->xWLr,VON->xVON}],{R0,100k}]


xR0=100k;


valVREF=VREF//.(Solve[xImossat==xIr0,{VREF}][[2]])//.Join[rval,{WLr->xWLr,R0->xR0}]
valI0=xImos//.Join[rval,{WLr->xWLr,R0->xR0,VREF->xVREF}]


latexVREF=xLaTeX["V_{REF}=\\SI{1.0018}{V}",0.7];


ff3=Evaluate[{xImos,xIr0}//.Join[rval,{WLr->xWLr,R0->xR0}]];
fig3=Plot[ff3,{VREF,0,2},PlotRange->{-0.005m,0.105m},PlotStyle->{Red,Blue},
FrameTicks->{xTicks[{0,2,0.2},{2,1}],xTicksProp[{0,0.1,0.02},{3,2},m]},
GridLines->{xGrids[{0,2,0.1}],xGrids[{0m,0.10m,0.01m}]},
Epilog->Join[
	{Black,Text[latexVREF,{valVREF+0.23,valI0-0.007m}]},
	xPoint[{{valVREF,valI0},{0,0}}]
	],
PlotLegends->Placed[LineLegend[{
		xLaTeX["I_{0,M1}",0.7],
		xLaTeX["I_{0,R0}",0.7]
	},
	LegendMarkerSize->{{22,1}},
	LegendFunction->"Panel",
	LegendMargins->1,
	LegendLayout -> {"Column",1}],{Left,Top}],
Evaluate@xFigLabel["V_{REF}~(\\si{V})","I_0~(\\si{mA})"],
Evaluate@xFigConfig]


fig={fig3}


figWithName=Transpose[{fig,nbPDFNameList[[1;;Length[fig]]]}];
ReleaseHold[(Hold[Export[#[[2]],#[[1]],ImageResolution->2000]]&)/@figWithName]
