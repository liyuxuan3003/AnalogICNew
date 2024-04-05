#!/usr/bin/env wolframscript
(* ::Package:: *)

ClearAll["Global`*"]


SetDirectory[NotebookDirectory[]];
Needs["MaTeX`"];


nbFileName="Chapter03A_04.fig.wls";
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


f=10^(-15);
p=10^(-12);
n=10^(-9);
u=10^(-6);
m=10^(-3);
k=10^(3);
M=10^(6);


xkx=110u;
xVT0=0.7;
xgamma=0.4;
xphifx=0.7;
xbeta=(W/L)kx;
xCGDO=220 10^-12;
xCox=24.7 10^-4;
xW=0.8u;
xL=0.8u;
xVH=5;
xVL=0;
xVS=1;
xVT=VT0+gamma(Sqrt[phifx-VBS]-Sqrt[phifx]);
xCL=200f;
xVHST=VH-VS-VT;
xVBS=-VS;


rval={
	kx->xkx,VT0->xVT0,gamma->xgamma,phifx->xphifx,beta->xbeta,CGDO->xCGDO,Cox->xCox,
	W->xW,L->xL,VH->xVH,VL->xVL,VS->xVS,VT->xVT,CL->xCL,
	VHST->xVHST,VBS->xVBS
};


xVerrorSlow=((W CGDO+W L Cox/2)/CL)((Pi U CL)/(2 beta))^(1/2)+(W CGDO)/CL (VH-VL+VT);
xVerrorFast=((W CGDO+W L Cox/2)/CL)(VHST-(beta VHST^3)/(6 U CL))+(W CGDO)/CL (VH-VL+VT);


xtdown=(VH-VL)/U;


xUmid=(beta VHST^2)/(2 CL);


cy={0.2n,10n,xtdown//.{U->xUmid}//.rval};
cx=(U//.(Solve[xtdown==tdown,U]//.rval))[[1]]//.{tdown->cy};


cx


cxv


cyv


cxv=Join[cx,{Last[cx]}];
cyv={
	xVerrorFast//.Join[rval,{U->cxv[[1]]}],
	xVerrorSlow//.Join[rval,{U->cxv[[2]]}],
	xVerrorFast//.Join[rval,{U->cxv[[3]]}],
	xVerrorSlow//.Join[rval,{U->cxv[[4]]}]};


fig1=LogLinearPlot[Evaluate[{xVerrorSlow,xVerrorFast}//.rval],{U,10^7,10^11},PlotRange->{-2m,42m},
FrameTicks->{xTicksLog[{7,11,1}],xTicksProp[{0.00,0.04,0.005},{4,3},1]},
GridLines->{xGridsLog[{7,11,1},{10^7,10^11}],xGrids[{0.00,0.04,0.005}]},
PlotStyle->{Blue,Red},
Epilog->{
	PointSize[0.015],Black,
	Point[Inner[List,Log[cxv[[1;;4]]],cyv[[1;;4]],List]],
	PointSize[0.010],White,
	Point[Inner[List,Log[cxv[[1;;4]]],cyv[[1;;4]],List]],
	Text[xLaTeX["\\SI{0.01095}{V}",0.6],{Log[2.7 10^8],0.012}],
	Text[xLaTeX["\\SI{0.01969}{V}",0.6],{Log[4.6 10^10],0.0215}]
},
Prolog->
{
	Thickness[0.001],Black,Dashed,
	Line/@(Outer[{#1,#2}&,Log[cxv[[1;;2]]],{-0.01,0.05}]),
	Line/@(Outer[{#2,#1}&,cyv[[1;;2]],Log[{10^6,10^12}]]),
	Thickness[0.001],Red,
	Line/@(Outer[{#1,#2}&,Log[cxv[[3;;3]]],{-0.01,0.05}]),
	Line/@(Outer[{#2,#1}&,cyv[[3;;3]],Log[{10^6,10^12}]])
},
PlotLegends->Placed[LineLegend[{
	xLaTeX["v_{\\phi,slow}",0.8],
	xLaTeX["v_{\\phi,fast}",0.8]},
	LegendFunction->"Panel",LegendLayout -> {"Column",1},
	LegendMargins->1,
	LegendMarkerSize->{22,1}],{Left,Top}],
Evaluate@xFigLabel["U~(\\si{V\\cdot s^{-1}})","V_{error}~(\\si{V})"],
Evaluate@xFigConfig]


ff2=Evaluate[xtdown//.rval];
fig2=LogLogPlot[ff2,{U,10^7,10^11},PlotRange->{0.5 10^-11,2 10^-6},
FrameTicks->{xTicksLog[{7,11,1}],xTicksLog[{-11,-6,1}]},
GridLines->{xGridsLog[{7,11,1},{10^7,10^11}],xGridsLog[{-11,-6,1},{10^-11,10^-6}]},
PlotStyle->Black,
Epilog->{
	PointSize[0.015],Black,
	Point[Inner[List,Log[cx[[1;;3]]],Log[cy[[1;;3]]],List]],
	PointSize[0.010],White,
	Point[Inner[List,Log[cx],Log[cy],List]],
	Text[xLaTeX["t=\\SI{10}{ns}~~U=\\SI{5.0e8}{}",0.5],{Log[1.7 10^8],Log[0.6 10^-8]}],
	Text[xLaTeX["slow",0.6],{Log[1.7 10^8],Log[0.3 10^-8]}],
	Text[xLaTeX["t=\\SI{0.2}{ns}~~U=\\SI{2.5e10}{}",0.5],{Log[8.2 10^9],Log[1.2 10^-10]}],
	Text[xLaTeX["fast",0.6],{Log[8.2 10^9],Log[0.6 10^-10]}]
},
Prolog->
{
	Thickness[0.001],Black,Dashed,
	Line/@(Outer[{#1,#2}&,Log[cx[[1;;2]]],Log[{10^-12,10^-5}]]),
	Line/@(Outer[{#2,#1}&,Log[cy[[1;;2]]],Log[{10^6,10^12}]]),
	Thickness[0.001],Red,
	Line/@(Outer[{#1,#2}&,Log[cx[[3;;3]]],Log[{10^-12,10^-5}]]),
	Line/@(Outer[{#2,#1}&,Log[cy[[3;;3]]],Log[{10^6,10^12}]])
},
Evaluate@xFigLabel["U~(\\si{V\\cdot s^{-1}})","t_{HL}~(\\si{s})"],
Evaluate@xFigConfig]


fig={fig1,fig2}


figWithName=Transpose[{fig,nbPDFNameList[[1;;Length[fig]]]}];
ReleaseHold[(Hold[Export[#[[2]],#[[1]],ImageResolution->2000]]&)/@figWithName]
