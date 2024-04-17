#!/usr/bin/env wolframscript
(* ::Package:: *)

ClearAll["Global`*"]


SetDirectory[NotebookDirectory[]];
Needs["MaTeX`"];


nbFileName="Chapter03D_05.fig.wls";
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


xKx=110 10^(-6);
xVT=0.7;
xlambda=0.04;


xperr=(1+dKx/(2Kx))/(1-dKx/(2Kx)) (1-dVT/(2(vGS-VT)))^2/(1+dVT/(2(vGS-VT)))^2-1;


xperraprr=dKx/Kx-(2dVT)/(vGS-VT);


xperrvd=(1+lambda vDS2)/(1+lambda vDS1)-1;


xvG=Sqrt[2 iD/Kx]+VT;


xiDlist={1u,10u,100u};
xvGlist=Table[xvG,{iD,xiDlist}]//.{VT->xVT,Kx->xKx};


xvGlist


ffK1={xperraprr}//.{Kx->xKx,VT->xVT,dVT->0};
fig1=Plot[Evaluate@ffK1,{dKx,-5 10^-6,5 10^-6},PlotRange->{-0.17,0.17},PlotStyle->Black,
FrameTicks->{xTicksProp[{-5,5,1},{2,1},u],xTicks[{-0.14999,0.15001,0.05},{3,2}]},
Evaluate@xFigLabel["\\Delta K'~(\\si{uA\\cdot V^{-1}})","(i_{OUT}/i_{REF})-1"],
Evaluate@xFigConfig]


ffV1=Table[{xperraprr}//.{Kx->xKx,VT->xVT,dKx->0},{vGS,xvGlist}];
fig2=Plot[Evaluate@ffV1,{dVT,-0.01,0.01},PlotRange->{-0.17,0.17},
FrameTicks->{xTicksProp[{-10,10,2},{2,0},m],xTicks[{-0.14999,0.15001,0.05},{3,2}]},
PlotLegends->Placed[LineLegend[{
	xLaTeX["I_{REF}=\\SI{1}{uA}",0.7],
	xLaTeX["I_{REF}=\\SI{10}{uA}",0.7],
	xLaTeX["I_{REF}=\\SI{100}{uA}",0.7]},
	LegendFunction->"Panel",LegendLayout -> {"Column",1},
	LegendMargins->0.3,
	LegendMarkerSize->{22,1}],{Left,Bottom}],
Evaluate@xFigLabel["\\Delta V_T~(\\si{mV})","(i_{OUT}/i_{REF})-1"],
Evaluate@xFigConfig]


ffD1=Table[{xperrvd}//.{lambda->xlambda,vDS2->(vDS1+dvDS)},{vDS1,xvGlist}];
fig3=Plot[Evaluate@ffD1,{dvDS,-xVT,4},PlotRange->{-0.17,0.17},
FrameTicks->{Join[xTicks[{0,5,1},{2,1}],xTicks[{-0.7,-0.7,1},{2,1}]],xTicks[{-0.14999,0.15001,0.05},{3,2}]},
PlotLegends->Placed[LineLegend[{
	xLaTeX["I_{REF}=\\SI{1}{uA}",0.7],
	xLaTeX["I_{REF}=\\SI{10}{uA}",0.7],
	xLaTeX["I_{REF}=\\SI{100}{uA}",0.7]},
	LegendFunction->"Panel",LegendLayout -> {"Column",1},
	LegendMargins->0.3,
	LegendMarkerSize->{22,1}],{Right,Bottom}],
Evaluate@xFigLabel["v_{DS2}-v_{DS1}~(\\si{V})","(i_{OUT}/I_{REF})-1"],
Evaluate@xFigConfig]


fig={fig1,fig2,fig3}


figWithName=Transpose[{fig,nbPDFNameList[[1;;Length[fig]]]}];
ReleaseHold[(Hold[Export[#[[2]],#[[1]],ImageResolution->2000]]&)/@figWithName]
