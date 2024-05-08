#!/usr/bin/env wolframscript
(* ::Package:: *)

ClearAll["Global`*"]


SetDirectory[NotebookDirectory[]];
Needs["MaTeX`"];


nbFileName="Chapter03E_06.fig.wls";
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
xT0=300;
xVG0=1.206;
xTCFR0=1500u;


xVt=kB T/q;


rval={kB->xkB,q->xq,m0->xm0,T0->xT0,Vt->xVt,VG0->xVG0,TCFR0->xTCFR0};


xIbjt=Is (Exp[VREF/Vt]-1);
xIbjtappr=Is Exp[VREF/Vt];
xIr0=(VDD-VREF)/R0;
xIr0appr=VDD/R0;


xTCF=1/(VREF T) (VREF-VG0-3Vt)-Vt/VREF TCFR0;


xR0temp=R0 Exp[(T-T0)TCFR0];


xIs=Is0 (T^3 Exp[-q VG0/(kB T)])/(T0^3 Exp[-q VG0/(kB T0)]);


xVREF=VREF//.Solve[xIbjt==xIr0,{VREF}][[1]];
xVREFappr=Vt Log[VDD/(R0 Is)];


xVREF


xSVV=D[xVREF,VDD](VREF/VDD)^(-1);
xSVVappr=D[xVREFappr,VDD](VREF/VDD)^(-1);


xVREFappr


xTCFMMA=(1/xVREFappr)D[xVREFappr/.{Is->xIs,Vt->xVt,R0->xR0temp},T];


xIs0=10^(-15);
xR0=100k;
xVDD=5;


xVDDlist={4.5,5.0,5.5};


style={Directive[Gray,Dashed,Thickness[0.003]],Black};


valVREF=xVREF//.Join[rval,{Is->xIs0,R0->xR0,VDD->xVDD,T->T0}]
valSVV=xSVV//.Join[rval,{Is->xIs0,R0->xR0,VDD->xVDD,VREF->xVREF,T->T0}]
valTCF=xTCF/u//.Join[rval,{Is->xIs0,R0->xR0,VDD->xVDD,VREF->xVREF,T->T0}]


latexVREF=xLaTeX["V_{REF}=\\SI{0.6334}{V}",0.7];
latexVREFWithVDD=xLaTeX["V_{REF}=\\SI{0.6334}{V}|_{V_{DD}=\\SI{5.0}{V}}",0.7];
latexSVV=xLaTeX["S_{V_{DD}}^{V_{REF}}=0.0465",0.7];
latexTCF=xLaTeX["TC_F(V_{REF})=\\SI{-3483}{ppm\\cdot K^{-1}}",0.7];


ff1=Evaluate[{xVREFappr,xVREF}//.Join[rval,{Is->xIs0,R0->xR0,T->T0}]];
fig1=Plot[ff1,{VDD,0,8},PlotRange->{-0.05,1.05},PlotStyle->style,
FrameTicks->{xTicks[{0,8,1},{2,1}],xTicks[{0,1,0.2},{2,1}]},
GridLines->{Automatic,xGrids[{0,1,0.1}]},
Epilog->Join[
	{Blue,Thickness[0.004],Dotted,{Line[{{xVDD,-10},{xVDD,10}}]}},
	{Black,Text[latexVREF,{xVDD+1,valVREF+0.05}]},
	xPoint[{xVDD,valVREF}]],
PlotLegends->Placed[LineLegend[{
		xLaTeX["V_{REF}~\\approx",0.7],
		xLaTeX["V_{REF}",0.7]
	},
	LegendMarkerSize->{{22,1}},
	LegendFunction->"Panel",
	LegendMargins->1,
	LegendLayout -> {"ReversedColumn",1}],{Left,Top}],
Evaluate@xFigLabel["V_{DD}~(\\si{V})","V_{REF}~(\\si{V})"],
Evaluate@xFigConfig]


ff2=Evaluate[{xVREFappr,xVREF}//.Join[rval,{VDD->xVDD,R0->xR0,T->T0}]];
fig2=LogLinearPlot[ff2,{Is,10^(-18),10^(-12)},PlotRange->{-0.05,1.05},PlotStyle->style,
FrameTicks->{xTicksLog[{-18,-12,1}],xTicks[{0,1,0.2},{2,1}]},
GridLines->{xGridsLog[{-18,-12,1},{10^-18,10^-12}],xGrids[{0,1,0.1}]},
Epilog->Join[
	{Blue,Thickness[0.004],Dotted,{Line[{{Log[xIs0],-10},{Log[xIs0],10}}]}},
	{Black,Text[latexVREF,{Log[xIs0 5.5],valVREF+0.05}]},
	xPoint[{Log[xIs0],valVREF}]],
PlotLegends->Placed[LineLegend[{
		xLaTeX["V_{REF}~\\approx",0.7],
		xLaTeX["V_{REF}",0.7]
	},
	LegendMarkerSize->{{22,1}},
	LegendFunction->"Panel",
	LegendMargins->1,
	LegendLayout -> {"ReversedColumn",1}],{Right,Top}],
Evaluate@xFigLabel["I_{s}~(\\si{A})","V_{REF}~(\\si{V})"],
Evaluate@xFigConfig]


valVREFlist=VREF//.Table[FindRoot[(xIbjt==xIr0)//.Join[rval,{Is->xIs0,R0->xR0,T->T0}],{VREF,valVREF}],{VDD,xVDDlist}];
valI0list=Table[xIbjt//.Join[rval,{Is->xIs0,R0->xR0,T->T0}],{VREF,valVREFlist}];


ff3PointList=Inner[List,valVREFlist,valI0list,List];


ff3=Evaluate[Join[{xIbjt},Table[xIr0,{VDD,xVDDlist}]]//.Join[rval,{Is->xIs0,R0->xR0,T->T0}]];
fig3=Plot[ff3,{VREF,0,1},PlotRange->{-0.005m,0.105m},PlotStyle->{Red,Lighter[Lighter[Blue]],Blue,Darker[Darker[Blue]]},
FrameTicks->{xTicks[{0,1,0.2},{2,1}],xTicksProp[{0,0.1,0.02},{3,2},m]},
GridLines->{xGrids[{0,1,0.1}],xGrids[{0m,0.10m,0.01m}]},
Epilog->Join[
	{Black,Text[latexVREFWithVDD,{valVREFlist[[2]]+0.19,valI0list[[2]]+0.01m}]},
	xPoint[ff3PointList]],
PlotLegends->Placed[LineLegend[{
		xLaTeX["I_{Q1}",0.7],
		xLaTeX["I_{R0}~V_{DD}=\\SI{4.5}{V}",0.7],
		xLaTeX["I_{R0}~V_{DD}=\\SI{5.0}{V}",0.7],
		xLaTeX["I_{R0}~V_{DD}=\\SI{5.5}{V}",0.7]
	},
	LegendMarkerSize->{{22,1}},
	LegendFunction->"Panel",
	LegendMargins->1,
	LegendLayout -> {"Column",1}],{Left,Top}],
Evaluate@xFigLabel["V_{REF}~(\\si{V})","I_0~(\\si{mA})"],
Evaluate@xFigConfig]


ff4=Evaluate[{xSVVappr,xSVV}//.Join[rval,{Is->xIs0,R0->xR0,VREF->xVREF,T->T0}]];
fig4=Plot[ff4,{VDD,0,8},PlotRange->{-0.05,1.05},
PlotStyle->style,
FrameTicks->{xTicks[{0,8,1},{2,1}],xTicks[{0,1,0.1},{2,1}]},
GridLines->{Automatic,xGrids[{0,1,0.1}]},
Epilog->Join[
	{Blue,Thickness[0.004],Dotted,{Line[{{xVDD,-10},{xVDD,10}}]}},
	{Black,Text[latexSVV,{xVDD+1,valSVV+0.07}]},
	xPoint[{xVDD,valSVV}]],
PlotLegends->Placed[LineLegend[{
		xLaTeX["S_{V_{DD}}^{V_{REF}}~\\approx",0.7],
		xLaTeX["S_{V_{DD}}^{V_{REF}}",0.7]
	},
	LegendMarkerSize->{{22,1}},
	LegendFunction->"Panel",
	LegendMargins->1,
	LegendLayout -> {"ReversedColumn",1}],{Right,Top}],
Evaluate@xFigLabel["V_{DD}~(\\si{V})","S_{V_{DD}}^{V_{REF}}"],
Evaluate@xFigConfig]


ff5=Evaluate[{xSVVappr,xSVV}//.Join[rval,{VDD->xVDD,R0->xR0,VREF->xVREF,T->T0}]];
fig5=LogLinearPlot[ff5,{Is,10^(-18),10^(-12)},PlotRange->{-0.05,1.05},PlotStyle->style,
FrameTicks->{xTicksLog[{-18,-12,1}],xTicks[{0,1,0.1},{2,1}]},
GridLines->{xGridsLog[{-18,-12,1},{10^-18,10^-12}],xGrids[{0,1,0.1}]},
Epilog->Join[
	{Blue,Thickness[0.004],Dotted,{Line[{{Log[xIs0],-10},{Log[xIs0],10}}]}},
	{Black,Text[latexSVV,{Log[xIs0 5.5],valSVV+0.07}]},
	xPoint[{Log[xIs0],valSVV}]],
PlotLegends->Placed[LineLegend[{
		xLaTeX["S_{V_{DD}}^{V_{REF}}~\\approx",0.7],
		xLaTeX["S_{V_{DD}}^{V_{REF}}",0.7]
	},
	LegendMarkerSize->{{22,1}},
	LegendFunction->"Panel",
	LegendMargins->1,
	LegendLayout -> {"ReversedColumn",1}],{Right,Top}],
Evaluate@xFigLabel["I_{s}~(\\si{A})","S_{V_{DD}}^{V_{REF}}"],
Evaluate@xFigConfig]


ff6=Evaluate[{xTCF}/.{R0->xR0temp}//.Join[rval,{VDD->xVDD,R0->xR0,VREF->xVREF,Is->xIs,Is0->xIs0}]];
fig6=Plot[ff6/u,{T,250,500},PlotRange->{-14500,-1500},
FrameTicks->{xTicks[{250,500,50},{3,0}],xTicks[{-14000,-2000,2000},{4,0}]},
PlotStyle->Black,
Epilog->Join[
	xPoint[{xT0,valTCF}],
	{Text[latexTCF,{xT0+52,valTCF+500}]}
],
Evaluate@xFigLabel["T~(\\si{K})","TC_F(V_{REF})~~(\\si{ppm\\cdot K^{-1}})"],
Evaluate@xFigConfig]


ff7=Evaluate[{xVREF}/.{R0->xR0temp}//.Join[rval,{VDD->xVDD,R0->xR0,VREF->xVREF,Is->xIs,Is0->xIs0}]];
fig7=Plot[ff7,{T,250,500},PlotRange->{-0.02,0.82},
FrameTicks->{xTicks[{250,500,50},{3,0}],xTicks[{0.0,0.8,0.1},{2,1}]},
PlotStyle->Black,
Epilog->Join[
	xPoint[{xT0,valVREF}],
	{Text[latexVREF,{xT0+32,valVREF+0.01}]}
],
Evaluate@xFigLabel["T~(\\si{K})","V_{REF}~(\\si{V})"],
Evaluate@xFigConfig]


fig={fig1,fig2,fig3,fig4,fig5,fig6,fig7}


figWithName=Transpose[{fig,nbPDFNameList[[1;;Length[fig]]]}];
ReleaseHold[(Hold[Export[#[[2]],#[[1]],ImageResolution->2000]]&)/@figWithName]
