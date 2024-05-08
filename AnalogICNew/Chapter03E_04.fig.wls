#!/usr/bin/env wolframscript
(* ::Package:: *)

ClearAll["Global`*"]


SetDirectory[NotebookDirectory[]];
Needs["MaTeX`"];


nbFileName="Chapter03E_04.fig.wls";
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


rval={kB->xkB,q->xq,m0->xm0,T->xT,Vt->xVt};


xIr1=VREF/R1;
xIr0=(VDD-VREF)/R0;


xVREF=VREF//.Solve[xIr1==xIr0,{VREF}][[1]];
xVREFappr=Vt Log[VDD/(R0 Is)];


xSVV=D[xVREF,VDD](VREF/VDD)^(-1);
xSVVappr=D[xVREFappr,VDD](VREF/VDD)^(-1);


xR0=100k;
xR1=25k;
xVDD=5;


xVDDlist={4.5,5.0,5.5};


style={Black};


valVREF=xVREF//.Join[rval,{R1->xR1,R0->xR0,VDD->xVDD}]
valSVV=xSVV//.Join[rval,{R1->xR1,R0->xR0,VDD->xVDD,VREF->xVREF}]


latexVREF=xLaTeX["V_{REF}=\\SI{1.0000}{V}",0.7];
latexVREFWithVDD=xLaTeX["V_{REF}=\\SI{1.0000}{V}|_{V_{DD}=\\SI{5.0}{V}}",0.7];
latexSVV=xLaTeX["S_{V_{DD}}^{V_{REF}}=1.0000",0.7];


ff1=Evaluate[{xVREF}//.Join[rval,{R1->xR1,R0->xR0}]];
fig1=Plot[ff1,{VDD,0,8},PlotRange->{-0.1,2.1},PlotStyle->style,
FrameTicks->{xTicks[{0,8,1},{2,1}],xTicks[{0,2,0.2},{2,1}]},
GridLines->{Automatic,xGrids[{0,2,0.1}]},
Epilog->Join[
	{Blue,Thickness[0.004],Dotted,{Line[{{xVDD,-10},{xVDD,10}}]}},
	{Black,Text[latexVREF,{xVDD+1,valVREF-0.05}]},
	xPoint[{xVDD,valVREF}]],
PlotLegends->Placed[LineLegend[{
		xLaTeX["V_{REF}",0.7]
	},
	LegendMarkerSize->{{22,1}},
	LegendFunction->"Panel",
	LegendMargins->1,
	LegendLayout -> {"ReversedColumn",1}],{Left,Top}],
Evaluate@xFigLabel["V_{DD}~(\\si{V})","V_{REF}~(\\si{V})"],
Evaluate@xFigConfig]


ff2=Evaluate[{xVREFappr,xVREF}//.Join[rval,{VDD->xVDD,R0->xR0}]];
fig2=LogLinearPlot[ff2,{R1,1k,1M},PlotRange->{-0.1,2.1},PlotStyle->style,
FrameTicks->{xTicksLog[{3,6,1}],xTicks[{0,2,0.2},{2,1}]},
GridLines->{xGridsLog[{3,6,1},{10^3,10^6}],xGrids[{0,2,0.1}]},
Epilog->Join[
	{Blue,Thickness[0.004],Dotted,{Line[{{Log[xR1],-10},{Log[xR1],10}}]}},
	{Black,Text[latexVREF,{Log[xR1 2.5],valVREF-0.05}]},
	xPoint[{Log[xR1],valVREF}]],
PlotLegends->Placed[LineLegend[{
		xLaTeX["V_{REF}",0.7]
	},
	LegendMarkerSize->{{22,1}},
	LegendFunction->"Panel",
	LegendMargins->1,
	LegendLayout -> {"ReversedColumn",1}],{Right,Top}],
Evaluate@xFigLabel["R_{1}~(\\si{\\ohm})","V_{REF}~(\\si{V})"],
Evaluate@xFigConfig]


valVREFlist=VREF//.Table[FindRoot[(xIr1==xIr0)//.Join[rval,{R1->xR1,R0->xR0}],{VREF,valVREF}],{VDD,xVDDlist}];
valI0list=Table[xIr1//.Join[rval,{R1->xR1,R0->xR0}],{VREF,valVREFlist}];


ff3PointList=Inner[List,valVREFlist,valI0list,List];


ff3=Evaluate[Join[{xIr1},Table[xIr0,{VDD,xVDDlist}]]//.Join[rval,{R1->xR1,R0->xR0}]];
fig3=Plot[ff3,{VREF,0,2},PlotRange->{-0.005m,0.105m},PlotStyle->{Red,Lighter[Lighter[Blue]],Blue,Darker[Darker[Blue]]},
FrameTicks->{xTicks[{0,2,0.2},{2,1}],xTicksProp[{0,0.1,0.02},{3,2},m]},
GridLines->{xGrids[{0,2,0.1}],xGrids[{0m,0.10m,0.01m}]},
Epilog->Join[
	{Black,Text[latexVREFWithVDD,{valVREFlist[[2]]+0.19,valI0list[[2]]-0.015m}]},
	xPoint[ff3PointList]],
PlotLegends->Placed[LineLegend[{
		xLaTeX["I_{R1}",0.7],
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


ff4=Evaluate[{xSVV}//.Join[rval,{Is->xIs,R0->xR0,VREF->xVREF}]];
fig4=Plot[ff4,{VDD,0,8},PlotRange->{-0.05,1.05},
PlotStyle->style,
FrameTicks->{xTicks[{0,8,1},{2,1}],xTicks[{0,1,0.1},{2,1}]},
GridLines->{Automatic,xGrids[{0,1,0.1}]},
Epilog->Join[
	{Blue,Thickness[0.004],Dotted,{Line[{{xVDD,-10},{xVDD,10}}]}},
	{Black,Text[latexSVV,{xVDD+1,valSVV-0.05}]},
	xPoint[{xVDD,valSVV}]],
PlotLegends->Placed[LineLegend[{
		xLaTeX["S_{V_{DD}}^{V_{REF}}",0.7]
	},
	LegendMarkerSize->{{22,1}},
	LegendFunction->"Panel",
	LegendMargins->1,
	LegendLayout -> {"ReversedColumn",1}],{Right,Bottom}],
Evaluate@xFigLabel["V_{DD}~(\\si{V})","S_{V_{DD}}^{V_{REF}}"],
Evaluate@xFigConfig]


ff5=Evaluate[{xSVV}//.Join[rval,{VDD->xVDD,R0->xR0,VREF->xVREF}]];
fig5=LogLinearPlot[ff5,{R1,1k,1M},PlotRange->{-0.05,1.05},PlotStyle->style,
FrameTicks->{xTicksLog[{3,6,1}],xTicks[{0,1,0.1},{2,1}]},
GridLines->{xGridsLog[{3,6,1},{10^3,10^6}],xGrids[{0,1,0.1}]},
Epilog->Join[
	{Blue,Thickness[0.004],Dotted,{Line[{{Log[xR1],-10},{Log[xR1],10}}]}},
	{Black,Text[latexSVV,{Log[xR1 2.5],valSVV-0.05}]},
	xPoint[{Log[xR1],valSVV}]],
PlotLegends->Placed[LineLegend[{
		xLaTeX["S_{V_{DD}}^{V_{REF}}",0.7]
	},
	LegendMarkerSize->{{22,1}},
	LegendFunction->"Panel",
	LegendMargins->1,
	LegendLayout -> {"ReversedColumn",1}],{Right,Bottom}],
Evaluate@xFigLabel["R_1~(\\si{\\ohm})","S_{V_{DD}}^{V_{REF}}"],
Evaluate@xFigConfig]


fig={fig1,fig2,fig3,fig4,fig5}


figWithName=Transpose[{fig,nbPDFNameList[[1;;Length[fig]]]}];
ReleaseHold[(Hold[Export[#[[2]],#[[1]],ImageResolution->2000]]&)/@figWithName]
