#!/usr/bin/env wolframscript
(* ::Package:: *)

ClearAll["Global`*"]


SetDirectory[NotebookDirectory[]];
Needs["MaTeX`"];


nbFileName="Chapter03E_05.fig.wls";
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
xalpha=-2.3m;
xTCFR0=1500u;


xVt=kB T/q;
xKx=110 10^(-6);
xVT=0.7;


xVTtemp=VT+alpha(T-T0);
xKxtemp=Kx (T/T0)^(-3/2);
xR0temp=R0 Exp[(T-T0)TCFR0];


D[xVTtemp,T]
D[xKxtemp,T]
D[xR0temp,T]


DSolve[y'[x]==a y[x],y,x]


rval={kB->xkB,q->xq,m0->xm0,T0->xT0,Vt->xVt,Kx->xKx,VT->xVT,alpha->xalpha,TCFR0->xTCFR0};


xImos=Piecewise[{{Kx WLr (VREF-VT)^2/2,VREF>VT},{0,VREF<=VT}}];
xImossat=Kx WLr (VREF-VT)^2/2;
xIr0=(VDD-VREF)/R0;


xVREFsat=FullSimplify[VREF//.Solve[xImossat==xIr0,{VREF}][[2]]];
xVREF=Piecewise[{{xVREFsat,VDD>VT},{VDD,VDD<=VT}}];


xVREFsat


xdVREFGiven=(alpha+(Sqrt[(VDD-VREF)/(2 Kx WLr R0)])(1.5/T-TCFR0))/(1+1/Sqrt[2 Kx WLr R0(VDD-VREF)]);


xTCF=(1/VREF)(alpha-(VREF-VT)(-(3/2)(1/T)+TCFR0)+(VDD-VT)/(Kx WLr R0 (VREF-VT)+1) (-(3/2)(1/T)+TCFR0)-alpha/(Kx WLr R0 (VREF-VT)+1));


xTCFs1=(1/VREF)(alpha-1/(Kx^2 WLr^2 R0^2) (-(3/2)Kx WLr(1/T) R0+Kx WLr R0 TCFR0)(Sqrt[1+2 Kx WLr R0(VDD-VT)]-1)+1/(2 Kx WLr R0) 2/Sqrt[1+2 Kx WLr R0(VDD-VT)]((-(3/2)Kx WLr(1/T) R0+Kx WLr R0 TCFR0)(VDD-VT)-alpha Kx WLr R0));


xTCFs2=(1/VREF)(alpha-1/(Kx WLr R0) (-(3/2)(1/T)+TCFR0)(Sqrt[1+2 Kx WLr R0(VDD-VT)]-1)+1/Sqrt[1+2 Kx WLr R0(VDD-VT)] ((-(3/2)(1/T)+TCFR0)(VDD-VT)-alpha));


xTCFGiven=Simplify[((1/VREF)xdVREFGiven)/.{VREF->xVREFsat}/.{VT->xVTtemp,Kx->xKxtemp,R0->xR0temp}];


xTCFDtMMAPart1=xVREFsat/.{VT->xVTtemp,Kx->xKxtemp,R0->xR0temp};
xTCFDtMMA=(1/xTCFDtMMAPart1)D[xTCFDtMMAPart1,T];


xTCFDtHAND=Simplify[xTCF/.{VREF->xVREFsat}/.{VT->xVTtemp,Kx->xKxtemp,R0->xR0temp}];
xTCFDtHANDs1=Simplify[xTCFs1/.{VREF->xVREFsat}/.{VT->xVTtemp,Kx->xKxtemp,R0->xR0temp}];
xTCFDtHANDs2=Simplify[xTCFs2/.{VREF->xVREFsat}/.{VT->xVTtemp,Kx->xKxtemp,R0->xR0temp}];


xSVV=D[xVREF,VDD](VREF/VDD)^(-1);


xWLr=2;
xR0=100k;
xVDD=5;


xVDDlist={4.5,5.0,5.5};


style={Black};


valVREF=xVREF//.Join[rval,{WLr->xWLr,R0->xR0,VDD->xVDD}]
valSVV=xSVV//.Join[rval,{WLr->xWLr,R0->xR0,VDD->xVDD,VREF->xVREF}]


valTCF=xTCFDtMMA/u//.Join[rval,{WLr->xWLr,R0->xR0,VDD->xVDD,T->xT0}]


latexVREF=xLaTeX["V_{REF}=\\SI{1.2814}{V}",0.7];
latexVREFWithVDD=xLaTeX["V_{REF}=\\SI{1.2814}{V}|_{V_{DD}=\\SI{5.0}{V}}",0.7];
latexSVV=xLaTeX["S_{V_{DD}}^{V_{REF}}=0.2829",0.7];
latexTCF=xLaTeX["TC_F(V_{REF})=\\SI{-928}{ppm\\cdot K^{-1}}",0.7];


ff1=Evaluate[{xVREF}//.Join[rval,{WLr->xWLr,R0->xR0,T->T0}]];
fig1=Plot[ff1,{VDD,0,8},PlotRange->{-0.1,2.1},PlotStyle->style,
FrameTicks->{xTicks[{0,8,1},{2,1}],xTicks[{0,2,0.2},{2,1}]},
GridLines->{Automatic,xGrids[{0,2,0.1}]},
Epilog->Join[
	{Blue,Thickness[0.004],Dotted,{Line[{{xVDD,-10},{xVDD,10}}]}},
	{Black,Text[latexVREF,{xVDD+1,valVREF+0.2}]},
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


ff2=Evaluate[{xVREFappr,xVREF}//.Join[rval,{VDD->xVDD,R0->xR0,T->T0}]];
fig2=LogLinearPlot[ff2,{WLr,1,1000},PlotRange->{-0.1,2.1},PlotStyle->style,
FrameTicks->{xTicksLog[{0,3,1}],xTicks[{0,2,0.2},{2,1}]},
GridLines->{xGridsLog[{0,3,1},{1,1000}],xGrids[{0,2,0.1}]},
Epilog->Join[
	{Blue,Thickness[0.004],Dotted,{Line[{{Log[xWLr],-10},{Log[xWLr],10}}]}},
	{Black,Text[latexVREF,{Log[xWLr 2.5],valVREF+0.2}]},
	xPoint[{Log[xWLr],valVREF}]],
PlotLegends->Placed[LineLegend[{
		xLaTeX["V_{REF}",0.7]
	},
	LegendMarkerSize->{{22,1}},
	LegendFunction->"Panel",
	LegendMargins->1,
	LegendLayout -> {"ReversedColumn",1}],{Right,Top}],
Evaluate@xFigLabel["(W/L)","V_{REF}~(\\si{V})"],
Evaluate@xFigConfig]


valVREFlist=VREF//.Table[FindRoot[(xImos==xIr0)//.Join[rval,{WLr->xWLr,R0->xR0}],{VREF,valVREF}],{VDD,xVDDlist}];
valI0list=Table[xImos//.Join[rval,{WLr->xWLr,R0->xR0}],{VREF,valVREFlist}];


valVREFlist


valI0list


ff3PointList=Inner[List,valVREFlist,valI0list,List];


ff3=Evaluate[Join[{xImos},Table[xIr0,{VDD,xVDDlist}]]//.Join[rval,{WLr->xWLr,R0->xR0,T->T0}]];
fig3=Plot[ff3,{VREF,0,2},PlotRange->{-0.005m,0.105m},PlotStyle->{Red,Lighter[Lighter[Blue]],Blue,Darker[Darker[Blue]]},
FrameTicks->{xTicks[{0,2,0.2},{2,1}],xTicksProp[{0,0.1,0.02},{3,2},m]},
GridLines->{xGrids[{0,2,0.1}],xGrids[{0m,0.10m,0.01m}]},
Epilog->Join[
	{Black,Text[latexVREFWithVDD,{valVREFlist[[2]]+0.3,valI0list[[2]]-0.014m}]},
	xPoint[ff3PointList]],
PlotLegends->Placed[LineLegend[{
		xLaTeX["I_{M1}",0.7],
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


FullSimplify[xSVV//.{VREF->xVREF},Assumptions->{VDD>VT}]


ff4=Evaluate[{xSVV}//.Join[rval,{WLr->xWLr,R0->xR0,VREF->xVREF,T->T0}]];
fig4=Plot[ff4,{VDD,0,8},PlotRange->{-0.05,1.05},MaxRecursion->10,
PlotStyle->style,
FrameTicks->{xTicks[{0,8,1},{2,1}],xTicks[{0,1,0.1},{2,1}]},
GridLines->{Automatic,xGrids[{0,1,0.1}]},
Epilog->Join[
	{Blue,Thickness[0.004],Dotted,{Line[{{xVDD,-10},{xVDD,10}}]}},
	{Black,Text[latexSVV,{xVDD+1,valSVV+0.07}]},
	xPoint[{xVDD,valSVV}]],
PlotLegends->Placed[LineLegend[{
		xLaTeX["S_{V_{DD}}^{V_{REF}}",0.7]
	},
	LegendMarkerSize->{{22,1}},
	LegendFunction->"Panel",
	LegendMargins->1,
	LegendLayout -> {"ReversedColumn",1}],{Right,Top}],
Evaluate@xFigLabel["V_{DD}~(\\si{V})","S_{V_{DD}}^{V_{REF}}"],
Evaluate@xFigConfig]


ff5=Evaluate[{xSVV}//.Join[rval,{VDD->xVDD,R0->xR0,VREF->xVREF,T->T0}]];
fig5=LogLinearPlot[ff5,{WLr,1,1000},PlotRange->{-0.05,1.05},PlotStyle->style,
FrameTicks->{xTicksLog[{0,3,1}],xTicks[{0,1,0.1},{2,1}]},
GridLines->{xGridsLog[{0,3,1},{1,1000}],xGrids[{0,1,0.1}]},
Epilog->Join[
	{Blue,Thickness[0.004],Dotted,{Line[{{Log[xWLr],-10},{Log[xWLr],10}}]}},
	{Black,Text[latexSVV,{Log[xWLr 2.5],valSVV+0.07}]},
	xPoint[{Log[xWLr],valSVV}]],
PlotLegends->Placed[LineLegend[{
		xLaTeX["S_{V_{DD}}^{V_{REF}}",0.7]
	},
	LegendMarkerSize->{{22,1}},
	LegendFunction->"Panel",
	LegendMargins->1,
	LegendLayout -> {"ReversedColumn",1}],{Right,Top}],
Evaluate@xFigLabel["(W/L)","S_{V_{DD}}^{V_{REF}}"],
Evaluate@xFigConfig]


ff6=Evaluate[({xTCFDtHAND}/u)//.Join[rval,{WLr->xWLr,VDD->xVDD,R0->xR0}]];
fig6=Plot[Evaluate[ff6],{T,250,500},PlotRange->{-1650,-750},
FrameTicks->{xTicks[{250,500,50},{3,0}],xTicks[{-1600,-800,200},{4,0}]},
PlotStyle->Black,
Epilog->Join[
	xPoint[{xT0,valTCF}],
	{Text[latexTCF,{xT0+52,valTCF+20}]}
],
Evaluate@xFigLabel["T~(\\si{K})","TC_F(V_{REF})~~(\\si{ppm\\cdot K^{-1}})"],
Evaluate@xFigConfig]


ff7=Evaluate[xVREFsat/.{Kx->xKxtemp,R0->xR0temp,VT->xVTtemp}//.Join[rval,{VDD->xVDD,WLr->xWLr,R0->xR0}]];
fig7=Plot[ff7,{T,250,500},PlotRange->{0.98,1.37},
FrameTicks->{xTicks[{250,500,50},{3,0}],xTicks[{1.0,1.35,0.05},{3,2}]},
PlotStyle->Black,
Epilog->Join[
	xPoint[{xT0,valVREF}],
	{Text[latexVREF,{xT0+32,valVREF+0.005}]}
],
Evaluate@xFigLabel["T~(\\si{K})","V_{REF}~(\\si{V})"],
Evaluate@xFigConfig]


fig={fig1,fig2,fig3,fig4,fig5,fig6,fig7}


figWithName=Transpose[{fig,nbPDFNameList[[1;;Length[fig]]]}];
ReleaseHold[(Hold[Export[#[[2]],#[[1]],ImageResolution->2000]]&)/@figWithName]
