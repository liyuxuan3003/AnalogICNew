#!/usr/bin/env wolframscript
(* ::Package:: *)

ClearAll["Global`*"]


SetDirectory[NotebookDirectory[]];
Needs["MaTeX`"];


nbFileName="Chapter03F_01.fig.wls";
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


xLegend=(Placed[LineLegend[#1,
	LegendMarkerSize->{{22,1}},
	LegendFunction->"Panel",
	LegendMargins->1,
	LegendLayout->#2],#3])&;


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
xVt0=kB T0/q;
xT0=300;
xgamma=3;


xIs0=10^(-15);
xVG0=1.206;


rval={kB->xkB,q->xq,m0->xm0,T0->xT0,Vt0->xVt0,Is0->xIs0,VG0->xVG0,gamma->xgamma};


xVt=kB T/q;


xIs=Is0 (T/T0)^gamma Exp[VG0/Vt0-VG0/Vt];


rfunc={Is->xIs,Vt->xVt};


xIsTCF=(1/xIs)D[xIs//.rfunc,T];


xiD=Is Exp[vD/Vt];
xiDTCF=(1/xiD)D[xiD//.rfunc,T];


listT={300,350,400};
listvD={0.6,0.5,0.4};


valvD=0.6;
valT=300;


legendT={xLaTeX["T=\\SI{300}{K}",0.7],xLaTeX["T=\\SI{350}{K}",0.7],xLaTeX["T=\\SI{400}{K}",0.7]};
legendvD={xLaTeX["v_D=\\SI{0.6}{V}",0.7],xLaTeX["v_D=\\SI{0.5}{V}",0.7],xLaTeX["v_D=\\SI{0.4}{V}",0.7]};
legendIS={xLaTeX["I_s",0.7],xLaTeX["I_{s0}",0.7]};


valiDTCF=xiDTCF//.Join[rval,rfunc]//.{T->valT,vD->valvD}
valIsTCF=xIsTCF//.Join[rval,rfunc]//.{T->valT,vD->valvD}


labelTCFiD=xLaTeX["TC_F(i_D)=\\SI{0.0881}{K^{-1}}",0.7];
labelTCFIs=xLaTeX["TC_F(i_D)=\\SI{0.1655}{K^{-1}}",0.7];


labelIs0=xLaTeX["I_{s0}=\\SI{e-15}{A}",0.7];


ff1=Table[xiD//.Join[rfunc,rval],{T,listT}];
fig1=Plot[Evaluate@ff1,{vD,0,1},PlotRange->{-0.05m,0.55m},
FrameTicks->{xTicks[{0,1,0.1},{2,1}],xTicksProp[{0,0.5,0.1},{2,1},m]},
PlotLegends->xLegend[legendT,{"Column",1},{Left,Top}],
Evaluate@xFigLabel["v_D~(\\si{V})","i_D~(\\si{mA})"],
Evaluate@xFigConfig]


ff2=Evaluate@Table[xiD//.Join[rfunc,rval],{vD,listvD}];
fig2=Plot[Evaluate@ff2,{T,250,500},PlotRange->{-0.05m,0.55m},
FrameTicks->{xTicks[{250,500,50},{3,0}],xTicksProp[{0,0.5,0.1},{2,1},m]},
PlotLegends->xLegend[legendvD,{"Column",1},{Left,Top}],
Evaluate@xFigLabel["T~(\\si{K})","i_D~(\\si{mA})"],
Evaluate@xFigConfig]


ff3=Table[xiDTCF//.Join[rfunc,rval],{T,listT}];
fig3=Plot[Evaluate@ff3,{vD,0,1},PlotRange->{-0.01,0.26},
FrameTicks->{xTicks[{0,1,0.1},{2,1}],xTicksProp[{0,0.25,0.05},{3,2},1]},
PlotLegends->xLegend[legendT,{"Column",1},{Right,Top}],
Epilog->Join[
	xPoint[{valvD,valiDTCF}],
	{Text[labelTCFiD,{valvD+0.05,valiDTCF+0.02}]}
],
Evaluate@xFigLabel["v_D~(\\si{V})","TC_F(i_D)~~(\\si{K^{-1}})"],
Evaluate@xFigConfig]


ff4=Table[xiDTCF//.Join[rfunc,rval],{vD,listvD}];
ff3PntText=xLaTeX["TC_F(i_D)=\\SI{0.0881}{K^{-1}}",0.7];
ff3PntY=First[ff3]//.{vD->valvD}
fig4=Plot[Evaluate@ff4,{T,250,500},PlotRange->{-0.01,0.26},
FrameTicks->{xTicks[{250,500,50},{3,0}],xTicksProp[{0,0.25,0.05},{3,2},1]},
PlotLegends->xLegend[legendvD,{"Column",1},{Right,Top}],
Epilog->Join[
	xPoint[{valT,valiDTCF}],
	{Text[labelTCFiD,{valT-13,valiDTCF-0.022}]}
],
Evaluate@xFigLabel["T~(\\si{K})","TC_F(i_D)~~(\\si{K^{-1}})"],
Evaluate@xFigConfig]


ff5=Evaluate[{xIs//.Join[rfunc,rval],xIs0}];
fig5=LogPlot[Evaluate@ff5,{T,250,500},PlotRange->{0.5 10^-21,2 10^-6},
PlotStyle->{Black,Directive[Black,Dashed]},
FrameTicks->{xTicks[{250,500,50},{3,0}],xTicksLog[{-21,-6,3}]},
GridLines->{Automatic,xGridsLog[{-21,-6,1},{10^-21,10^-6}]},
PlotLegends->xLegend[legendIS,{"Column",1},{Left,Top}],
Epilog->Join[
	xPoint[{valT,Log[xIs0]}],
	{Text[labelIs0,{valT-22,Log[xIs0 5]}]}
],
Evaluate@xFigLabel["T~(\\si{K})","I_s~(\\si{A})"],
Evaluate@xFigConfig]


ff6=Evaluate[xIsTCF//.Join[rfunc,rval]];
fig6=Plot[Evaluate@ff6,{T,250,500},PlotRange->{-0.01,0.26},
PlotStyle->Black,
FrameTicks->{xTicks[{250,500,50},{3,0}],xTicksProp[{0,0.25,0.05},{3,2},1]},
Epilog->Join[
	xPoint[{valT,valIsTCF}],
	{Text[labelTCFIs,{valT+30,valIsTCF+0.015}]}
],
Evaluate@xFigLabel["T~(\\si{K})","TC_F(I_s)~~(\\si{K^{-1}})"],
Evaluate@xFigConfig]


Log[1.1655,2]
Log[1.0881,2]


fig={fig1,fig2,fig3,fig4,fig5,fig6}


figWithName=Transpose[{fig,nbPDFNameList[[1;;Length[fig]]]}];
ReleaseHold[(Hold[Export[#[[2]],#[[1]],ImageResolution->2000]]&)/@figWithName]
