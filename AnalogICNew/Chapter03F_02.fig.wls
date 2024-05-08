#!/usr/bin/env wolframscript
(* ::Package:: *)

ClearAll["Global`*"]


SetDirectory[NotebookDirectory[]];
Needs["MaTeX`"];


nbFileName="Chapter03F_02.fig.wls";
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
xT0=300;


xVt=kB T/q;
xW=1u;
xL=1u;
xT0=300;
xVT0=0.7;
xalpha=-2.3m;
xKx0=110u;


rval={kB->xkB,q->xq,m0->xm0,T0->xT0,Vt->xVt,W->xW,L->xL,T0->xT0,VT0->xVT0,alpha->xalpha,Kx0->xKx0};


xKx=Kx0 (T/T0)^(-1.5);
xVT=VT0+alpha(T-T0);
xbeta=Kx (W/L);
rfunc={Kx->xKx,VT->xVT,beta->xbeta};


xiDsat=(1/2)beta (vGS-VT)^2;
xiDoff=0;
xiD=Piecewise[{{xiDsat,vGS>=VT},{xiDoff,vGS<VT}}];


xvGSZTC=vGS//.FullSimplify@Solve[D[xiDsat//.rfunc,T]==0,{vGS}][[2]];
xvGSZTCval=xvGSZTC//.Join[rval,{T->300}]


xTCFiDsat=D[xiDsat//.rfunc,T]/xiDsat;
xTCFiD=Piecewise[{{xTCFiDsat,vGS>=VT},{Indeterminate,True}}];


xTCFKx=FullSimplify[D[xKx//.rfunc,T]/xKx];
xTCFVT=FullSimplify[D[xVT//.rfunc,T]/xVT];


valTCFKx=xTCFKx//.Join[rval,rfunc]//.{T->xT0}
valTCFVT=xTCFVT//.Join[rval,rfunc]//.{T->xT0}


listT={300,400,500};
listvGS={1,1.5,2,xvGSZTCval};


legendT={xLaTeX["T=\\SI{300}{K}",0.7],xLaTeX["T=\\SI{350}{K}",0.7],xLaTeX["T=\\SI{400}{K}",0.7]};
legendvGS={xLaTeX["v_{GS}=\\SI{1.0}{V}",0.7],xLaTeX["v_{GS}=\\SI{1.5}{V}",0.7],xLaTeX["v_{GS}=\\SI{2.0}{V}",0.7],xLaTeX["v_{GS}=V_{GS}(\\text{ZTC})",0.7]};


legendVT={xLaTeX["V_{T0}",0.7],xLaTeX["V_T",0.7],xLaTeX["V_{GS}(\\text{ZTC})",0.7]};
legendKx={xLaTeX["K_{0}'",0.7],xLaTeX["K'",0.7]};


labelTCFKx=xLaTeX["TC_{F}(K')=\\SI{-0.0050}{K^{-1}}",0.7];
labelTCFVT=xLaTeX["TC_{F}(V_T)=\\SI{-0.0033}{K^{-1}}",0.7];


labelVT0=xLaTeX["V_{T0}=\\SI{0.70}{V}",0.7];
labelKx0=xLaTeX["K'_0=\\SI{110}{uA\\cdot V^{-2}}",0.7];
labelVGSZTC=xLaTeX["V_{GS}(\\text{ZTC})=\\SI{1.62}{V}",0.7];
labelZTC=xLaTeX["\\text{ZTC}",0.7];


stylevGS=Join[Table[ColorData[10][[4]][x],{x,1,3}],{Directive[Gray,Dashed]}];
styleVT={Directive[Black,Dashed],Black,Directive[Gray,Dashed]};
styleKx={Directive[Black,Dashed],Black};


ff1=Evaluate@Table[xiD//.Join[rval,rfunc],{T,listT}];
ff1ZTCy=First[ff1]//.{vGS->xvGSZTCval};
fig1=Plot[Evaluate@ff1,{vGS,0,2.5},PlotRange->{-0.005m,0.125m},
FrameTicks->{xTicks[{0,2.5,0.5},{2,1}],xTicksProp[{0,0.12,0.02},{3,2},m]},
Epilog->Join[
	xPoint[{xvGSZTCval,ff1ZTCy}],
	{Text[labelZTC,{xvGSZTCval+0.12,ff1ZTCy-0.003m}]}
],
PlotLegends->xLegend[legendT,{"Column",1},{Left,Top}],
Evaluate@xFigLabel["v_{GS}~(\\si{V})","i_D~(\\si{mA})"],
Evaluate@xFigConfig]


ff2=Evaluate@Table[xiD//.Join[rval,rfunc],{vGS,listvGS}];
ff2ZTCy=Last[ff2]//.{T->xT0};
fig2=Plot[Evaluate@ff2,{T,250,500},PlotRange->{-0.005m,0.125m},
PlotStyle->stylevGS,
FrameTicks->{xTicks[{250,500,50},{3,0}],xTicksProp[{0,0.12,0.02},{3,2},m]},
Epilog->Join[
	xPoint[{xT0,ff2ZTCy}],
	{Text[labelZTC,{xT0+10,ff2ZTCy+0.005m}]}
],
PlotLegends->xLegend[legendvGS,{"Column",2},{Right,Top}],
Evaluate@xFigLabel["T~(\\si{K})","i_D~(\\si{mA})"],
Evaluate@xFigConfig]


ff3={xVT0,xVT//.Join[rval,rfunc],xvGSZTC//.Join[rval,rfunc]};
fig3=Plot[Evaluate@ff3,{T,250,500},PlotRange->{-0.1,2.1},
PlotStyle->styleVT,
FrameTicks->{xTicks[{250,500,50},{3,0}],xTicks[{0,2,0.5},{2,1}]},
Epilog->Join[
	xPoint[{{xT0,xVT0},{xT0,xvGSZTCval}}],
	{Text[labelVT0,{xT0+23,xVT0+0.1}],Text[labelVGSZTC,{xT0+20,xvGSZTCval+0.12}]}
],
PlotLegends->xLegend[legendVT,{"Column",2},{Left,Bottom}],
PlotStyle->{Directive[Black,Dashed],Black,Directive[Gray,Dashed]},
Evaluate@xFigLabel["T~(\\si{K})","V_T,V_{GS}(\\text{ZTC})~~(\\si{V})"],
Evaluate@xFigConfig]


ff4=Evaluate@Table[xTCFiD//.Join[rval,rfunc],{T,listT}];
ff4ZTCy=0;
fig4=Plot[Evaluate@ff4,{vGS,0,2.5},PlotRange->{-0.004,0.008},
FrameTicks->{xTicks[{0,2.5,0.5},{2,1}],xTicks[{-0.004,0.008,0.002},{4,3}]},
Epilog->Join[
	xPoint[{xvGSZTCval,ff4ZTCy}],
	{Text[labelZTC,{xvGSZTCval-0.1,ff4ZTCy-0.0005}]}
],
PlotLegends->xLegend[legendT,{"Column",1},{Left,Top}],
Evaluate@xFigLabel["v_{GS}~(\\si{V})","TC_F(i_D)~~(\\si{K^{-1}})"],
Evaluate@xFigConfig]


ff5=Evaluate@Table[xTCFiD//.Join[rval,rfunc],{vGS,listvGS}];
ff5ZTCy=0;
fig5=Plot[Evaluate@ff5,{T,250,500},PlotRange->{-0.004,0.008},
PlotStyle->stylevGS,
FrameTicks->{xTicks[{250,500,50},{3,0}],xTicks[{-0.004,0.008,0.002},{4,3}]},
Epilog->Join[
	xPoint[{xT0,ff5ZTCy}],
	{Text[labelZTC,{xT0,ff2ZTCy-0.0007}]}
],
PlotLegends->xLegend[legendvGS,{"Column",2},{Right,Bottom}],
Evaluate@xFigLabel["T~(\\si{K})","TC_F(i_D)~~(\\si{K^{-1}})"],
Evaluate@xFigConfig]


ff6={xKx0,xKx//.Join[rval,rfunc]};
fig6=Plot[Evaluate@ff6,{T,250,500},
PlotStyle->styleKx,PlotRange->{-5u,205u},
FrameTicks->{xTicks[{250,500,50},{3,0}],xTicksProp[{0,200,50},{3,0},u]},
Epilog->Join[
	xPoint[{xT0,xKx0}],
	{Text[labelKx0,{xT0+20,xKx0+10u}]}
],
PlotLegends->xLegend[legendKx,{"Column",1},{Left,Bottom}],
Evaluate@xFigLabel["T~(\\si{K})","K'~(\\si{uA\\cdot V^{-2}})"],
Evaluate@xFigConfig]


ff7=Evaluate[xTCFVT//.Join[rval,rfunc]];
fig7=Plot[Evaluate@ff7,{T,250,500},PlotRange->{-0.0105,-0.0015},
PlotStyle->Black,
FrameTicks->{xTicks[{250,500,50},{3,0}],xTicks[{-0.010,-0.002,0.002},{4,3}]},
Epilog->{xPoint[{xT0,valTCFVT}],Text[labelTCFVT,{xT0+43,valTCFVT+0.0003}]},
Evaluate@xFigLabel["T~(\\si{K})","TC_F(V_T)~~(\\si{K^{-1}})"],
Evaluate@xFigConfig]


ff8=Evaluate[xTCFKx//.Join[rval,rfunc]];
fig8=Plot[Evaluate@ff8,{T,250,500},PlotRange->{-0.0105,-0.0015},
PlotStyle->Black,
FrameTicks->{xTicks[{250,500,50},{3,0}],xTicks[{-0.010,-0.002,0.002},{4,3}]},
Epilog->{xPoint[{xT0,valTCFKx}],Text[labelTCFKx,{xT0+43,valTCFKx-0.0003}]},
Evaluate@xFigLabel["T~(\\si{K})","TC_F(K')~~(\\si{K^{-1}})"],
Evaluate@xFigConfig]


fig={fig1,fig2,fig3,fig4,fig5,fig6,fig7,fig8}


figWithName=Transpose[{fig,nbPDFNameList[[1;;Length[fig]]]}];
ReleaseHold[(Hold[Export[#[[2]],#[[1]],ImageResolution->2000]]&)/@figWithName]
