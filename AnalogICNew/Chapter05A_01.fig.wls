#!/usr/bin/env wolframscript
(* ::Package:: *)

ClearAll["Global`*"]


SetDirectory[NotebookDirectory[]];
Needs["MaTeX`"];


nbFileName="Chapter05A_01.fig.wls";
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


xLegend=(Placed[LineLegend[#1,
	LegendMarkerSize->{{22,1}},
	LegendFunction->"Panel",
	LegendMargins->1,
	LegendLayout->#2],#3])&;


xPoint=({Black,PointSize[0.015],Point[#],White,PointSize[0.010],Point[#],Black})&;


f=10^(-15);
p=10^(-12);
n=10^(-9);
u=10^(-6);
m=10^(-3);
k=10^(3);
M=10^(6);


MixColorTable=(Table[Blend[{#1,#2},x],{x,0+#4/2,1-#4/2,(1-#4)/(#3-1)}])&;


xAvfs=Avf/(1-2zeta s/wn+s^2/wn^2);


xAvfsmax=1/(2zeta Sqrt[1-zeta]);


ff1=Table[FullSimplify[Abs[(xAvfs/Avf)//.{s->I w}//.{w->x wn}]],{zeta,{0.1,0.2,0.5,Sqrt[2]/2,1,2}}];
fig1=LogLogPlot[Evaluate@ff1,{x,0.1,10},
PlotRange->{0.8 10^-2,1.2 10^1},
PlotStyle->Table[Directive[Thickness[0.003],x],{x,MixColorTable[Blue,Red,6,0.2]}],
FrameTicks->{xTicksLog[{-1,1,1}],xTicksLog[{-2,1,1}]},
GridLines->{xGridsLog[{-1,1,1},{10^-1,10^1}],xGridsLog[{-2,1,1},{10^-2,10^1}]},
PlotLegends->xLegend[{
	xLaTeX["\\zeta=0.10",0.7],
	xLaTeX["\\zeta=0.20",0.7],
	xLaTeX["\\zeta=0.50",0.7],
	xLaTeX["\\zeta=0.707",0.7],
	xLaTeX["\\zeta=1.00",0.7],
	xLaTeX["\\zeta=2.00",0.7]
},{"ReversedColumn",2},{Left,Bottom}],
Evaluate@xFigLabel["\\omega/\\omega_n","|A_{vf}(s)/A_{vf}|"],
Evaluate@xFigConfig]


(* ::InheritFromParent:: *)
(**)


ff2=Piecewise[{{xAvfsmax,zeta<1},{100,zeta>1}}];
fig2=LogPlot[ff2,{zeta,0,1},
	PlotRange->{0.8 10^-2,1.2 10^1},PlotStyle->Black,
	Prolog->{Thickness[0.003],Red,Dashed,Line[{{Sqrt[2]/2,Log[100]},{Sqrt[2]/2,Log[0.001]}}]},
	FrameTicks->{xTicks[{0,1,0.1},{2,1}],xTicksLog[{-2,1,1}]},
	GridLines->{Automatic,xGridsLog[{-2,1,1},{10^-2,10^1}]},
	Evaluate@xFigLabel["\\zeta","\\max |A_{vf}(s)/A_{vf}|"],
	Evaluate@xFigConfig
]


xvoutDelta=FullSimplify@InverseLaplaceTransform[xAvfs,s,t]


xvoutEpsi=FullSimplify@Integrate[xvoutDelta//.{t->tau},{tau,0,t}]


ff3=FullSimplify@Table[(xvoutEpsi/Avf)//.{t->-x/wn},{zeta,{2,0.999999,Sqrt[2]/2,0.5,0.2,0.1}}];
fig3=Plot[ff3,{x,0,25},PlotRange->{-0.05,2.05},
PlotStyle->Table[Directive[Thickness[0.003],x],{x,MixColorTable[Red,Blue,6,0.2]}],
FrameTicks->{xTicks[{0,25,5},{2,0}],xTicks[{0,2,0.5},{2,1}]},
PlotLegends->xLegend[{
	xLaTeX["\\zeta=2.00",0.7],
	xLaTeX["\\zeta=1.00",0.7],
	xLaTeX["\\zeta=0.707",0.7],
	xLaTeX["\\zeta=0.50",0.7],
	xLaTeX["\\zeta=0.20",0.7],
	xLaTeX["\\zeta=0.10",0.7]
},{"Column",2},{Right,Top}],
Evaluate@xFigLabel["\\omega_n t","v_{OUT}(t)/A_{vf}"],
Evaluate@xFigConfig]


xtp=-(Pi/(wn Sqrt[1-zeta^2]));


xvoutOver=Exp[(-Pi zeta)/Sqrt[1-zeta^2]];


ff4=FullSimplify[(-xtp wn)];
fig4=Plot[ff4,{zeta,0,1},PlotRange->{0,25},PlotStyle->Black,
FrameTicks->{xTicks[{0,1,0.1},{2,1}],xTicks[{0,25,5},{2,0}]},
Evaluate@xFigLabel["\\zeta","\\omega_n t_p"],
Evaluate@xFigConfig]


ff5=FullSimplify[xvoutOver];
fig5=Plot[ff5,{zeta,0,1},PlotRange->{-0.05,1.05},PlotStyle->Black,
FrameTicks->{xTicks[{0,1,0.1},{2,1}],xTicks[{0,1,0.2},{2,1}]},
Evaluate@xFigLabel["\\zeta","v_{OUT}(t_p)/A_{vf}-1"],
Evaluate@xFigConfig]


ff6=FullSimplify[xvoutOver];
fig6=LogPlot[ff6,{zeta,0,1},PlotRange->{0.5 10^-4,2},PlotStyle->Black,
FrameTicks->{xTicks[{0,1,0.1},{2,1}],xTicksLog[{-4,0,1}]},
GridLines->{Automatic,xGridsLog[{-4,0,1},{10^-4,10^0}]},
Evaluate@xFigLabel["\\zeta","v_{OUT}(t_p)/A_{vf}-1"],
Evaluate@xFigConfig]


xwc=wn Sqrt[Sqrt[4zeta^4+1]-2zeta^2];


ff7=FullSimplify[xwc/wn];
fig7=Plot[ff7,{zeta,0,1},PlotRange->{-0.05,1.05},PlotStyle->Black,
FrameTicks->{xTicks[{0,1,0.1},{2,1}],xTicks[{0,1,0.2},{2,1}]},
Evaluate@xFigLabel["\\zeta","\\omega_c/\\omega_n"],
Evaluate@xFigConfig]


xphim=ArcCos[Sqrt[4 zeta^4+1]-2 zeta^2]


ff8=FullSimplify[xphim/Pi];
fig8=Plot[ff8,{zeta,0,1},PlotRange->{-0.025,0.525},PlotStyle->Black,
FrameTicks->{xTicks[{0,1,0.1},{2,1}],xTicks[{0,0.5,0.1},{2,1}]},
Epilog->{Red,Thickness[0.004],Dashed,Line[{{-1,1/3},{2,1/3}}],Line[{{-1,1/4},{2,1/4}}],
Text[xLaTeX["\\pi/3",0.7],{0.95,1/3+0.02}],
Text[xLaTeX["\\pi/4",0.7],{0.95,1/4+0.02}]
},
Evaluate@xFigLabel["\\zeta","\\phi_m/\\pi"],
Evaluate@xFigConfig]


sol=s//.Solve[1+2zeta s/wn+s^2/wn^2==0,s]
wp1x=sol[[1]];
wp2x=sol[[2]];


xzeta=-(1/2)((w1+w2)/Sqrt[w1 w2 (1+Avf0 Fx)]);
xwn=Sqrt[w1 w2(1+Avf0 Fx)];


FullSimplify[{wp1x,wp2x}//.{zeta->xzeta,wn->xwn}]


fig={fig1,fig2,fig3,fig4,fig5,fig6,fig7,fig8}


figWithName=Transpose[{fig,nbPDFNameList[[1;;Length[fig]]]}];
ReleaseHold[(Hold[Export[#[[2]],#[[1]],ImageResolution->2000]]&)/@figWithName]
