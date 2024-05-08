#!/usr/bin/env wolframscript
(* ::Package:: *)

ClearAll["Global`*"]


SetDirectory[NotebookDirectory[]];
Needs["MaTeX`"];


nbFileName="Chapter03F_04.fig.wls";
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


xVt0=kB T0/q;
xT0=300;
xVT0=0.7;
xgamma=3;
xalpha=1;
xVG0=1.205;
xVD0=0.6;
xIs0[1]=10^-15;
xIs0[2]=10^-14;


rval={kB->xkB,q->xq,m0->xm0,T0->xT0,Vt0->xVt0,
gamma->xgamma,alpha->xalpha,VG0->xVG0,VD0->xVD0,Is0->xIs0
};


xVt=kB T/q;
xVPTAT=Vt Log[Is0[2]/Is0[1]];
xVCTAT=VG0-(VG0-VD0)(T/T0)-Vt(gamma-alpha)Log[T/T0];
xVCTATalt=VG0-(VG0-VD0)(T/T0);


xVPTAT0=xVPTAT//.{T->T0,Vt->Vt0};
xVCTAT0=xVCTAT//.{T->T0,Vt->Vt0};


xK=-((VCTAT0-VG0-Vt0(gamma-alpha))/VPTAT0);
xVREF=K VPTAT+VCTAT;


rfunc={Vt->xVt,VCTAT->xVCTAT,VPTAT->xVPTAT,VCTAT0->xVCTAT0,VPTAT0->xVPTAT0,K->xK,VREF->xVREF};


xTCF=(1/VREF)D[{xVREF}//.Join[rval,rfunc],T];


xDVPTAT=D[(xVPTAT//.Join[rval,rfunc]),T];
xDVCTAT=D[(xVCTAT//.Join[rval,rfunc]),T];
xKT=-xDVCTAT/xDVPTAT;


xVREF0


xVPTAT0val=xVPTAT0//.Join[rval,rfunc]
xVCTAT0val=xVCTAT0//.Join[rval,rfunc]
xVREF0val=xVREF//.Join[rval,rfunc,{T->300}]
xKVal=xK//.Join[rval,rfunc,{T->300}]


xVPTAT0LaTeX=xLaTeX["V_{PTAT}=\\SI{0.0595}{V}",0.7];
xVCTAT0LaTeX=xLaTeX["V_{CTAT}=\\SI{0.6000}{V}",0.7];
xVREF0LaTeX=xLaTeX["V_{REF}=\\SI{1.2567}{V}",0.7];
xVG0LaTeX=xLaTeX["V_{G0}=\\SI{1.2050}{V}",0.7];
xKLaTeX=xLaTeX["K=\\SI{11.0321}{}",0.7];
xTCFLaTeX=xLaTeX["TC_F(V_{REF})=\\SI{0}{ppm\\cdot K^{-1}}",0.7];


ff1={xVPTAT,xVCTAT,xVCTATalt}//.Join[rval,rfunc];
fig1=Plot[Evaluate@ff1,{T,250,500},PlotRange->{-0.02,0.82},
PlotStyle->{Red,Blue,Directive[Blue,Thickness[0.003],Dashed]},
FrameTicks->{xTicks[{250,500,50},{3,0}],xTicks[{0,0.8,0.1},{2,1}]},
Epilog->{
	xPoint[{{xT0,xVPTAT0val},{xT0,xVCTAT0val}}],
	Text[xVPTAT0LaTeX,{xT0+10,xVPTAT0val+0.04}],
	Text[xVCTAT0LaTeX,{xT0+30,xVCTAT0val+0.03}]
},
PlotLegends->xLegend[
	{xLaTeX["V_{PTAT}",0.7],xLaTeX["V_{CTAT}",0.7],xLaTeX["V_{CTAT}~(\\text{Linear})",0.7]},
	{"Column",1},
	{Right,Top}
],
Evaluate@xFigLabel["T~(\\si{K})","V_{PTAT},V_{CTAT}~(\\si{V})"],
Evaluate@xFigConfig]


ff2={xVREF,xVG0}//.Join[rval,rfunc];
fig2=Plot[Evaluate@ff2,{T,250,500},PlotRange->{1.198,1.262},
PlotStyle->{Black,Gray},
FrameTicks->{xTicks[{250,500,50},{3,0}],xTicks[{1.20,1.28,0.01},{3,2}]},
GridLines->{Automatic,xGrids[{1.20,1.26,0.005}]},
Epilog->{
	xPoint[{{xT0,xVG0},{xT0,xVREF0val}}],
	Text[xVREF0LaTeX,{xT0+20,xVREF0val-0.004}],
	Text[xVG0LaTeX,{xT0+5,xVG0+0.0028}]
},
PlotLegends->xLegend[
	{xLaTeX["V_{REF}",0.7],xLaTeX["V_{G0}",0.7]},{"Column",1},{Right,Bottom}
],
Evaluate@xFigLabel["T~(\\si{K})","V_{REF},V_{G0}~(\\si{V})"],
Evaluate@xFigConfig]


ff3=(1/u)xTCF//.Join[rval,rfunc];
fig3=Plot[Evaluate@ff3,{T,250,500},
PlotStyle->Black,
PlotRange->{-85,45},
FrameTicks->{xTicks[{250,500,50},{3,0}],xTicks[{-80,40,20},{2,0}]},
Epilog->{xPoint[{xT0,0}],Text[xTCFLaTeX,{xT0+43,5}]},
Evaluate@xFigLabel["T~(\\si{K})","TC_F(V_{REF})~~(\\si{ppm\\cdot K^{-1}})"],
Evaluate@xFigConfig]


ff4={xKT}//.Join[rval,rfunc];
fig4=Plot[Evaluate@ff4,{T,250,500},PlotRange->{10.75,11.55},
PlotStyle->Black,
FrameTicks->{xTicks[{250,500,50},{3,0}],xTicks[{10.8,11.5,0.1},{3,1}]},
Epilog->{xPoint[{xT0,xKVal}],Text[xKLaTeX,{xT0+24,xKVal-0.013}]},
Evaluate@xFigLabel["T~(\\si{K})","K(T)"],
Evaluate@xFigConfig]


fig={fig1,fig2,fig3,fig4}


figWithName=Transpose[{fig,nbPDFNameList[[1;;Length[fig]]]}];
ReleaseHold[(Hold[Export[#[[2]],#[[1]],ImageResolution->2000]]&)/@figWithName]
