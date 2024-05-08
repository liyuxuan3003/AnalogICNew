#!/usr/bin/env wolframscript
(* ::Package:: *)

ClearAll["Global`*"]


SetDirectory[NotebookDirectory[]];
Needs["MaTeX`"];


nbFileName="Chapter03F_03.fig.wls";
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
xVG0=1.205;
xIs0[1]=10^-15;
xIs0[2]=10^-14;
xID=0.1m;


rval={kB->xkB,q->xq,m0->xm0,T0->xT0,Vt0->xVt0,
gamma->xgamma,VG0->xVG0,ID->xID,Is0->xIs0
};


xVt=kB T/q;
xVD=(Vt Log[ID/Is[#]])&;
xIs=(Is0[#] (T/T0)^gamma Exp[VG0/Vt0-VG0/Vt])&;


rfunc={Vt->xVt,VD->xVD,Is->xIs};


xVDval[1]=xVD[1]//.Join[rval,rfunc,{T->300}]
xVDval[2]=xVD[2]//.Join[rval,rfunc,{T->300}]
xVPTATval=xVDval[1]-xVDval[2]


xVDvalLaTeX[1]=xLaTeX["V_{D1}=\\SI{0.6548}{V}",0.7];
xVDvalLaTeX[2]=xLaTeX["V_{D2}=\\SI{0.5953}{V}",0.7];
xVPTATLaTeX=xLaTeX["V_{PTAT}=\\SI{0.0595}{V}",0.7];


ff1={xVD[1],xVD[2]}//.Join[rval,rfunc];
fig1=Plot[Evaluate@ff1,{T,250,500},PlotRange->{-0.02,0.82},
PlotStyle->{Black,Gray},
FrameTicks->{xTicks[{250,500,50},{3,0}],xTicks[{0,0.8,0.1},{2,1}]},
PlotLegends->xLegend[
	{xLaTeX["V_{D1}",0.7],xLaTeX["V_{D2}",0.7]},
	{"Column",1},
	{Right,Top}
],
Epilog->{xPoint[{{xT0,xVDval[1]},{xT0,xVDval[2]}}],
	Text[xVDvalLaTeX[1],{xT0+30,xVDval[1]+0.017}],
	Text[xVDvalLaTeX[2],{xT0-12,xVDval[2]-0.058}]
},
Evaluate@xFigLabel["T~(\\si{K})","V_{D1},V_{D2}~(\\si{V})"],
Evaluate@xFigConfig]


((xVD[1]-xVD[2])//.Join[rval,rfunc]//.{T->300})/xID


ff2={xVD[1]-xVD[2]}//.Join[rval,rfunc];
fig2=Plot[Evaluate@ff2,{T,250,500},
PlotStyle->{Red},PlotRange->{-0.005,0.105},
FrameTicks->{xTicks[{250,500,50},{3,0}],xTicks[{0,0.1,0.02},{3,2}]},
Epilog->{xPoint[{xT0,xVPTATval}],Text[xVPTATLaTeX,{xT0+30,xVPTATval-0.005}]},
Evaluate@xFigLabel["T~(\\si{K})","V_{PTAT}~(\\si{V})"],
Evaluate@xFigConfig]


fig={fig1,fig2}


figWithName=Transpose[{fig,nbPDFNameList[[1;;Length[fig]]]}];
ReleaseHold[(Hold[Export[#[[2]],#[[1]],ImageResolution->2000]]&)/@figWithName]
