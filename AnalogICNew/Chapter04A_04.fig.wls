#!/usr/bin/env wolframscript
(* ::Package:: *)

ClearAll["Global`*"]


SetDirectory[NotebookDirectory[]];
Needs["MaTeX`"];


nbFileName="Chapter04A_04.fig.wls";
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


xKxN=110u;
xKxP=50u;
xlambdaN=0.04;
xlambdaP=0.05;


rval={kB->xkB,q->xq,m0->xm0,T0->xT0,
KxN->xKxN,KxP->xKxP,lambdaN->xlambdaN,lambdaP->xlambdaP};


xKx[1]=KxN;
xKx[2]=KxP;
xL[1]=1u;
xL[2]=1u;
xW[1]=1u;
xW[2]=2u;
xID=0.1m;
xlambda[1]=lambdaN;
xlambda[2]=lambdaP;
xbeta=(Kx[#] WLr[#])&;
xWLr=(W[#]/L[#])&;
xgm=(Sqrt[2 ID beta[#]])&;
xgds=(lambda[#] ID)&;


rfunc={Kx->xKx,L->xL,W->xW,lambda->xlambda,beta->xbeta,gm->xgm,gds->xgds};


revalWLr={WLr->xWLr};
revalId={ID->xID};


xAvDIO=-gm[1](gm[2]+gds[1]+gds[2])^(-1);
xAvDIOapr=-gm[1]gm[2]^(-1);
xAvCUR=-gm[1](gds[1]+gds[2])^(-1);
xAvINV=-(gm[1]+gm[2])(gds[1]+gds[2])^(-1);


xRoutDIO=(gm[2]+gds[1]+gds[2])^(-1);
xRoutDIOapr=gm[2]^(-1);
xRoutCUR=(gds[1]+gds[2])^(-1);
xRoutINV=(gds[1]+gds[2])^(-1);


valAv={xAvDIOapr,xAvDIO,xAvCUR,xAvINV}//.Join[rval,rfunc,revalWLr,revalId]
pointsAv=Table[{Log[xID],Log[-x]},{x,valAv[[2;;4]]}];


valRout={xRoutDIOapr,xRoutDIO,xRoutCUR,xRoutINV}//.Join[rval,rfunc,revalWLr,revalId]
pointsRout=Table[{Log[xID],Log[x]},{x,valRout[[2;;3]]}];


latexAv={
	xLaTeX["A_v=-0.9860(\\approx-1.0488)",0.7],
	xLaTeX["A_v=-16.4804",0.7],
	xLaTeX["A_v=-32.1939",0.7]};


latexRout={
	xLaTeX["R_{out}=\\SI{6.648}{k\\ohm}(\\approx\\SI{7.071}{k\\ohm})",0.7],
	xLaTeX["R_{out}=\\SI{111.1}{k\\ohm}",0.7]};


plsLoad={Directive[Lighter[Red],Thickness[0.002],Dashed],Red,Blue,Directive[Lighter[Gray],Dashed]};


legendLoad={
	xLaTeX["\\text{Diode}\\approx",0.7],xLaTeX["\\text{Diode}",0.7],
	xLaTeX["\\text{Current}",0.7],xLaTeX["\\text{Push-Pull}",0.7]};


ff1=-{xAvDIOapr,xAvDIO,xAvCUR,xAvINV}//.Join[rval,rfunc,revalWLr];
fig1=LogLogPlot[Evaluate@ff1,{ID,0.01m,1m},PlotRange->{(Sqrt[2]/2) 1,(Sqrt[2]) 100},
PlotStyle->plsLoad,
FrameTicks->{xTicksLog[{-5,-3,1}],xTicksLog[{0,2,1}]},
GridLines->{xGridsLog[{-5,-3,1},{10^-5,10^-3}],xGridsLog[{0,2,1},{1,100}]},
Epilog->Join[xPoint[pointsAv],{
	Text[latexAv[[1]],pointsAv[[1]]+{Log[1],Log[1.3]}],
	Text[latexAv[[2]],pointsAv[[2]]+{Log[0.8],Log[0.75]}],
	Text[latexAv[[3]],pointsAv[[3]]+{Log[1.7],Log[1.1]}]
}],
PlotLegends->xLegend[legendLoad,{"Column",2},{Right,Top}],
Evaluate@xFigLabel["I_D~(\\si{A})","|A_v|"],
Evaluate@xFigConfig]


ff2={xRoutDIOapr,xRoutDIO,xRoutCUR,xRoutINV}//.Join[rval,rfunc,revalWLr];
fig2=LogLogPlot[Evaluate@ff2,{ID,0.01m,1m},
PlotStyle->plsLoad,
FrameTicks->{xTicksLog[{-5,-3,1}],xTicksLog[{3,6,1}]},
GridLines->{xGridsLog[{-5,-3,1},{10^-5,10^-3}],xGridsLog[{3,6,1},{10^3,10^6}]},
Epilog->Join[xPoint[pointsRout],{
	Text[latexRout[[1]],pointsRout[[1]]+{Log[2],Log[1.3]}],
	Text[latexRout[[2]],pointsRout[[2]]+{Log[0.7],Log[0.75]}]
}],
PlotLegends->xLegend[legendLoad,{"Column",2},{Right,Top}],
Evaluate@xFigLabel["I_D~(\\si{A})","R_{out}~(\\si{\\ohm})"],
Evaluate@xFigConfig]


fig={fig1,fig2}


figWithName=Transpose[{fig,nbPDFNameList[[1;;Length[fig]]]}];
ReleaseHold[(Hold[Export[#[[2]],#[[1]],ImageResolution->2000]]&)/@figWithName]
