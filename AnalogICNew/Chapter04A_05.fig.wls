#!/usr/bin/env wolframscript
(* ::Package:: *)

ClearAll["Global`*"]


SetDirectory[NotebookDirectory[]];
Needs["MaTeX`"];


nbFileName="Chapter04A_05.fig.wls";
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


xCgsv=2.00f;
xCgdv=0.50f;
xCbdv=10.0f;
xCbsv=10.0f;
xCL=1p;


rval={kB->xkB,q->xq,m0->xm0,T0->xT0,
KxN->xKxN,KxP->xKxP,lambdaN->xlambdaN,lambdaP->xlambdaP,
Cgsv->xCgsv,Cgdv->xCgdv,Cbdv->xCbdv,Cbsv->xCbsv,CL->xCL};


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
xCgs=(Cgsv)&;
xCgd=(Cgdv)&;
xCbs=(Cbsv)&;
xCbd=(Cbdv)&;


rfunc={Kx->xKx,L->xL,W->xW,lambda->xlambda,beta->xbeta,gm->xgm,gds->xgds,
Cgs->xCgs,Cgd->xCgd,Cbs->xCbs,Cbd->xCbd
};


revalWLr={WLr->xWLr};
revalId={ID->xID};


xAvDIO=-gm[1](gm[2]+gds[1]+gds[2])^(-1);
xAvCUR=-gm[1](gds[1]+gds[2])^(-1);
xAvINV=-(gm[1]+gm[2])(gds[1]+gds[2])^(-1);


xRoutDIO=(gm[2]+gds[1]+gds[2])^(-1);
xRoutCUR=(gds[1]+gds[2])^(-1);
xRoutINV=(gds[1]+gds[2])^(-1);


xCmDIO=Cgd[1];
xCmCUR=Cgd[1];
xCmINV=Cgd[1]+Cgd[2];


xCoutDIO=Cbd[1]+Cbd[2]+Cgs[2]+CL;
xCoutCUR=Cbd[1]+Cbd[2]+Cgd[2]+CL;
xCoutINV=Cbd[1]+Cbd[2]+CL;


rDIO={Av->xAvDIO,Rout->xRoutDIO,Cm->xCmDIO,Cout->xCoutDIO};
rCUR={Av->xAvCUR,Rout->xRoutCUR,Cm->xCmCUR,Cout->xCoutCUR};
rINV={Av->xAvINV,Rout->xRoutINV,Cm->xCmINV,Cout->xCoutINV};


xp1=-(Cout+Cm)^(-1)Rout^(-1);
xz1=-(Cm)^(-1)Rout^(-1)Av;


xAvs=Av (1-s/z1)/(1-s/p1);
xs=2 Pi I fx;


plsLoad={Red,Blue,Directive[Lighter[Gray],Dashed]};
plsPZ={Black,Gray};


legendLoad={xLaTeX["\\text{Diode}",0.7],xLaTeX["\\text{Current}",0.7],xLaTeX["\\text{Push-Pull}",0.7]};
legendPZ={xLaTeX["-f_{p1}",0.7],xLaTeX["+f_{z1}",0.7]};


fPZ={xp1,xz1}/(2Pi)//.{rDIO,rCUR,rINV}//.Join[rval,rfunc,revalWLr,revalId]
pointsPZ=Table[{{Log[xID],Log[Abs[x[[1]]]]},{Log[xID],Log[Abs[x[[2]]]]}},{x,fPZ}];


latexP={
	xLaTeX["f_{p1}=-\\SI{23.41}{MHz}",0.7],
	xLaTeX["f_{p1}=-\\SI{1.40}{MHz}",0.7],
	xLaTeX["f_{p1}=-\\SI{1.40}{MHz}",0.7]};
latexZ={
	xLaTeX["f_{z1}=-\\SI{47.21}{GHz}",0.7],
	xLaTeX["f_{z1}=-\\SI{47.21}{GHz}",0.7],
	xLaTeX["f_{z1}=-\\SI{46.11}{GHz}",0.7]};


ff1=Abs[{xp1,xz1}/(2Pi)]//.{rDIO}//.Join[rval,rfunc,revalWLr];
fig1=LogLogPlot[Evaluate@ff1,{ID,0.01m,1m},PlotRange->{(Sqrt[2]/2)10^3,(Sqrt[2])10^12},
PlotStyle->plsPZ,
FrameTicks->{xTicksLog[{-5,-3,1}],xTicksLog[{3,12,3}]},
GridLines->{xGridsLog[{-5,-3,1},{10^-5,10^-3}],xGridsLog[{3,12,1},{10^3,10^12}]},
PlotLegends->xLegend[legendPZ,{"Column",1},{Right,Bottom}],
Epilog->Join[xPoint[pointsPZ[[1]]],{
	Text[latexP[[1]],pointsPZ[[1]][[1]]+{Log[1.8],Log[0.5]}],
	Text[latexZ[[1]],pointsPZ[[1]][[2]]+{Log[1.8],Log[0.5]}]
}],
Evaluate@xFigLabel["I_D~(\\si{A})","f_{p1},f_{z1}~(\\si{Hz})"],
Evaluate@xFigConfig]


ff2=Abs[{xp1,xz1}/(2Pi)]//.{rCUR}//.Join[rval,rfunc,revalWLr];
fig2=LogLogPlot[Evaluate@ff2,{ID,0.01m,1m},PlotRange->{(Sqrt[2]/2)10^3,(Sqrt[2])10^12},
PlotStyle->plsPZ,
FrameTicks->{xTicksLog[{-5,-3,1}],xTicksLog[{3,12,3}]},
GridLines->{xGridsLog[{-5,-3,1},{10^-5,10^-3}],xGridsLog[{3,12,1},{10^3,10^12}]},
PlotLegends->xLegend[legendPZ,{"Column",1},{Right,Bottom}],
Epilog->Join[xPoint[pointsPZ[[2]]],{
	Text[latexP[[2]],pointsPZ[[2]][[1]]+{Log[1.8],Log[0.5]}],
	Text[latexZ[[2]],pointsPZ[[2]][[2]]+{Log[1.8],Log[0.5]}]
}],
Evaluate@xFigLabel["I_D~(\\si{A})","f_{p1},f_{z1}~(\\si{Hz})"],
Evaluate@xFigConfig]


ff3=Abs[{xp1,xz1}/(2Pi)]//.{rINV}//.Join[rval,rfunc,revalWLr];
fig3=LogLogPlot[Evaluate@ff3,{ID,0.01m,1m},PlotRange->{(Sqrt[2]/2)10^3,(Sqrt[2])10^12},
PlotStyle->plsPZ,
FrameTicks->{xTicksLog[{-5,-3,1}],xTicksLog[{3,12,3}]},
GridLines->{xGridsLog[{-5,-3,1},{10^-5,10^-3}],xGridsLog[{3,12,1},{10^3,10^12}]},
PlotLegends->xLegend[legendPZ,{"Column",1},{Right,Bottom}],
Epilog->Join[xPoint[pointsPZ[[3]]],{
	Text[latexP[[3]],pointsPZ[[3]][[1]]+{Log[1.8],Log[0.5]}],
	Text[latexZ[[3]],pointsPZ[[3]][[2]]+{Log[1.8],Log[0.5]}]
}],
Evaluate@xFigLabel["I_D~(\\si{A})","f_{p1},f_{z1}~(\\si{Hz})"],
Evaluate@xFigConfig]


ff4=Abs[xAvs//.{s->xs,z1->xz1,p1->xp1}]//.{rDIO,rCUR,rINV}//.Join[rval,rfunc,revalWLr,revalId];
fig4=LogLogPlot[Evaluate@ff4,{fx,10^3,10^15},PlotRange->{(Sqrt[2]/2)10^-4,(Sqrt[2])10^2},
PlotStyle->plsLoad,
FrameTicks->{xTicksLog[{3,15,3}],xTicksLog[{-4,2,1}]},
GridLines->{xGridsLog[{3,15,1},{10^3,10^15}],xGridsLog[{-4,2,1},{10^-4,10^2}]},
PlotLegends->xLegend[legendLoad,{"Column",1},{Right,Top}],
Evaluate@xFigLabel["f~(\\si{Hz})","|A_v(s)|"],
Evaluate@xFigConfig]


ff5=Arg[(xAvs)//.{s->xs,z1->xz1,p1->xp1}]/Pi//.{rDIO,rCUR,rINV}//.Join[rval,rfunc,revalWLr,revalId];
fig5=LogLinearPlot[Evaluate@ff5,{fx,10^3,10^15},PlotRange->{-1.1,1.1},
PlotStyle->plsLoad,
FrameTicks->{xTicksLog[{3,15,3}],xTicks[{-1,1,0.25},{3,2}]},
GridLines->{xGridsLog[{3,15,1},{10^3,10^15}],xGrids[{-1,1,0.25}]},
PlotLegends->xLegend[legendLoad,{"Column",1},{Right,Bottom}],
Evaluate@xFigLabel["f~(\\si{Hz})","\\arg A_v(s)/\\pi"],
Evaluate@xFigConfig]


fig={fig1,fig2,fig3,fig4,fig5}


figWithName=Transpose[{fig,nbPDFNameList[[1;;Length[fig]]]}];
ReleaseHold[(Hold[Export[#[[2]],#[[1]],ImageResolution->2000]]&)/@figWithName]
