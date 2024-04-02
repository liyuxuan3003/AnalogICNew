#!/usr/bin/env wolframscript
(* ::Package:: *)

ClearAll["Global`*"]


SetDirectory[NotebookDirectory[]];
Needs["MaTeX`"];


nbFileName="Chapter02B_01.fig.wls";
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
xTicksLogX=(Table[{x,xLaTeX["10^{"<>ToString[x]<>"}",0.5]},{x,#1[[1]],#1[[2]],#1[[3]]}])&;


xGrids=(Table[x,{x,#1[[1]],#1[[2]],#1[[3]]}])&;
xGridsLog=(Flatten[Table[Table[{If[#2[[1]]<=(n 10^x)<=#2[[2]],(n 10^x),{}]},{n,1,9,1}],{x,#1[[1]],#1[[2]],#1[[3]]}]])&;


xText=Style[#1,FontSize->#2,FontFamily->"Noto Serif CJK SC"]&;
xTextConsolas=Style[#1,FontSize->#2,FontFamily->"Consolas"]&;


f=10^(-15);
p=10^(-12);
n=10^(-9);
u=10^(-6);
m=10^(-3);
k=10^(3);
M=10^(6);


(*-*)


xkB=1.38 10^(-23);
xni=1 45 10^(10);
xe0=8.85 10^(-14);
xesi=11.7 e0;
xT=300;
xq=1.602 10^(-19);
xVt=kB T/q;
lval={kB->xkB,ni->xni,e0->xe0,esi->xesi,T->xT,q->xq,Vt->xVt};


xphi0=Vt Log[(NA ND)/ni^2];


fig1=LogLinearPlot[Evaluate@Table[xphi0//.{vD->0}//.lval,{ND,{10^16,10^17}}],{NA,10^15,10^20},
PlotRange->{0,1},
Evaluate@xFigLabel["N_A~(\\si{cm^{-3}})","\\phi_0~(\\si{V})"],
FrameTicks->{xTicksLog[{15,20,1}],xTicks[{0,1,0.2},{2,1}]},
GridLines->{xGridsLog[{15,20,1},{10^15,10^20}],xGrids[{0,1,0.1}]},
PlotStyle->{Black,Gray},
PlotLegends->Placed[LineLegend[{
	xLaTeX["\\phi_0, N_D=\\SI{1e16}{cm^{-3}}",0.8],
	xLaTeX["\\phi_0, N_D=\\SI{1e17}{cm^{-3}}",0.8]},
	LegendFunction->"Panel",LegendLayout -> {"Column",1},
	LegendMarkerSize->{22,1}],{Right,Bottom}],
Evaluate@xFigConfig]


fig2=Plot3D[Evaluate[xphi0/.{NA->10^NAexp,ND->10^NDexp}//.lval],{NAexp,14.8,20.2},{NDexp,14.8,20.2},
Ticks->{xTicksLogX[{15,20,1}],xTicksLogX[{15,20,1}],xTicks[{0,1,0.2},{2,1}]},
Mesh->{Table[x,{x,15,20,1}],Table[x,{x,15,20,1}]},
PlotRange->{-0.1,1.1},
ColorFunction->(ColorData["DarkRainbow"][#3]&),
ColorFunctionScaling->{False},
Filling->0,
FillingStyle->{Directive[Black,Opacity[0.1]]},
Evaluate@xFig3DLabel["N_A~(\\si{cm^{-3}})","N_D~(\\si{cm^{-3}})","\\phi_0~(\\si{V})"],
Evaluate@xFig3DConfig]


(*-*)


xn=Sqrt[(2esi(phi0-vD)NA)/(q ND(NA+ND))];
xp=Sqrt[(2esi(phi0-vD)ND)/(q NA(NA+ND))];


fig3=LogLogPlot[Evaluate@Table[10^(-2){xn+xp,xp,xn}//.{vD->0,phi0->xphi0}//.lval,{ND,{10^16}}],{NA,10^15,10^20},
PlotRange->{1/1.5 10^-9,1.5 10^-6},
Evaluate@xFigLabel["N_A~(\\si{cm^{-3}})","x~(\\si{m})"],
FrameTicks->{xTicksLog[{15,20,1}],xTicksLog[{-9,-6,1}]},
GridLines->{xGridsLog[{15,20,1},{10^15,10^20}],xGridsLog[{-9,-6,1},{10^-9,10^-6}]},
PlotStyle->{Black,Red,Blue,Directive[Black,Dashed],Directive[Red,Dashed],Directive[Blue,Dashed]},
PlotLegends->Placed[LineLegend[{
	xLaTeX["x_d",0.8],
	xLaTeX["x_p",0.8],
	xLaTeX["x_n",0.8],
	xLaTeX["x_d",0.8],
	xLaTeX["x_p",0.8],
	xLaTeX["x_n",0.8]},
	LegendFunction->"Panel",LegendLayout -> {"Column",1},
	LegendMarkerSize->{22,0.5},LegendLabel->xLaTeX["N_D=10^{16}\\si{cm^{-3}}",0.7]],{Right,Bottom}],
Evaluate@xFigConfig]


fig5=LogPlot[10^{-20}+Re@Evaluate[10^(-2){xn+xp}//.{phi0->xphi0,NA->10^16,ND->10^18}//.lval],{vD,-3,1},
Evaluate@xFigLabel["v_D~(V)","x~(\\si{m})"],
PlotRange->{1/1.5 10^-9,1.5 10^-6},
FrameTicks->{xTicks[{-3,1,0.5},{2,1}],xTicksLog[{-9,-6,1}]},
GridLines->{xGrids[{-3,1,0.5}],xGridsLog[{-9,-6,1},{10^-9,10^-6}]},
PlotStyle->{Black},
PlotLegends->Placed[LineLegend[{
	xLaTeX["x_d",0.8]},
	LegendFunction->"Panel",LegendLayout -> {"Column",2},
	LegendMarkerSize->{22,1},
	LegendLabel->xLaTeX["N_D=10^{16}\\si{cm^{-3}}~~N_A=10^{18}\\si{cm^{-3}}",0.7]],{Left,Bottom}],
Evaluate@xFigConfig]


fxn=Evaluate[Log10[10^(-2)xn]//.{vD->0,phi0->xphi0}/.{NA->10^NAexp,ND->10^NDexp}//.lval];
fig6=Plot3D[fxn,{NAexp,14.8,20.2},{NDexp,14.8,20.2},
Mesh->{Table[x,{x,15,20,1}],Table[x,{x,15,20,1}]},
Ticks->{xTicksLogX[{15,20,1}],xTicksLogX[{15,20,1}],xTicksLogX[{-10,-5,1}]},
PlotRange->{-10,-5},
Filling->Bottom,
Evaluate@xFig3DLabel["N_A~(\\si{cm^{-3}})","N_D~(\\si{cm^{-3}})","x_n~(\\si{m})"],
Evaluate@xFig3DConfig]


fxp=Evaluate[Log10[10^(-2)xp]//.{vD->0,phi0->xphi0}/.{NA->10^NAexp,ND->10^NDexp}//.lval];
fig7=Plot3D[fxp,{NAexp,14.8,20.2},{NDexp,14.8,20.2},
Mesh->{Table[x,{x,15,20,1}],Table[x,{x,15,20,1}]},
Ticks->{xTicksLogX[{15,20,1}],xTicksLogX[{15,20,1}],xTicksLogX[{-10,-5,1}]},
PlotRange->{-10,-5},
Filling->Bottom,
Evaluate@xFig3DLabel["N_A~(\\si{cm^{-3}})","N_D~(\\si{cm^{-3}})","x_p~(\\si{m})"],
Evaluate@xFig3DConfig]


fxd=Evaluate[Log10[10^(-2)(xp+xn)]//.{vD->0,phi0->xphi0}/.{NA->10^NAexp,ND->10^NDexp}//.lval];
fig8=Plot3D[fxd,{NAexp,14.8,20.2},{NDexp,14.8,20.2},
Mesh->{Table[x,{x,15,20,1}],Table[x,{x,15,20,1}]},
Ticks->{xTicksLogX[{15,20,1}],xTicksLogX[{15,20,1}],xTicksLogX[{-10,-5,1}]},
PlotRange->{-10,-5},
Filling->Bottom,
Evaluate@xFig3DLabel["N_A~(\\si{cm^{-3}})","N_D~(\\si{cm^{-3}})","x_d~(\\si{m})"],
Evaluate@xFig3DConfig]


(*-*)


xQj=Sqrt[(2esi q NA ND)/(NA+ND)]Sqrt[phi0-vD];
xE0=Sqrt[(2q NA ND)/(esi(NA+ND))] Sqrt[phi0-vD];
xCj=Sqrt[(esi q NA ND)/(2(NA+ND))]/Sqrt[phi0-vD];
xCj0=xCj//.{vD->0};
xCjm=xCj0/(1-vD/phi0)^(mx);


fig9=Plot3D[Log10[xQj]//.{vD->0,phi0->xphi0}/.{NA->10^NAexp,ND->10^NDexp}//.lval,{NAexp,14.8,20.2},{NDexp,14.8,20.2},
Mesh->{Table[x,{x,15,20,1}],Table[x,{x,15,20,1}]},
Ticks->{xTicksLogX[{15,20,1}],xTicksLogX[{15,20,1}],xTicksLogX[{-8,-6,1}]},
Filling->Bottom,
Evaluate@xFig3DLabel["N_A~(\\si{cm^{-3}})","N_D~(\\si{cm^{-3}})","Q_j~(\\si{C\\cdot cm^{-2}})"],
Evaluate@xFig3DConfig]


fig10=Plot3D[Log10[xE0]//.{vD->0,phi0->xphi0}/.{NA->10^NAexp,ND->10^NDexp}//.lval,{NAexp,14.8,20.2},{NDexp,14.8,20.2},
Mesh->{Table[x,{x,15,20,1}],Table[x,{x,15,20,1}]},
Ticks->{xTicksLogX[{15,20,1}],xTicksLogX[{15,20,1}],xTicksLogX[{4,6,1}]},Exclusions->{(Log10[xE0]//.{vD->0,phi0->xphi0}/.{NA->10^NAexp,ND->10^NDexp}//.lval)==Log10[3 10^5]},
Filling->Bottom,
Evaluate@xFig3DLabel["N_A~(\\si{cm^{-3}})","N_D~(\\si{cm^{-3}})","E_0~(\\si{V\\cdot cm^{-1}})"],
Evaluate@xFig3DConfig]


fig11=Plot3D[Log10[xCj]//.{vD->0,phi0->xphi0}/.{NA->10^NAexp,ND->10^NDexp}//.lval,{NAexp,14.8,20.2},{NDexp,14.8,20.2},
Mesh->{Table[x,{x,15,20,1}],Table[x,{x,15,20,1}]},
Ticks->{xTicksLogX[{15,20,1}],xTicksLogX[{15,20,1}],xTicksLogX[{-8,-6,1}]},
Filling->Bottom,
Evaluate@xFig3DLabel["N_A~(\\si{cm^{-3}})","N_D~(\\si{cm^{-3}})","C_j~(\\si{F\\cdot cm^{-2}})"],
Evaluate@xFig3DConfig]


fig12=LogPlot[Evaluate[{Piecewise[{{xQj,vD<phi0},{10^-10,vD>=phi0}}]}//.{phi0->xphi0,NA->10^16,ND->10^18}//.lval],{vD,-3,1},
Evaluate@xFigLabel["v_D~(V)","Q_j~(\\si{C\\cdot cm^{-2}})"],PlotRange->{(Sqrt[2]/2) 10^-9,(Sqrt[2])10^-6},
PlotPoints->500,
FrameTicks->{xTicks[{-3,1,0.5},{2,1}],xTicksLog[{-9,-6,1}]},
GridLines->{xGrids[{-3,1,0.5}],xGridsLog[{-9,-6,1},{10^-9,10^-6}]},
PlotStyle->{Black},
PlotLegends->Placed[LineLegend[{
	xLaTeX["Q_j",0.8]},
	LegendFunction->"Panel",LegendLayout -> {"Column",2},
	LegendMarkerSize->{22,1},
	LegendLabel->xLaTeX["N_D=10^{16}\\si{cm^{-3}}~~N_A=10^{18}\\si{cm^{-3}}",0.7]],{Left,Bottom}],
Epilog->{
	Gray,Dashed,Thickness[0.002],
	Line[{{xphi0,Log[10^-5]},{xphi0,Log[10^-10]}}]//.{phi0->xphi0,NA->10^16,ND->10^18}//.lval,
	Text[xLaTeX["v_D=\\phi_0",0.8],{0.3,Log[5 10^-7]}]
},
Evaluate@xFigConfig]


(*-*)


xEcirc=3 10^5
xEval0=xE0//.{phi0->xphi0,NA->10^16,ND->10^18}//.lval//.{vD->0}
xvDcirc=vD//.Solve[xE0==xEcirc,{vD}][[1]]//.{phi0->xphi0,NA->10^16,ND->10^18}//.lval


fig13=LogPlot[Evaluate[Piecewise[{{xE0,vD<phi0},{10^7,vD>phi0}}]//.{phi0->xphi0,NA->10^16,ND->10^18}//.lval],{vD,-40,5},
Evaluate@xFigLabel["v_D~(V)","-E_0~(\\si{V\\cdot cm^{-1}})"],PlotRange->{Sqrt[2]10^4/2,Sqrt[2]10^6},
FrameTicks->{xTicks[{-40,5,5},{2,1}],xTicksLog[{4,6,1}]},
GridLines->{xGrids[{-40,5,5}],xGridsLog[{4,6,1},{10^4,10^6}]},
PlotStyle->{Black},
PlotLegends->Placed[LineLegend[{
	xLaTeX["E_0",0.8]},
	LegendFunction->"Panel",LegendLayout -> {"Column",2},
	LegendMarkerSize->{22,1},
	LegendLabel->xLaTeX["N_D=10^{16}\\si{cm^{-3}}~~N_A=10^{18}\\si{cm^{-3}}",0.7]],{Left,Bottom}],
Prolog->{
	Gray,Dashed,Thickness[0.002],
	Line[{{xphi0,Log[10^3]},{xphi0,Log[10^7]}}]//.{phi0->xphi0,NA->10^16,ND->10^18}//.lval,
	Line[{{xvDcirc,Log[10^3]},{xvDcirc,Log[10^7]}}],
	Line[{{-65,Log[xEcirc]},{12,Log[xEcirc]}}],
	Text[xLaTeX["v_{D}=\\phi_0",0.8],{-3,Log[7 10^5]}],
	Text[xLaTeX["v_{D}=-V_{BR}",0.8],{-33,Log[7 10^5]}],
	Text[xLaTeX["E_\\text{max}=3\\times 10^5\\si{V.cm^{-1}}",0.8],{-21,Log[xEcirc*1.2]}]
},
Epilog->{
	PointSize[0.015],Black,
	Point[{xvDcirc,Log[xEcirc]}],
	PointSize[0.010],White,
	Point[{xvDcirc,Log[xEcirc]}]
},
Evaluate@xFigConfig]


fig14=LogPlot[Evaluate[Table[Piecewise[{{xCjm,vD<phi0},{10^-10,vD>phi0}}],{mx,{1/2,1/3}}]//.{phi0->xphi0,NA->10^16,ND->10^18}//.lval],{vD,-3,1},
Evaluate@xFigLabel["v_D~(V)","C_j~(\\si{F\\cdot cm^{-2}})"],PlotRange->{(Sqrt[2]/2)10^-8,Sqrt[2]10^-6},
FrameTicks->{xTicks[{-3,1,0.5},{2,1}],xTicksLog[{-8,-6,1}]},
GridLines->{xGrids[{-3,1,0.5}],xGridsLog[{-8,-6,1},{10^-8,10^-6}]},
PlotStyle->{Black,Gray},
PlotLegends->Placed[LineLegend[{
	xLaTeX["Q_j,m=1/2",0.8],xLaTeX["Q_j,m=1/3",0.8]},
	LegendFunction->"Panel",LegendLayout -> {"Column",1},
	LegendMarkerSize->{22,1},
	LegendLabel->xLaTeX["N_D=10^{16}\\si{cm^{-3}}~~N_A=10^{18}\\si{cm^{-3}}",0.7]],{Left,Top}],
Epilog->{
	Gray,Dashed,Thickness[0.002],
	Line[{{xphi0,Log[10^-4]},{xphi0,Log[10^-11]}}]//.{phi0->xphi0,NA->10^16,ND->10^18}//.lval,
	Text[xLaTeX["v_D=\\phi_0",0.8],{0.3,Log[2 10^-8]}]
},
Evaluate@xFigConfig]


xVBR=(esi (NA+ND))/(2q NA ND) Ecirc^2;


fig15=Plot3D[-xVBR//.{Ecirc->xEcirc,phi0->xphi0}//.{NA->10^NAexp,ND->10^NDexp}//.lval,{NAexp,14.8,20.2},{NDexp,14.8,20.2},
Mesh->{Table[x,{x,15,20,1}],Table[x,{x,15,20,1}]},PlotRange->{-80,0},
Ticks->{xTicksLogX[{15,20,1}],xTicksLogX[{15,20,1}],xTicks[{-80,0,20},{3,0}]},
Filling->Bottom,
Evaluate@xFig3DLabel["N_A~(\\si{cm^{-3}})","N_D~(\\si{cm^{-3}})","-V_{BR}~(\\si{V})"],
Evaluate@xFig3DConfig]


fig16=LogLinearPlot[Evaluate@Table[-xVBR//.{Ecirc->xEcirc,phi0->xphi0}//.lval,{ND,{10^16,10^17}}],{NA,10^15,10^20},
Evaluate@xFigLabel["N_A~(\\si{cm^{-3}})","-V_{BR}~(\\si{V})"],PlotRange->{-85,5},
FrameTicks->{xTicksLog[{15,20,1}],xTicks[{-80,0,10},{3,0}]},
GridLines->{xGridsLog[{15,20,1},{10^15,10^20}],Automatic},
PlotStyle->{Black,Gray},
PlotLegends->Placed[LineLegend[{
	xLaTeX["-V_{BR},~N_D=\\SI{1e16}{cm^{-3}}",0.8],
	xLaTeX["-V_{BR},~N_D=\\SI{1e17}{cm^{-3}}",0.8]},
	LegendFunction->"Panel",LegendLayout -> {"Column",1},
	LegendMarkerSize->{22,1}],{Right,Bottom}],
Evaluate@xFigConfig]


xiD=Is (Exp[vD/Vt]-1);
xM=Piecewise[{{(1/(1+(vD/VBR)^nx)),vD<=0},{1,vD>0}}];
xIs=10^-15;


fig17=Plot[Evaluate[(xiD xM)//.{VBR->xVBR,Is->xIs,ND->10^16,NA->10^18,nx->3,Ecirc->xEcirc}//.lval],{vD,-30,0},
PlotRange->{-8.5f,0.5f},PlotStyle->Black,
FrameTicks->{xTicks[{-30,0,5},{2,0}],xTicksProp[{-8,0,1},{1,0},f]},
GridLines->{xGrids[{-30,0,1}],Automatic},
Evaluate@xFigLabel["v_{D}~(\\si{V})","i_D~(\\si{fA})"],
Prolog->{
	Gray,Dashed,Thickness[0.002],
	Line[{{xvDcirc,1},{xvDcirc,-10}}],
	Text[xLaTeX["v_{D}=-V_{BR}",0.8],{-25.5,-0.5f}]
},
Evaluate@xFigConfig]


fig18=Plot[Evaluate[(xiD xM)//.{VBR->xVBR,Is->xIs,ND->10^16,NA->10^18,nx->3,Ecirc->xEcirc}//.lval],{vD,-1,1},
PlotRange->{-0.05m,1.05m},PlotStyle->Black,
FrameTicks->{xTicks[{-1,1,0.2},{2,1}],xTicksProp[{0,1,0.1},{2,1},m]},
GridLines->{xGrids[{-1,1,0.1}],Automatic},
Prolog->{
	Gray,Dashed,Thickness[0.002],
	Line[{{xphi0,-1m},{xphi0,2m}}]//.{phi0->xphi0,NA->10^16,ND->10^18}//.lval,
	Text[xLaTeX["v_{D}=\\phi_0",0.8],{0.45,0.85m}]
},
Evaluate@xFigLabel["v_{D}~(\\si{V})","i_D~(\\si{mA})"],
Evaluate@xFigConfig]


fig={fig1,fig2,fig3,fig5,fig6,fig7,fig8,fig9,fig10,fig11,fig12,fig13,fig14,fig15,fig16,fig17,fig18}


figWithName=Transpose[{fig,nbPDFNameList[[1;;Length[fig]]]}];
ReleaseHold[(Hold[Export[#[[2]],#[[1]],ImageResolution->2000]]&)/@figWithName]


(* *)
