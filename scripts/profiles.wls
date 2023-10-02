(* Profile functions in the steady flow (fully resolved) *)

(* This script plots the profile functions for u,w and q obtained from the full system at two different x values *)

(* clear the scope *)
Quiet[Remove["Global`*"];Remove["FreeSurf`*"]];

(* set current directory and directory of the engine folder *)
currentDirectory=DirectoryName[$InputFileName];
engineDirectory=FileNameJoin[{ParentDirectory[currentDirectory],"engine"}];

(* load engine *)
Get[FileNameJoin[{engineDirectory,"FreeSurf.wl"}]];

(* set input parameters here *)
Fr=.6;n=90;eta=.09;sigma=.25;
xMin=-.2;xMax=.6;

(* create a legend *)
legend=Framed[LineLegend[ColorData[106,#]&/@{1,2,3},{"0.2","0.4","0.6"},LegendLabel->"H"],Background->White,FrameMargins->{{0,0},{0,0}} ]; 

(* calculate results *)
{freeSurf0p2,freeSurf0p4,freeSurf0p6}=FreeSurf[#,Fr,eta,sigma,n]&/@{.2,.4,.6};

(* extract data for u*)
{u0p2,u0p4,u0p6}={u[x,\[Zeta]] /.#,\[Zeta]}&/@{freeSurf0p2,freeSurf0p4,freeSurf0p6};

(* create plots for u *)
{pUxMin,pUxMax}=ParametricPlot[
Evaluate[{u0p2,u0p4,u0p6} /. x->#],{\[Zeta],0,1},
PlotRange->{{.8,1.6},{0,1}},
AspectRatio->.5,
GridLines->{{.9,{1,Directive[Black,Thick]},1.1,1.2,1.3,1.4,1.5},Automatic},
GridLinesStyle->{Black,Opacity[1]},
FrameStyle->Black,
FrameLabel->{"u[\[Zeta]]","\[Zeta]"},
PlotTheme->{"Frame","BoldColor","MediumLines"},
LabelStyle->{Directive[14],Black}
]&/@{xMin,xMax};

(* extract data for w*)
{w0p2,w0p4,w0p6}={w[x,\[Zeta]] /.#,\[Zeta]}&/@{freeSurf0p2,freeSurf0p4,freeSurf0p6};

(* create plots for u *)
{pWxMin,pWxMax}=ParametricPlot[
Evaluate[{w0p2,w0p4,w0p6} /. x->#],{\[Zeta],0,1},
PlotRange->{{-.5,.7},{0,1}},
AspectRatio->.5,
Axes->False,
FrameStyle->Black,
FrameLabel->{"w[\[Zeta]]","\[Zeta]"},
GridLines->{{-.4,-.3,-.2,-.1,{0,Directive[Black,Thick]},.1,.2,.3,.4,.5,.6},Automatic},
GridLinesStyle->{Opacity[1],Black},
PlotTheme->{"Frame","BoldColor","MediumLines"},
LabelStyle->{Directive[14],Black}
]& /@ {xMin,xMax};

(* extract data for q*)
{q0p2,q0p4,q0p6}={q[x,\[Zeta]] /.#,\[Zeta]}&/@{freeSurf0p2,freeSurf0p4,freeSurf0p6};

(* create plots for q *)
{pQxMin,pQxMax}=ParametricPlot[
Evaluate[{q0p2,q0p4,q0p6} /. x->#],{\[Zeta],0,1},
PlotRange->{{-.06,0.05},{0,1}},
AspectRatio->.5,
Axes->False,
GridLines->{{-.04,-.02,{0,Directive[Black,Thick]},.02,.04},Automatic},
GridLinesStyle->{Black,Opacity[1]},
FrameStyle->Black,
FrameLabel->{"q[\[Zeta]]","\[Zeta]"},
PlotTheme->{"Frame","BoldColor","MediumLines"},
LabelStyle->{Directive[14],Black}
]& /@ {xMin,xMax};

(* set output directory and save plots *)
   
outputDir=FileNameJoin[{ParentDirectory[currentDirectory],"plots","profile"}];

Quiet[CreateDirectory[outputDir]];

destinationFiles=FileNameJoin[{outputDir,#}]& /@{"legend.pdf","pUxMin.pdf","pUxMax.pdf","pWxMin.pdf","pWxMax.pdf","pQxMin.pdf","pQxMax.pdf"};

MapThread[Export,{destinationFiles, {legend,pUxMin,pUxMax,pWxMin,pWxMax,pQxMin,pQxMax}}];

Print["Plots saved to " <> outputDir <> ". Terminating..."]
