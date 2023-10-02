(* 2D Steady flow over bottom topography (fully resolved) *)

(* This script uses the fully resolved system to simulate the steady flow *)
(* Then the fluid height is plotted over the x axis.*)

(* clear the scope *)
Quiet[Remove["Global`*"];Remove["FreeSurf`*"]];

(* set current directory and directory of the engine folder *)
currentDirectory=DirectoryName[$InputFileName];
engineDirectory=FileNameJoin[{ParentDirectory[currentDirectory],"engine"}];

(* load engine *)
Get[FileNameJoin[{engineDirectory,"FreeSurf.wl"}]];

(* set parameters *)
H=.4;Fr=.6;n=90;\[Eta]=.09;\[Sigma]=.25;hb[x_]=\[Eta]/H E^(-x^2/(2\[Sigma]^2));dxhb[x_]=-((E^(-(x^2/(2 \[Sigma]^2))) x \[Eta])/(\[Sigma]^2 H));

(*initialize plotting parameters*)
plotHeight=170;plotWidth=500;labelHeight=17;
defaultPlotOptions=Options[#]&/@{Plot,DensityPlot,ParametricPlot,StreamPlot}; (*save default for restoring later*)
SetOptions[#,
ImageSize->{plotWidth,Automatic},
ImagePadding->{{30,13},{25,10}},
AspectRatio->H/2,
PlotRangePadding->None,
Frame->False,
GridLines->None,
Axes->False,
LabelStyle->Directive[14],
PlotRangeClipping->True
]& /@ {Plot,DensityPlot,ParametricPlot,StreamPlot};
SetOptions[#,PlotRange->{0,1.2}]&/@{Plot,DensityPlot,ParametricPlot};

(* load data *) 
freeSurf=FreeSurf[H,Fr,\[Eta],\[Sigma],n];

(* create plot *)
plotH=Plot[{(h[x] /.freeSurf)+hb[x],hb[x]},{x,-1,1},Filling->{1->{2}},PlotRange->{0,1.2},PlotLabel->(tag /.freeSurf)<>" Plotted (Free surface): h[x], hb[x]",ImageSize->{Automatic,plotHeight+labelHeight},FrameLabel->{"x","z"},Frame->True];

(* set output directory and save plots *)
identifier = "H" <> StringReplace[ToString[H], "." -> "p"] <> "Fr" <> 
   StringReplace[ToString[Fr], "." -> "p"];
   
outputDir=FileNameJoin[{ParentDirectory[currentDirectory],"plots","height",identifier}];

Quiet[CreateDirectory[outputDir]];

destinationFile=FileNameJoin[{outputDir,"h.pdf"}];
Export[destinationFile, plotH];

Print["Plot saved to " <> outputDir <> ". Terminating..."]

