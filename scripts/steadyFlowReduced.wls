#!/usr/bin/env wolframscript
(* ::Package:: *)

(* 2D Steady flow over bottom topography (depth projected) *)

(* This script repeats a single experiment and uses a different reduced system for the solution in each repetition. *)
(* The parameters shallowness and Froude number are fixed. The plots allow for a direct comparison of how flow variables are modeled by the different systems. *)
(* The variables we study are the height, the means of velocities and the nonhydrostatic pressure and the respective moments up to third order. *)


(********************* Preliminaries **************************************************************)

(* clear the scope *)
contexts = (# <> "*") & /@ {"Global`", "FreeSurf`", "DSM0`", "DSM1`", 
    "DSM2`", "DSM3`", "GN`"};
Quiet[Remove /@ contexts];

(* set current directory and directory of the engine folder *)
currentDirectory=DirectoryName[$InputFileName];
engineDirectory=FileNameJoin[{ParentDirectory[currentDirectory],"engine"}];

(* load engines *)
Get[FileNameJoin[{engineDirectory, "FreeSurf.wl"}]] 
Get[FileNameJoin[{engineDirectory, "DSM0.wl"}]] 
Get[FileNameJoin[{engineDirectory, "DSM1.wl"}]] 
Get[FileNameJoin[{engineDirectory, "DSM2.wl"}]] 
Get[FileNameJoin[{engineDirectory, "DSM3.wl"}]] 
Get[FileNameJoin[{engineDirectory, "GN.wl"}]] 


(********************* Set parameters *************************************************************)

(* setup parameters *)
Fr = .6; H = .4; n = 90; eta=.09; sigma=.25;

(********************* Calculate results **********************************************************)

{dsm0, dsm1, dsm2, dsm3, full, 
   gn} = (#[H, Fr, eta, sigma, n]) & /@ {DSM0, DSM1, DSM2, DSM3, FreeSurf, GN};

(********************* Set up plotting environment ************************************************)

plotHeight = 170; labelHeight = 17;
Clear[pRange];
defaultPlotOptions = 
  Options[#] & /@ {Plot, 
    ListPlot}; (*save default for restoring later*)
SetOptions[#,
     ImageSize -> {Automatic, plotHeight},
     ImagePadding -> {{70, 50}, {40, 25}},
     AspectRatio -> 1/2,
     PlotRangePadding -> None,
     Frame -> None,
     Axes -> None,
     LabelStyle -> Directive[Black, 14],
     PlotRangeClipping -> True,
     GridLines -> None,
     PlotRange -> pRange
     ]; & /@ {Plot, ListPlot};


SetOptions[ListPlot, PlotStyle -> Black, PlotMarkers -> {"\[SmallCircle]", 28}];
SetOptions[Plot, PlotTheme -> {"BoldColor"}];

(* create legend *)
legend = Framed[
  LineLegend[
   Join[{Black}, ColorData[106, #] & /@ {1, 2, 3, 4, 5}], {"Full", 
    "DSM-0", "GN", "DSM-1", "DSM-2", "DSM-3"}], Background -> White, 
  FrameMargins -> {{0, 0}, {0, 0}} ];

niceBlue = ColorData[97, 1]; fillingBlue = 
 Directive[Opacity[0.2], niceBlue]; niceGreen = ColorData[97, 3];
 
(********************* Plot height and means ******************************************************)

pRange = Automatic; (* you can play with the plot range here or use wolfram's automatic range determination *)

plotHFull = 
  ListPlot[Table[{x, h[x] /. full}, {x, -1, 1, .05}], 
   PlotLegends -> Placed[tag /. full, Right],
   FrameLabel -> {x, h}, Frame -> True, GridLines -> Automatic];
plotHMoments = 
  Plot[Evaluate[h[x] /. # & /@ {dsm0, gn, dsm1, dsm2, dsm3}], {x, -1, 1}, 
   Filling -> Axis, FillingStyle -> fillingBlue, 
   PlotLegends -> Placed[tag /. # & /@ {dsm0, gn, dsm1, dsm2, dsm3}, Right]];
plotH = Show[plotHFull, plotHMoments];

pRange = Automatic;

plotUmFull = 
  ListPlot[Table[{x, um[x] /. full}, {x, -1, 1, .05}], 
   GridLines -> Automatic, Frame -> True, 
   PlotLegends -> Placed[tag /. full, Right], 
   FrameLabel -> {x, Subscript[u, m]}];
plotUmMoments = 
  Plot[Evaluate[um[x] /. # & /@ {dsm0, gn, dsm1, dsm2, dsm3}], {x, -1, 1}, 
   PlotLegends -> Placed[tag /. # & /@ {dsm0, gn, dsm1, dsm2, dsm3}, Right]];
plotUm = Show[plotUmFull, plotUmMoments];
pRange = Automatic;

plotWmFull = 
  ListPlot[Table[{x, wm[x] /. full}, {x, -1, 1, .05}], 
   GridLines -> Automatic, Frame -> True, 
   PlotLegends -> Placed[tag /. full, Right], 
   FrameLabel -> {x, Subscript[w, m]}];
plotWmMoments = 
  Plot[Evaluate[wm[x] /. # & /@ {dsm0, gn, dsm1, dsm2, dsm3}], {x, -1,
     1}, PlotLegends -> 
    Placed[tag /. # & /@ {dsm0, gn, dsm1, dsm2, dsm3}, Right]];
plotWm = Show[plotWmFull, plotWmMoments];

pRange = Automatic;

plotQmFull = 
  ListPlot[Table[{x, qm[x] /. full}, {x, -1, 1, .05}], 
   GridLines -> Automatic, Frame -> True, 
   PlotLegends -> Placed[tag /. full, Right], 
   FrameLabel -> {x, Subscript[q, m]}];
plotQmMoments = 
  Plot[Evaluate[qm[x] /. # & /@ {dsm0, gn, dsm1, dsm2, dsm3}], {x, -1,
     1}, PlotLegends -> 
    Placed[tag /. # & /@ {dsm0, gn, dsm1, dsm2, dsm3}, Right]];
plotQm = Show[plotQmFull, plotQmMoments];


(********************* Plot first order moments ***************************************************)
pRange = Automatic;

plotAlpha1Full = 
  ListPlot[Table[{x, alpha1[x] /. full}, {x, -1, 1, .05}], 
   GridLines -> Automatic, Frame -> True, 
   PlotLegends -> Placed[tag /. full, Right], 
   FrameLabel -> {x, Subscript[\[Alpha], 1]}];
plotAlpha1Moments = 
  Plot[Evaluate[
    alpha1[x] /. # & /@ {dsm0, gn, dsm1, dsm2, dsm3}], {x, -1, 1}, 
   PlotLegends -> 
    Placed[tag /. # & /@ {dsm0, gn, dsm1, dsm2, dsm3}, Right]];
plotAlpha1 = Show[plotAlpha1Full, plotAlpha1Moments];

pRange = Automatic;

plotGamma1Full = 
  ListPlot[Table[{x, gamma1[x] /. full}, {x, -1, 1, .05}], 
   GridLines -> Automatic, Frame -> True, 
   PlotLegends -> Placed[tag /. full, Right], 
   FrameLabel -> {x, Subscript[\[Gamma], 1]}];
plotGamma1Moments = 
  Plot[Evaluate[
    gamma1[x] /. # & /@ {dsm0, gn, dsm1, dsm2, dsm3}], {x, -1, 1}, 
   PlotLegends -> 
    Placed[tag /. # & /@ {dsm0, gn, dsm1, dsm2, dsm3}, Right]];
plotGamma1 = Show[plotGamma1Full, plotGamma1Moments];

pRange = Automatic;

plotKappa1Full = 
  ListPlot[Table[{x, kappa1[x] /. full}, {x, -1, 1, .05}], 
   GridLines -> Automatic, Frame -> True, 
   PlotLegends -> Placed[tag /. full, Right], 
   FrameLabel -> {x, Subscript[\[Kappa], 1]}];
plotKappa1Moments = 
  Plot[Evaluate[
    kappa1[x] /. # & /@ {dsm0, gn, dsm1, dsm2, dsm3}], {x, -1, 1}, 
   PlotLegends -> 
    Placed[tag /. # & /@ {dsm0, gn, dsm1, dsm2, dsm3}, Right]];
plotKappa1 = Show[plotKappa1Full, plotKappa1Moments];

(********************* Plot second order moments **************************************************)

pRange = Automatic;

plotAlpha2Full = 
  ListPlot[Table[{x, alpha2[x] /. full}, {x, -1, 1, .05}], 
   GridLines -> Automatic, Frame -> True, 
   PlotLegends -> Placed[tag /. full, Right], 
   FrameLabel -> {x, Subscript[\[Gamma], 2]}];
plotAlpha2Moments = 
  Plot[Evaluate[
    alpha2[x] /. # & /@ {dsm0, gn, dsm1, dsm2, dsm3}], {x, -1, 1}, 
   PlotLegends -> 
    Placed[tag /. # & /@ {dsm0, gn, dsm1, dsm2, dsm3}, Right]];
plotAlpha2 = Show[plotAlpha2Full, plotAlpha2Moments];

pRange = Automatic;

plotGamma2Full = 
  ListPlot[Table[{x, gamma2[x] /. full}, {x, -1, 1, .05}], 
   GridLines -> Automatic, Frame -> True, 
   PlotLegends -> Placed[tag /. full, Right], 
   FrameLabel -> {x, Subscript[\[Gamma], 2]}];
plotGamma2Moments = 
  Plot[Evaluate[
    gamma2[x] /. # & /@ {dsm0, gn, dsm1, dsm2, dsm3}], {x, -1, 1}, 
   PlotLegends -> 
    Placed[tag /. # & /@ {dsm0, gn, dsm1, dsm2, dsm3}, Right]];
plotGamma2 = Show[plotGamma2Full, plotGamma2Moments];

pRange = Automatic;

plotKappa2Full = 
  ListPlot[Table[{x, kappa2[x] /. full}, {x, -1, 1, .05}], 
   GridLines -> Automatic, Frame -> True, 
   PlotLegends -> Placed[tag /. full, Right], 
   FrameLabel -> {x, Subscript[\[Kappa], 2]}];
plotKappa2Moments = 
  Plot[Evaluate[
    kappa2[x] /. # & /@ {dsm0, gn, dsm1, dsm2, dsm3}], {x, -1, 1}, 
   PlotLegends -> 
    Placed[tag /. # & /@ {dsm0, gn, dsm1, dsm2, dsm3}, Right]];
plotKappa2 = Show[plotKappa2Full, plotKappa2Moments];


(********************* Plot third order moments ***************************************************)

pRange = Automatic;

plotAlpha3Full = 
  ListPlot[Table[{x, alpha3[x] /. full}, {x, -1, 1, .05}], 
   GridLines -> Automatic, Frame -> True, 
   PlotLegends -> Placed[tag /. full, Right], 
   FrameLabel -> {x, Subscript[\[Gamma], 3]}];
plotAlpha3Moments = 
  Plot[Evaluate[ 
    alpha3[x] /. # & /@ {dsm0, gn, dsm1, dsm2, dsm3}], {x, -1, 1}, 
   PlotLegends -> 
    Placed[tag /. # & /@ {dsm0, gn, dsm1, dsm2, dsm3}, Right]];
plotAlpha3 = Show[plotAlpha3Full, plotAlpha3Moments];

pRange = Automatic;

plotGamma3Full = 
  ListPlot[Table[{x, gamma3[x] /. full}, {x, -1, 1, .05}], 
   GridLines -> Automatic, Frame -> True, 
   PlotLegends -> Placed[tag /. full, Right], 
   FrameLabel -> {x, Subscript[\[Gamma], 3]}];
plotGamma3Moments = 
  Plot[Evaluate[ 
    gamma3[x] /. # & /@ {dsm0, gn, dsm1, dsm2, dsm3}], {x, -1, 1}, 
   PlotLegends -> 
    Placed[tag /. # & /@ {dsm0, gn, dsm1, dsm2, dsm3}, Right]];
plotGamma3 = Show[plotGamma3Full, plotGamma3Moments];

pRange = Automatic;

plotKappa3Full = 
  ListPlot[Table[{x, kappa3[x] /. full}, {x, -1, 1, .05}], 
   GridLines -> Automatic, Frame -> True, 
   PlotLegends -> Placed[tag /. full, Right], 
   FrameLabel -> {x, Subscript[\[Kappa], 3]}];
plotKappa3Moments = 
  Plot[Evaluate[
    kappa3[x] /. # & /@ {dsm0, gn, dsm1, dsm2, dsm3}], {x, -1, 1}, 
   PlotLegends -> 
    Placed[tag /. # & /@ {dsm0, gn, dsm1, dsm2, dsm3}, Right]];
plotKappa3 = Show[plotKappa3Full, plotKappa3Moments];

(********************* Save Plots *****************************************************************)

(* set output directory *)
identifier = "H" <> StringReplace[ToString[H], "." -> "p"] <> "Fr" <> 
   StringReplace[ToString[Fr], "." -> "p"];
   
outputDir=FileNameJoin[{ParentDirectory[currentDirectory],"plots","compare",identifier}];

Quiet[CreateDirectory[outputDir]];

destinationFiles=FileNameJoin[{outputDir,#}]& /@{
"legend.pdf",
"h.pdf",
"um.pdf",
"wm.pdf",
"qm.pdf",
"alpha1.pdf",
"gamma1.pdf",
"kappa1.pdf",
"alpha2.pdf",
"gamma2.pdf",
"kappa2.pdf",
"alpha3.pdf",
"gamma3.pdf",
"kappa3.pdf"
};

MapThread[Export,{destinationFiles, 
{
legend,
plotH,
plotUm,
plotWm,
plotQm,
plotAlpha1,
plotGamma1,
plotKappa1,
plotAlpha2,
plotGamma2,
plotKappa2,
plotAlpha3,
plotGamma3,
plotKappa3
}
}];

Print["Plots saved to " <> outputDir <> ". Terminating..."]

