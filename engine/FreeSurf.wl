(* ::Package:: *)

BeginPackage["FreeSurf`"]


FreeSurf::usage="FreeSurf[H_,Fr_,\[Eta]_,\[Sigma]_,n_] gives the steady flow simulation for free surface flow, vertically resolved";


Begin["`Private`"];


FreeSurf[H_,Fr_,\[Eta]_,\[Sigma]_,n_]:=(

(* create multipurpose identifier tag *)
tagFull="FreeSurf;H="<>ToString[H]<>";Fr="<>ToString[Fr];

(* mesh size *)
\[CapitalDelta]x=2/(n+1);\[CapitalDelta]\[Zeta]=1/(n+1);

(* create grid axes *)
xs[i_]=-1+i \[CapitalDelta]x;
\[Zeta]s[j_]=j \[CapitalDelta]\[Zeta];

(* define bottom topography *)
hb[x_]=\[Eta] E^(-x^2/(2\[Sigma]^2));dxhb[x_]=-((E^(-(x^2/(2 \[Sigma]^2))) x \[Eta])/\[Sigma]^2);

(* define all unknowns *)
varsFull=Join[Table[h[i],{i,0,n+1}],Flatten[Table[{u[i,j],w[i,j],q[i,j]},{i,0,n+1},{j,0,n+1}]]];

(* set up horizontal boundary conditions with finite differences *)
bcx=Join[{h[0],h[n+1]-(h[-2+n]+3 (-h[-1+n]+h[n]))},
Table[{
(*left*)
{
u[0,j]-3 u[1,j]+3 u[2,j]-u[3,j], (*interpolating polynomial*)
w[0,j],
q[0,j]
},
(*right*)
{
u[n+1,j],
w[n+1,j]-3w[n,j]+3 w[n-1,j]-w[n-2,j],
q[n+1,j]-3q[n,j]+3q[n-1,j]-q[n-2,j]
}
},{j,0,n+1}]];

(* set up vertical boundary conditions with finite differences *)
bc\[Zeta]=Table[{
(*bottom*)
{
u[i,0]-3 u[i,1]+3 u[i,2]-u[i,3],
w[i,0]-  1/H dxhb[xs[i]],
q[i,0]-3 q[i,1]+3 q[i,2]-q[i,3]
},
(*top*)
{
u[i,n+1]-3u[i,n]+3 u[i,n-1]-u[i,n-2],
w[i,n+1]-3w[i,n]+3w[i,n-1]-w[i,n-2],
q[i,n+1]
}
},{i,1,n}];

(* set up linear system with finite differences *)
eqns=Join[
Table[(h[i+1]-h[i-1])/(2\[CapitalDelta]x)+\[CapitalDelta]\[Zeta]/(2\[CapitalDelta]x) ((u[i+1,0]-u[i-1,0])/2 +\!\(
\*UnderoverscriptBox[\(\[Sum]\), \(j = 1\), \(n\)]\((u[i + 1, j] - u[i - 1, j])\)\)+(u[i+1,n+1]-u[i-1,n+1])/2), {i,1,n}],(*incompress. condition*)
Flatten[
Table[
{
Fr^2 (u[i+1,j]-u[i-1,j])/(2\[CapitalDelta]x)+(h[i+1]-h[i-1])/(2\[CapitalDelta]x)+(q[i+1,j]-q[i-1,j])/(2\[CapitalDelta]x)+1/H dxhb[xs[i]], (*u*)
H^2 Fr^2 (w[i+1,j]-w[i-1,j])/(2\[CapitalDelta]x)+(q[i,j+1]-q[i,j-1])/(2\[CapitalDelta]\[Zeta]),(*w*) 
(w[i,j+1]-w[i,j-1])/(2\[CapitalDelta]\[Zeta])+(u[i+1,j]-u[i-1,j])/(2\[CapitalDelta]x)(*q*)
}
,{i,1,n},{j,1,n}]
],
Flatten[bcx],
Flatten[bc\[Zeta]]
];

(* switch to matrix-vector formulation and solve system *)
{b,A}=CoefficientArrays[eqns,varsFull];
sol=LinearSolve[A,-b];

(*define the legendre polynomials *)
poly[i_,var_]:=LegendreP[i,-2var+1];

(* extract point values for the height and the other variables*)
heights=sol[[1;;n+2]];
rest=sol[[n+3;;]];

(* define helper function for the creation of a look up table of the values for later interpolation*)
valueList[vals_]:=Transpose[{
Flatten[Table[xs[i],{i,0,n+1},{j,0,n+1}]],
Flatten[Table[\[Zeta]s[j],{i,0,n+1},{j,0,n+1}]],
Flatten[vals]}]; 

(* define helper function for the creation of a table with the polynomial values on the grid*)
repeatedPolyVals[degree_]:=Table[poly[degree,\[Zeta]s[Mod[j,n+2]]],{j,0,(n+2)^2-1}];

(* create look up table for height *)
listH[tagFull]=Transpose[{Table[xs[i],{i,0,n+1}],heights+1}];

(* create look up table for horizontal velocity variables um and alpha 1-3 *)
(* shift by one because of the linearization around um=1 *)
{listU,listUp1,listUp2,listUp3}=Table[valueList[rest[[1;;;;3]]*repeatedPolyVals[i]],{i,0,3}];
listU[[All,3]]=listU[[All,3]]+1;

(* create look up table for vertical velocity variables wm and gamma 1-3 *)
{listW,listWp1,listWp2,listWp3}=Table[valueList[rest[[2;;;;3]]*repeatedPolyVals[i]],{i,0,3}];

(* create look up table for nonhydrostatic pressure variables qm and kappa 1-3*)
{listQ,listQp1,listQp2,listQp3}=Table[valueList[rest[[3;;;;3]]*repeatedPolyVals[i]],{i,0,3}];

(* Create the interpolation functions for h, u, w and q *)
hN[tagFull]=Interpolation[listH[tagFull],InterpolationOrder->1];

uN[tagFull]=Interpolation[listU,InterpolationOrder->1];
up1N[tagFull]=Interpolation[listUp1,InterpolationOrder->1];
up2N[tagFull]=Interpolation[listUp2,InterpolationOrder->1];
up3N[tagFull]=Interpolation[listUp3,InterpolationOrder->1];

wN[tagFull]=Interpolation[listW,InterpolationOrder->1];
wp1N[tagFull]=Interpolation[listWp1,InterpolationOrder->1];
wp2N[tagFull]=Interpolation[listWp2,InterpolationOrder->1];
wp3N[tagFull]=Interpolation[listWp3,InterpolationOrder->1];

qN[tagFull]=Interpolation[listQ,InterpolationOrder->1];
qp1N[tagFull]=Interpolation[listQp1,InterpolationOrder->1];
qp2N[tagFull]=Interpolation[listQp2,InterpolationOrder->1];
qp3N[tagFull]=Interpolation[listQp3,InterpolationOrder->1];

(* determine means from vertically resolved functions *)
umN[tagFull][x_]=Integrate[uN[tagFull][x,\[Zeta]],{\[Zeta],0,1}];
wmN[tagFull][x_]=Integrate[wN[tagFull][x,\[Zeta]],{\[Zeta],0,1}];
qmN[tagFull][x_]=Integrate[qN[tagFull][x,\[Zeta]],{\[Zeta],0,1}];

(* determine momentes from the vertically resolved functions. Scale with 3,5,7 because the Legendre polynomials are not normalized *)
alpha1N[tagFull][x_]=3Integrate[up1N[tagFull][x,\[Zeta]],{\[Zeta],0,1}];
gamma1N[tagFull][x_]=3Integrate[wp1N[tagFull][x,\[Zeta]],{\[Zeta],0,1}];
kappa1N[tagFull][x_]=3Integrate[qp1N[tagFull][x,\[Zeta]],{\[Zeta],0,1}];

alpha2N[tagFull][x_]=5Integrate[up2N[tagFull][x,\[Zeta]],{\[Zeta],0,1}];
gamma2N[tagFull][x_]=5Integrate[wp2N[tagFull][x,\[Zeta]],{\[Zeta],0,1}];
kappa2N[tagFull][x_]=5Integrate[qp2N[tagFull][x,\[Zeta]],{\[Zeta],0,1}];

alpha3N[tagFull][x_]=7Integrate[up3N[tagFull][x,\[Zeta]],{\[Zeta],0,1}];
gamma3N[tagFull][x_]=7Integrate[wp3N[tagFull][x,\[Zeta]],{\[Zeta],0,1}];
kappa3N[tagFull][x_]=7Integrate[qp3N[tagFull][x,\[Zeta]],{\[Zeta],0,1}];

(* define the z to zeta transformation *)
\[Zeta][tagFull][x_,z_]=(z-1/H hb[x])/hN[tagFull][x];
z[tagFull][x_,\[Zeta]_]=\[Zeta] hN[tagFull][x]+1/H hb[x];

(* reshift the variables to the original coordiate system with z axis *)
uM[tagFull][x_,z_]:=If[0<=\[Zeta][tagFull][x,z]<=1,uN[tagFull][x,\[Zeta][tagFull][x,z]],Undefined];
wM[tagFull][x_,z_]:=If[0<=\[Zeta][tagFull][x,z]<=1,wN[tagFull][x,\[Zeta][tagFull][x,z]],Undefined];
qM[tagFull][x_,z_]:=If[0<=\[Zeta][tagFull][x,z]<=1,qN[tagFull][x,\[Zeta][tagFull][x,z]],Undefined];

(* return the variables *)
{
Symbol["h"]->hN[tagFull],
Symbol["utilde"]->uM[tagFull],Symbol["u"]->uN[tagFull],Symbol["um"]->umN[tagFull],Symbol["alpha1"]->alpha1N[tagFull],Symbol["alpha2"]->alpha2N[tagFull],Symbol["alpha3"]->alpha3N[tagFull],
Symbol["wtilde"]->wM[tagFull],Symbol["w"]->wN[tagFull],Symbol["wm"]->wmN[tagFull],Symbol["gamma1"]->gamma1N[tagFull],Symbol["gamma2"]->gamma2N[tagFull],Symbol["gamma3"]->gamma3N[tagFull],
Symbol["qtilde"]->qM[tagFull],Symbol["q"]->qN[tagFull],Symbol["qm"]->qmN[tagFull],Symbol["kappa1"]->kappa1N[tagFull],Symbol["kappa2"]->kappa2N[tagFull],Symbol["kappa3"]->kappa3N[tagFull],
Symbol["tag"]->tagFull
}

)


End[];


EndPackage[];
