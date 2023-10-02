(* ::Package:: *)

BeginPackage["GN`"]


GN::usage="GN[H_,Fr_,\[Eta]_,\[Sigma]_,n_] gives the steady flow simulation with Green-Naghdi model";


Begin["Private`"]


GN[H_,Fr_,\[Eta]_,\[Sigma]_,n_]:=(

(* create multipurpose identifier tag *)
tag="Green-Naghdi;H="<>ToString[H]<>";Fr="<>ToString[Fr];

(* set mesh size and create grid *)
\[CapitalDelta]=2/(n+1);
xs[i_]=-1+i \[CapitalDelta];

(* define bottom topography *)
hb[x_]=\[Eta]/H E^(-x^2/(2\[Sigma]^2));dxhb[x_]=-((E^(-(x^2/(2 \[Sigma]^2))) x \[Eta]/H)/\[Sigma]^2);

(* define all the unknowns *)
varsFull=Flatten[Table[{h[i],u[i],w[i],q[i],\[Gamma]1[i],\[Kappa]1[i]},{i,0,n+1}]];

(* set up horizontal boundary conditions with finite differences *)
bcx={
(*left*)
h[0],
u[0]-(3 u[1]-3 u[2]+u[3]),
w[0],
\[Gamma]1[0],
q[0],
\[Kappa]1[0],
(*right*)
h[n+1]-(h[-2+n]+3 (-h[-1+n]+h[n])),
u[n+1],
w[n+1]-(w[-2+n]+3 (-w[-1+n]+w[n])),
\[Gamma]1[n+1]-(\[Gamma]1[-2+n]+3 (-\[Gamma]1[-1+n]+\[Gamma]1[n])),
q[n+1]-(q[-2+n]+3 (-q[-1+n]+q[n])),
\[Kappa]1[n+1]-(\[Kappa]1[-2+n]+3 (-\[Kappa]1[-1+n]+\[Kappa]1[n]))
};

(* set up linear equation system with finite differences *)
eqns=Join[Flatten[
Table[{
(h[i+1]-h[i-1])/(2\[CapitalDelta])+(u[i+1]-u[i-1])/(2\[CapitalDelta]), (*incompress. condition*)

Fr^2 (u[i+1]-u[i-1])/(2\[CapitalDelta])+(h[i+1]-h [i-1])/(2\[CapitalDelta])+(q[i+1]-q[i-1])/(2\[CapitalDelta])+dxhb[xs[i]], (*u momentum balance*)
H^2 Fr^2 (w[i+1]-w[i-1])/(2\[CapitalDelta])-2\[Kappa]1[i],(*w momentum balance*)
H^2 Fr^2 (\[Gamma]1[i+1]-\[Gamma]1[i-1])/(2\[CapitalDelta])-6(\[Kappa]1[i]-q[i]),
(u[i+1]-u[i-1])/(2\[CapitalDelta])-2\[Gamma]1[i](*pressure*),
\[Gamma]1[i]-dxhb[xs[i]]+w[i]
},{i,1,n}
]],
Flatten[bcx]
];

(* switch to matrix-vector formulation and solve system *)
{b,A}=CoefficientArrays[eqns,varsFull];
sol=LinearSolve[A,-b];

(* extract point values from solution *)
listH=Transpose[{Table[xs[i],{i,0,n+1}],sol[[1;; ;; 6]]+1}];
listUm=Transpose[{Table[xs[i],{i,0,n+1}],sol[[2;; ;; 6]]+1}];
listWm=Transpose[{Table[xs[i],{i,0,n+1}],sol[[3;; ;; 6]]}];
listQm=Transpose[{Table[xs[i],{i,0,n+1}],sol[[4;; ;; 6]]}];
list\[Gamma]1=Transpose[{Table[xs[i],{i,0,n+1}],sol[[5;; ;; 6]]}];
list\[Kappa]1=Transpose[{Table[xs[i],{i,0,n+1}],sol[[6;; ;; 6]]}];

(* interpolate the point values *)
hN[tag]=Interpolation[listH,InterpolationOrder->1];
umN[tag]=Interpolation[listUm,InterpolationOrder->1];
wmN[tag]=Interpolation[listWm,InterpolationOrder->1];
qmN[tag]=Interpolation[listQm,InterpolationOrder->1];
gamma1N[tag]=Interpolation[list\[Gamma]1,InterpolationOrder->1];
kappa1N[tag]=Interpolation[list\[Kappa]1,InterpolationOrder->1];

(* return the variables *)
{
Symbol["h"]->hN[tag],
Symbol["um"]->umN[tag],
Symbol["wm"]->wmN[tag],Symbol["gamma1"]->gamma1N[tag],
Symbol["qm"]->qmN[tag],Symbol["kappa1"]->kappa1N[tag],
Symbol["tag"]->tag
}

)


End[];


EndPackage[];
