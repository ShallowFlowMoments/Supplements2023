(* ::Package:: *)

BeginPackage["DSM2`"]


DSM2::usage="DSM2[H_,Fr_,\[Eta]_,\[Sigma]_,n_] gives the steady flow simulation with DSM-2 model";


Begin["Private`"]


DSM2[H_,Fr_,\[Eta]_,\[Sigma]_,n_]:=(

(* create multipurpose identifier tag *)
tag="DSM-2;H="<>ToString[H]<>";Fr="<>ToString[Fr];

(* set mesh size and create grid *)
\[CapitalDelta]=2/(n+1);
xs[i_]=-1+i \[CapitalDelta];

(* define bottom topography *)
hb[x_]=\[Eta]/H E^(-x^2/(2\[Sigma]^2));dxhb[x_]=-((E^(-(x^2/(2 \[Sigma]^2))) x \[Eta]/H)/\[Sigma]^2);

(* define all the unknowns *)
varsFull=Flatten[Table[{h[i],u[i],w[i],q[i],\[Alpha]1[i],\[Gamma]1[i],\[Kappa]1[i],\[Alpha]2[i],\[Gamma]2[i],\[Kappa]2[i]},{i,0,n+1}]];

(* set up horizontal boundary conditions with finite differences *)
bcx={
(*left*)
h[0],
u[0]-(3 u[1]-3 u[2]+u[3]),
\[Alpha]1[0]-(3 \[Alpha]1[1]-3 \[Alpha]1[2]+\[Alpha]1[3]),
\[Alpha]2[0]-(3 \[Alpha]2[1]-3 \[Alpha]2[2]+\[Alpha]2[3]),
w[0],
\[Gamma]1[0],
\[Gamma]2[0],
q[0],
\[Kappa]1[0],
\[Kappa]2[0],
(*right*)
h[n+1]-(h[-2+n]+3 (-h[-1+n]+h[n])),
u[n+1],
\[Alpha]1[n+1],
\[Alpha]2[n+1],
w[n+1]-(w[-2+n]+3 (-w[-1+n]+w[n])),
\[Gamma]1[n+1]-(\[Gamma]1[-2+n]+3 (-\[Gamma]1[-1+n]+\[Gamma]1[n])),
\[Gamma]2[n+1]-(\[Gamma]2[-2+n]+3 (-\[Gamma]2[-1+n]+\[Gamma]2[n])),
q[n+1]-(q[-2+n]+3 (-q[-1+n]+q[n])),
\[Kappa]1[n+1]-(\[Kappa]1[-2+n]+3 (-\[Kappa]1[-1+n]+\[Kappa]1[n])),
\[Kappa]2[n+1]-(\[Kappa]2[-2+n]+3 (-\[Kappa]2[-1+n]+\[Kappa]2[n]))
};

(* set up linear equation system with finite differences *)
eqns=Join[Flatten[
Table[{
(h[i+1]-h[i-1])/(2\[CapitalDelta])+(u[i+1]-u[i-1])/(2\[CapitalDelta]), (*incompress. condition*)

Fr^2 (u[i+1]-u[i-1])/(2\[CapitalDelta])+(h[i+1]-h [i-1])/(2\[CapitalDelta])+(q[i+1]-q[i-1])/(2\[CapitalDelta])+dxhb[xs[i]], (*u momentum balance*)
Fr^2 (\[Alpha]1[i+1]-\[Alpha]1[i-1])/(2\[CapitalDelta])+(\[Kappa]1[i+1]-\[Kappa]1[i-1])/(2\[CapitalDelta]),
Fr^2 (\[Alpha]2[i+1]-\[Alpha]2[i-1])/(2\[CapitalDelta])+(\[Kappa]2[i+1]-\[Kappa]2[i-1])/(2\[CapitalDelta]),
H^2 Fr^2 (w[i+1]-w[i-1])/(2\[CapitalDelta])-(q[i]+\[Kappa]1[i]+\[Kappa]2[i]),(*w momentum balance*)
H^2 Fr^2 (\[Gamma]1[i+1]-\[Gamma]1[i-1])/(2\[CapitalDelta])-3(-q[i]+\[Kappa]1[i]+\[Kappa]2[i]),
H^2 Fr^2 (\[Gamma]2[i+1]-\[Gamma]2[i-1])/(2\[CapitalDelta])-5(q[i]-\[Kappa]1[i]+\[Kappa]2[i]),
(u[i+1]-u[i-1])/(2\[CapitalDelta])-(dxhb[xs[i]]-w[i]+\[Gamma]1[i]-\[Gamma]2[i]),(*pressure*)
(\[Alpha]1[i+1]-\[Alpha]1[i-1])/(2\[CapitalDelta])-3(dxhb[xs[i]]-w[i]-\[Gamma]1[i]+\[Gamma]2[i]),
(\[Alpha]2[i+1]-\[Alpha]2[i-1])/(2\[CapitalDelta])-5(dxhb[xs[i]]-w[i]-\[Gamma]1[i]-\[Gamma]2[i])
},{i,1,n}
]],
Flatten[bcx]
];

(* switch to matrix-vector formulation and solve system *)
{b,A}=CoefficientArrays[eqns,varsFull];
sol=LinearSolve[A,-b];

(* extract point values from solution *)
listH=Transpose[{Table[xs[i],{i,0,n+1}],sol[[1;; ;; 10]]+1}];
listUm=Transpose[{Table[xs[i],{i,0,n+1}],sol[[2;; ;; 10]]+1}];
listWm=Transpose[{Table[xs[i],{i,0,n+1}],sol[[3;; ;; 10]]}];
listQm=Transpose[{Table[xs[i],{i,0,n+1}],sol[[4;; ;; 10]]}];
list\[Alpha]1=Transpose[{Table[xs[i],{i,0,n+1}],sol[[5;; ;; 10]]}];
list\[Gamma]1=Transpose[{Table[xs[i],{i,0,n+1}],sol[[6;; ;; 10]]}];
list\[Kappa]1=Transpose[{Table[xs[i],{i,0,n+1}],sol[[7;; ;; 10]]}];
list\[Alpha]2=Transpose[{Table[xs[i],{i,0,n+1}],sol[[8;; ;; 10]]}];
list\[Gamma]2=Transpose[{Table[xs[i],{i,0,n+1}],sol[[9;; ;; 10]]}];
list\[Kappa]2=Transpose[{Table[xs[i],{i,0,n+1}],sol[[10;; ;; 10]]}];

(* interpolate the point values *)
hN[tag]=Interpolation[listH,InterpolationOrder->1];
umN[tag]=Interpolation[listUm,InterpolationOrder->1];
wmN[tag]=Interpolation[listWm,InterpolationOrder->1];
qmN[tag]=Interpolation[listQm,InterpolationOrder->1];
alpha1N[tag]=Interpolation[list\[Alpha]1,InterpolationOrder->1];
gamma1N[tag]=Interpolation[list\[Gamma]1,InterpolationOrder->1];
kappa1N[tag]=Interpolation[list\[Kappa]1,InterpolationOrder->1];
alpha2N[tag]=Interpolation[list\[Alpha]2,InterpolationOrder->1];
gamma2N[tag]=Interpolation[list\[Gamma]2,InterpolationOrder->1];
kappa2N[tag]=Interpolation[list\[Kappa]2,InterpolationOrder->1];

(* return the variables *)
{
Symbol["h"]->hN[tag],
Symbol["um"]->umN[tag],Symbol["alpha1"]->alpha1N[tag],Symbol["alpha2"]->alpha2N[tag],
Symbol["wm"]->wmN[tag],Symbol["gamma1"]->gamma1N[tag],Symbol["gamma2"]->gamma2N[tag],
Symbol["qm"]->qmN[tag],Symbol["kappa1"]->kappa1N[tag],Symbol["kappa2"]->kappa2N[tag],
Symbol["tag"]->tag
}

)


End[];


EndPackage[];
