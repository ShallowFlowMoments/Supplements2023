(* ::Package:: *)

BeginPackage["DSM0`"]


DSM0::usage="DSM0[H_,Fr_,\[Eta]_,\[Sigma]_,n_] gives the steady flow simulation with DSM-0 model";


Begin["Private`"]


DSM0[H_,Fr_,\[Eta]_,\[Sigma]_,n_]:=(

(* create multipurpose identifier tag *)
tag="DSM-0;H="<>ToString[H]<>";Fr="<>ToString[Fr];

(* set mesh size and create grid *)
\[CapitalDelta]=2/(n+1);
xs[i_]=-1+i \[CapitalDelta];

(* define bottom topography *)
hb[x_]=\[Eta]/H E^(-x^2/(2\[Sigma]^2));dxhb[x_]=-((E^(-(x^2/(2 \[Sigma]^2))) x \[Eta]/H)/\[Sigma]^2);

(* define all the unknowns *)
varsFull=Flatten[Table[{h[i],u[i],w[i],q[i]},{i,0,n+1}]];

(* set up horizontal boundary conditions with finite differences *)
bcx={
(*left*)
h[0],
u[0]-(3 u[1]-3 u[2]+u[3]),
w[0],
q[0],
(*right*)
h[n+1]-(h[-2+n]+3 (-h[-1+n]+h[n])),
u[n+1],
w[n+1]-(w[-2+n]+3 (-w[-1+n]+w[n])),
q[n+1]-(q[-2+n]+3 (-q[-1+n]+q[n]))
};

(* set up linear equation system with finite differences *)
eqns=Join[Flatten[
Table[{
(h[i+1]-h[i-1])/(2\[CapitalDelta])+(u[i+1]-u[i-1])/(2\[CapitalDelta]), (*incompress. condition*)

Fr^2 (u[i+1]-u[i-1])/(2\[CapitalDelta])+(h[i+1]-h [i-1])/(2\[CapitalDelta])+(q[i+1]-q[i-1])/(2\[CapitalDelta])+dxhb[xs[i]], (*u momentum balance*)
H^2 Fr^2 (w[i+1]-w[i-1])/(2\[CapitalDelta])-q[i],(*w momentum balance*)
w[i]-dxhb[xs[i]]+(u[i+1]-u[i-1])/(2\[CapitalDelta])(*pressure*)
},{i,1,n}
]],
Flatten[bcx]
];

(* switch to matrix-vector formulation and solve system *)
{b,A}=CoefficientArrays[eqns,varsFull];
sol2=LinearSolve[A,-b];

(* extract point values from solution *)
listH=Transpose[{Table[xs[i],{i,0,n+1}],sol2[[1;; ;; 4]]+1}];
listUm=Transpose[{Table[xs[i],{i,0,n+1}],sol2[[2;; ;; 4]]+1}];
listWm=Transpose[{Table[xs[i],{i,0,n+1}],sol2[[3;; ;; 4]]}];
listQm=Transpose[{Table[xs[i],{i,0,n+1}],sol2[[4;; ;; 4]]}];

(* interpolate the point values *)
hN[tag]=Interpolation[listH,InterpolationOrder->1];
umN[tag]=Interpolation[listUm,InterpolationOrder->1];
wmN[tag]=Interpolation[listWm,InterpolationOrder->1];
qmN[tag]=Interpolation[listQm,InterpolationOrder->1];

(* return the variables *)
{
Symbol["h"]->hN[tag],
Symbol["um"]->umN[tag],
Symbol["wm"]->wmN[tag],
Symbol["qm"]->qmN[tag],
Symbol["tag"]->tag
}

)


End[];


EndPackage[];
