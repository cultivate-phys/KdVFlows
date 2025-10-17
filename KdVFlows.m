BeginPackage["KdVFlows`"];

KdV::usage = "KdV[u_] := u u_x + u_{xxx}. Use as KdV[u[x]].";
FrechetD::usage = "FrechetD[F_, h_] computes the algebraic Fréchet derivative dF[u].h = Sum_k (\[PartialD]F/\[PartialD]u^{(k)}) \[PartialD]_x^k h.";
CommuteQ::usage = "CommuteQ[F_] returns True if FrechetApply[KdV[u],F] - FrechetApply[F,KdV[u]] == 0.";
BasisKdV::usage = "BasisKdV[w, maxOrd] returns a list of homogeneous differential monomials of weight w with highest derivative order \[LessEqual] maxOrd.";
SymmetrySolve::usage = "SymmetrySolve[w, maxOrd] solves [KdV, F]=0 within the homogeneous basis of weight w and returns a normalized flow F.";
LenardFlow::usage = "LenardFlow[n] returns the (2n+1)-th KdV flow by a Lenard-style recursion (assumes decaying or zero-mean boundary).";
FlowsUpTo::usage = "FlowsUpTo[n] returns {F1,F3,...,F_{2n+1}}.";
SetWeighting::usage = "SetWeighting[wU, wDx] sets wt(u)=wU and wt(\[PartialD]x)=wDx (default 2,1).";
SetBoundary::usage = "SetBoundary[mode] sets boundary handling for Dinv: \"Decaying\" (default) or \"PeriodicZeroMean\".";

u::usage = "Dependent variable u[x]";
x::usage = "Independent variable x";

(* ---------- 0. Weighting rules ---------- *)
wU = 2; wDx = 1;

SetWeighting[wUin_:2, wDXin_:1] := (wU = wUin; wDx = wDXin);

(* ---------- 1. Core KdV and algebraic Fréchet derivative ---------- *)
ClearAll[KdV, FrechetD, MaxDerivOrder];

KdV[uu_] := uu*D[uu, x] + D[uu, {x, 3}];

MaxDerivOrder[expr_] := Max@Join[{0}, Cases[expr, Derivative[k_][u][x] :> k, All]];

FrechetD[F_, h_] := Module[{K = MaxDerivOrder[F]},
  Sum[D[F, Derivative[k][u][x]]*D[h, {x, k}], {k, 0, K}] // Expand
];

(* commutator on expressions in x *)
ClearAll[CommExpr, CommuteQ];
CommExpr[F_] := Expand[FrechetD[KdV[u[x]], F] - FrechetD[F, KdV[u[x]]]];
CommuteQ[F_] := Simplify[CommExpr[F] === 0];

(* ---------- 2. Basis generation (homogeneous differential polynomials) ---------- *)

ClearAll[MonomialFromData];
MonomialFromData[m0_Integer, ks_List] := Module[{tally, term, i},
  tally = Tally[ks];
  term = If[m0 > 0, u[x]^m0, 1];
  Do[
    With[{k = tally[[i, 1]], mult = tally[[i, 2]]},
      term = term * Power[Derivative[k][u][x], mult]
    ],
    {i, Length[tally]}
  ];
  term
];

ClearAll[BasisKdV];
BasisKdV[w_Integer, maxOrd_Integer] := Module[
  {res = {}, m0, target, parts, ks},
  Do[
    target = w - m0*wU;
    If[target < 0, Continue[]];
    parts = IntegerPartitions[target, All, Range[2, 2 + maxOrd]];
    Do[
      ks = parts[[i]] - 2; (* convert parts to derivative orders *)
      If[Max[If[ks === {}, {0}, ks]] <= maxOrd,
        AppendTo[res, MonomialFromData[m0, ks]]
      ]
      , {i, Length[parts]}
    ]
    , {m0, 0, Floor[w/wU]}
  ];
  res = DeleteDuplicates[res // Expand];
  res = SortBy[res, {
      -Total[(Cases[#, Derivative[k_][u][x] :> k, All] /. {} -> {0})] &,
      -Exponent[#, u[x]] &
    }
  ];
  res
];

(* ---------- 3. Coefficient solving via commutator = 0 ---------- *)

(* Replace jet variables by commuting symbols U0, U1, ... *)
ClearAll[toJetSymbols, jetVars];
toJetSymbols[expr_, maxOrd_Integer] := Module[{rep},
  rep = Join[{u[x] -> Symbol["U0"]},
              Table[Derivative[k][u][x] -> Symbol["U" <> ToString[k]], {k, 1, maxOrd}]];
  expr /. rep
];
jetVars[maxOrd_Integer] := Array[Symbol["U" <> ToString[#]] &, maxOrd + 1, 0];

ClearAll[SymmetrySolve];
Options[SymmetrySolve] = {NormalizeLeading -> True};
SymmetrySolve[w_Integer, maxOrd_Integer:100, OptionsPattern[]] := Module[
  {B, n, ci, F, CE, ord, CEc, U, rules, exprs, mat, vec, sol},
  B = BasisKdV[w, maxOrd];
  n = Length[B];
  If[n == 0, Return[$Failed]];
  ci = Array[c, n];
  F = Plus @@ (ci*B) // Expand;
  CE = CommExpr[F] // Expand;
  ord = MaxDerivOrder[CE];
  CEc = toJetSymbols[CE, ord];
  U = jetVars[ord];
  rules = CoefficientRules[CEc, U];
  exprs = Values[rules]; (* linear in ci *)
  mat = Table[Coefficient[exprs[[i]], ci[[j]]], {i, Length[exprs]}, {j, n}] // Normal;
  vec = - (exprs /. Thread[ci -> 0]) // Normal;
  If[TrueQ[OptionValue[NormalizeLeading]],
    mat = Join[mat, {UnitVector[n, 1]}];
    vec = Join[vec, {1}];
  ];
  sol = Quiet@LinearSolve[mat, vec];
  If[Head[sol] =!= List && Head[sol] =!= Vector, Return[$Failed]];
  (F /. Thread[ci -> sol]) // Expand
];

EndPackage[];