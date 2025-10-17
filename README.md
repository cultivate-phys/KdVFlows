# KdVFlows

A Mathematica package for computing symmetries and flows of the KdV equation.

## Overview

KdVFlows provides tools for symbolic computation with the KdV equation and its hierarchy of commuting flows. The package implements:

- KdV equation: $u_t = u u_x + u_xxx$
- Basis generation for homogeneous differential polynomials
- Fréchet derivatives
- Automated solving for symmetries at any weight
- Symmetry analysis via commutator equations

## Installation

1. Download or clone this repository
2. Place `KdVFlows.m` in a location accessible to Mathematica
3. Load the package:

```mathematica
Get["path/to/KdVFlows.m"]
```

Or set your notebook directory and use:

```mathematica
SetDirectory@NotebookDirectory[];
Get["KdVFlows.m"];
```

## Quick Start

```mathematica
(* Load the package *)
Get["KdVFlows.m"];

(*Get basis elements of weight 5 with max derivative order 10*)
BasisKdV[5, 10]
(*Returns: {(u^(3))[x],u[x] (u^\[Prime])[x]}*)

(*Compute Fréchet derivative*)
FrechetD[u[x]^2*u'[x], h]
(*Returns: 2 h u[x] (u^\[Prime])[x]*)

(*Compute a symmetry of weight 5*)
F5 = SymmetrySolve[5]
(*Returns: u[x] (u^\[Prime])[x]+(u^(3))[x]*)

(*Check if it commutes with KdV*)
CommuteQ[F5]
(*Returns: True*)
```

## Main Functions

- **`KdV[u]`** - Computes the KdV equation: $u u_x + u_xxx$
- **`BasisKdV[w, maxOrd]`** - Generates basis of differential monomials of weight `w` with maximum derivative order `maxOrd`
- **`FrechetD[F, h]`** - Computes the Fréchet derivative of `F` in direction `h`
- **`SymmetrySolve[w]`** - Solves for a symmetry of weight `w`
- **`CommuteQ[F]`** - Tests if `F` commutes with the KdV flow