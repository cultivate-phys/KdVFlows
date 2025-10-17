# KdVFlows

A Mathematica package for computing symmetries and flows of the KdV equation. For mathematical background, see:

```bibtex
@book{miwa1991solitons,
    author = "Miwa, Tetsuji and Date, Etsuro and Jimbo, Michio",
    title = "{Solitons: Differential Equations, Symmetries and Infinite Dimensional Algebras}",
    isbn = "978-1-107-40419-9",
    series = "Cambridge Tracts in Mathematics"
}
```

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
See the [examples](examples/) directory for examples.

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
- **`SetWeighting[wU, wDx]`** - Sets weight of u and ∂x (default: 2, 1)

## Requirements

- Wolfram Mathematica (tested with version 14.1)

## License

MIT License - see [LICENSE](LICENSE) file for details.