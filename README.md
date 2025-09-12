# QuantumFCS.jl

[![CI](https://github.com/marcelojbp/QuantumFCS/actions/workflows/CI.yml/badge.svg)](https://github.com/marcelojbp/QuantumFCS/actions/workflows/CI.yml)
[![Docs: dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://marcelojbp.github.io/QuantumFCS/dev/)
[![License: MIT](https://img.shields.io/badge/license-MIT-green.svg)](LICENSE)
[![Cite](https://img.shields.io/badge/cite-CITATION.bib-informational.svg)](CITATION.bib)

*A Julia package for computing Full Counting Statistics (FCS) of open quantum systems.*

## Description

**QuantumFCS.jl** provides tools to study **Full Counting Statistics** (FCS) of quantum transport and quantum optics models based on Lindblad master equations. 
It implements a recursive method in which the $n+1$-th cumulant is computed using the $n$-th and the application of the Drazin inverse of the Liouvillean.

- The package follows approach introduced in [Flindt et al., Phys. Rev. B 82, 155407 (2010)](https://arxiv.org/abs/1002.4506).  
- It is designed for efficient numerical calculations using sparse or dense linear algebra.

## Features

- Compute current cumulants to arbitrary order
- Support for **sparse** and **dense** matrix backends  
- Methods for integration with `QuantumOptics.jl` operators, states.
- Does not need `QuantumOptics.jl`; you can build your Liouvillean however you like and use `QuantumFCS.jl`.


### Main function
- `fcscumulants_recursive(L, mJ, nC, rho_ss; ...)` – compute FCS cumulants up to order `nC`, taking as input 
a (vectorised) Liouvillean, `L`, the steady-state `rho_ss`, and a vector of monitored jump `mJ`. An optional argument, `nu` can be passed to attribute sign 
and unitful weighs to the monitored jumps
- `fcscumulants_recursive(H,J, mJ, nC, rho_ss; ...)` – Method for quantum optics---you can pass, instead of the Liouvellean,
the Hamiltonian and the vector of all jump operators.

*(See the [API docs](https://marcelojbp.github.io/QuantumFCS/dev/) for the full list.)*

## Installation

Until registration, install directly from GitHub:

```julia
] add https://github.com/marcelojbp/QuantumFCS
```
### Quickstart
```julia
using QuantumFCS

L = #Vectorised Liouvillean
rho_ss = #Steady-state (0-th eigenvalue of the L)

# compute first 3 cumulants
avg, var, skw = fcscumulants_recursive(L, J, 3, rho_ss)
```



<!-- [![Build Status](https://github.com/aarondanielphysics/QuantumFCS.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/aarondanielphysics/QuantumFCS.jl/actions/workflows/CI.yml?query=branch%3Amain) -->
