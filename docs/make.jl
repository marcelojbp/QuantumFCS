# docs/make.jl  — classic Documenter deploy to gh-pages

import Pkg
Pkg.activate(@__DIR__)    
Pkg.instantiate()

using Documenter
using QuantumFCS

DocMeta.setdocmeta!(QuantumFCS, :DocTestSetup, :(using QuantumFCS); recursive=true)

makedocs(
    sitename = "QuantumFCS.jl",
    modules  = [QuantumFCS],
    format   = Documenter.HTML(),
    doctest  = true,
    pages    = [
        "Home"       => "index.md",
        "Quickstart" => "quickstart.md",
        "API"        => "api.md",
    ],
)

deploydocs(
    repo      = "github.com/marcelojbp/QuantumFCS",
    devbranch = "main",
)
