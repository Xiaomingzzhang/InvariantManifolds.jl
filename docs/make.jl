using InvariantManifolds
using Documenter

DocMeta.setdocmeta!(InvariantManifolds, :DocTestSetup, :(using InvariantManifolds); recursive=true)

makedocs(;
    modules=[InvariantManifolds],
    authors="XiaomingZhang",
    sitename="InvariantManifolds.jl",
    format=Documenter.HTML(;
        canonical="https://Xiaomingzzhang.github.io/InvariantManifolds.jl",
        edit_link="master",
        assets=String[],
    ),
    pages = [
        "Home"=>"index.md",
        "State of art"=>"state.md",
        "Smooth mapping"=>"smooth.md",
        "Time-T-map of non-smooth ODE"=>"nonsmooth.md",
        "API"=>"api.md"
    ]
)

deploydocs(;
    repo="github.com/Xiaomingzzhang/InvariantManifolds.jl",
    devbranch="gh-pages",
)
