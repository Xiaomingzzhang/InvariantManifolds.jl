using InvariantManifolds
using Documenter

DocMeta.setdocmeta!(InvariantManifolds, :DocTestSetup, :(using InvariantManifolds); recursive=true)

makedocs(;
    modules=[InvariantManifolds],
    authors="XiaomingZhang",
    sitename="InvariantManifolds.jl",
    warnonly = true,
    format=Documenter.HTML(;
        canonical="https://Xiaomingzzhang.github.io/InvariantManifolds.jl",
        edit_link="master",
        assets=String[],
        prettyurls = get(ENV, "CI", nothing) == "true",
        collapselevel = 1
    ),
    pages = [
        "Home"=>"index.md",
        "Getting Started: One-Dimensional Smooth Manifolds"=>"smooth_one.md",
        "Two-Dimensional Smooth Manifolds"=>"smooth_two.md",
        "One-Dimensional Non-Smooth Manifolds"=>"non_smooth_one.md",
        "Two-dimensional Non-Smooth Manifolds"=>"non_smooth_two.md",
        "API"=>"api.md"
    ]
)

deploydocs(;
    repo="github.com/Xiaomingzzhang/InvariantManifolds.jl",
    devbranch = "master"
)
