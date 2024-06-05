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
        "State of art"=>"art.md",
        "Examples"=>"smooth.md",
        "Extensions"=>"nonsmooth.md",
        "API"=>"api.md"
    ]
)

deploydocs(;
    repo="github.com/Xiaomingzzhang/InvariantManifolds.jl",
    devbranch="master",tag_prefix="v0.2.1"
)
